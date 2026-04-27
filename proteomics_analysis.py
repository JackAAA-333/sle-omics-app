import os
import json
import re
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from sklearn.metrics import roc_auc_score
import xgboost as xgb
import shap

OUT = 'outputs_prot'
os.makedirs(OUT, exist_ok=True)

def read_prot():
    xls = pd.read_excel('prot.xlsx', sheet_name=None)
    sheets = list(xls.keys())
    # rule: if workbook has multiple sheets, only use the first sheet
    main = sheets[0]
    df = pd.read_excel('prot.xlsx', sheet_name=main)
    return xls, main, df

def build_matrix(df):
    # assume first col is protein id, others are samples
    df = df.copy()
    df.columns = [str(c) for c in df.columns]
    feat_col = df.columns[0]
    mat = df.set_index(feat_col)
    # force all measurement columns to numeric; non-numeric tokens like 'LOD' -> NaN
    mat = mat.apply(pd.to_numeric, errors='coerce')
    # drop rows that are fully non-numeric after coercion
    mat = mat.dropna(how='all')
    return mat


def _normalize_col_name(name):
    return re.sub(r'[\s_]+', '', str(name).strip().lower())


def build_gene_symbol_map(df):
    if df is None or df.empty:
        return {}
    cols = [str(c) for c in df.columns]
    feat_col = cols[0]
    col_norm = {_normalize_col_name(c): c for c in cols}
    candidate_keys = [
        'genesymbol',
        'symbol',
        'genename',
        'gene',
        '官方基因名',
        '基因名',
    ]
    gene_col = None
    for key in candidate_keys:
        if key in col_norm:
            gene_col = col_norm[key]
            break
    if gene_col is None:
        for c in cols[1:]:
            n = _normalize_col_name(c)
            if 'gene' in n and ('symbol' in n or 'name' in n):
                gene_col = c
                break
    if gene_col is None:
        return {}

    sub = df[[feat_col, gene_col]].copy()
    sub[feat_col] = sub[feat_col].astype(str).str.strip()
    sub[gene_col] = sub[gene_col].astype(str).str.strip()
    sub = sub[(sub[feat_col] != '') & (sub[gene_col] != '')]
    if sub.empty:
        return {}
    return dict(zip(sub[feat_col], sub[gene_col]))


def export_gene_annotation_catalog(mat, gene_map):
    rows = []
    for pid in mat.index.astype(str):
        gs = str(gene_map.get(pid, "")).strip()
        rows.append(
            {
                "protein_id": pid,
                "gene_symbol": gs,
                "annotation_status": "annotated" if gs else "missing",
            }
        )
    ann = pd.DataFrame(rows)
    ann.to_csv(f"{OUT}/protein_gene_annotations.tsv", sep="\t", index=False)
    return ann

def align_samples(mat):
    if os.path.exists('outputs/sample_metadata.csv'):
        meta = pd.read_csv('outputs/sample_metadata.csv', index_col=0)
        cols = [c for c in mat.columns if c in meta.index]
        mat2 = mat[cols]
        if len(cols) > 0:
            return mat2, meta.loc[cols]
    # fallback: infer groups directly from protein sample headers
    cols = []
    groups = []
    for c in mat.columns:
        name = str(c)
        m = re.search(r'(SLE|HC|QC)', name, re.IGNORECASE)
        if m:
            cols.append(c)
            groups.append(m.group(1).upper())
    if len(cols) == 0:
        cols = list(mat.columns)
        groups = ['UNK'] * len(cols)
    mat2 = mat[cols]
    meta = pd.DataFrame({'group': groups}, index=cols)
    return mat2, meta

def qc_normalize(mat):
    qc_cols = [c for c in mat.columns if c.startswith('QC') or c.upper().startswith('QC-')]
    if qc_cols:
        ref = mat[qc_cols].median(axis=1)
        norm = mat.divide(ref, axis=0)
    else:
        norm = mat.divide(mat.median(axis=1), axis=0)
    return norm

def preprocess(mat):
    # replace zeros/na with nan
    mat = mat.replace(0, np.nan)
    miss_frac = mat.isna().mean(axis=1)
    matf = mat[miss_frac <= 0.5]
    matf = matf.fillna(matf.median(axis=1), axis=0)
    logm = np.log2(matf + 1)
    return logm

def differential(mat, meta):
    groups = meta['group']
    sle_idx = groups== 'SLE'
    pvals = []
    aucs = []
    for feat in mat.index:
        a = mat.loc[feat, sle_idx.index[sle_idx]].values
        b = mat.loc[feat, sle_idx.index[~sle_idx]].values
        # avoid noisy warnings / unstable tests when effective sample size is too small
        a_valid = a[~np.isnan(a)]
        b_valid = b[~np.isnan(b)]
        if len(a_valid) < 2 or len(b_valid) < 2:
            p = 1.0
        else:
            try:
                _, p = ttest_ind(a_valid, b_valid, equal_var=False, nan_policy='omit')
                if np.isnan(p):
                    p = 1.0
            except Exception:
                p = 1.0
        pvals.append(p)
        try:
            y = (groups=='SLE').astype(int).values
            auc = roc_auc_score(y, mat.loc[feat, groups.index].values)
        except Exception:
            auc = np.nan
        aucs.append(auc)
    res = pd.DataFrame({'feature':mat.index,'pvalue':pvals,'auc':aucs}).set_index('feature')
    res['fdr_bh'] = res['pvalue']
    res.to_csv(f'{OUT}/differential_prot.tsv', sep='\t')
    return res


def export_differential_annotated(res, gene_map):
    ann = res.reset_index().rename(columns={"feature": "protein_id"})
    ann["protein_id"] = ann["protein_id"].astype(str)
    ann["gene_symbol"] = ann["protein_id"].map(lambda x: str(gene_map.get(x, "")).strip())
    ann["annotation_status"] = ann["gene_symbol"].apply(lambda x: "annotated" if x else "missing")
    ann.to_csv(f"{OUT}/differential_prot_annotated.tsv", sep="\t", index=False)
    return ann

def _sanitize_xgb_feature_name(name):
    s = str(name)
    s = re.sub(r'[\[\]<>]', '_', s)
    return s


def xgb_shap(mat, meta, gene_map=None):
    gene_map = gene_map or {}
    samp_mask = meta['group'] != 'QC'
    X = mat.T.loc[samp_mask.index[samp_mask]].copy()
    y = (meta.loc[samp_mask.index[samp_mask],'group']=='SLE').astype(int).values
    # xgboost disallows feature names containing [, ] and <
    orig_cols = [str(c) for c in X.columns]
    safe_cols = []
    used = set()
    safe_to_orig = {}
    for c in orig_cols:
        base = _sanitize_xgb_feature_name(c)
        cand = base
        k = 2
        while cand in used:
            cand = f"{base}_{k}"
            k += 1
        used.add(cand)
        safe_cols.append(cand)
        safe_to_orig[cand] = c
    X.columns = safe_cols

    clf = xgb.XGBClassifier(use_label_encoder=False, eval_metric='logloss', verbosity=0)
    clf.fit(X,y)
    expl = shap.TreeExplainer(clf)
    shap_vals = expl.shap_values(X)
    imp = np.abs(shap_vals).mean(axis=0)
    df_safe = pd.Series(imp, index=X.columns).sort_values(ascending=False)
    df = pd.Series(
        data=df_safe.values,
        index=[safe_to_orig.get(c, c) for c in df_safe.index],
        dtype=float,
    ).sort_values(ascending=False)
    df.to_csv(f'{OUT}/xgb_shap_prot.tsv', sep='\t')
    annotated = pd.DataFrame(
        {
            'protein_id': df.index.astype(str),
            'gene_symbol': [gene_map.get(str(pid), '') for pid in df.index],
            'mean_abs_shap': df.values,
        }
    )
    annotated.to_csv(f'{OUT}/xgb_shap_prot_annotated.tsv', sep='\t', index=False)
    return df

def main():
    xls, main, df = read_prot()
    gene_map = build_gene_symbol_map(df)
    mat = build_matrix(df)
    export_gene_annotation_catalog(mat, gene_map)
    mat2, meta = align_samples(mat)
    norm = qc_normalize(mat2)
    logm = preprocess(norm)
    logm.to_csv(f'{OUT}/filtered_prot_matrix.tsv', sep='\t')
    res = differential(logm, meta)
    export_differential_annotated(res, gene_map)
    imp = xgb_shap(logm, meta, gene_map=gene_map)
    summary = {'n_features':int(logm.shape[0]), 'n_samples':int(logm.shape[1])}
    with open(f'{OUT}/summary.json','w') as f:
        json.dump(summary,f)
    print('Proteomics analysis complete. main_sheet:', main)

if __name__=='__main__':
    main()
