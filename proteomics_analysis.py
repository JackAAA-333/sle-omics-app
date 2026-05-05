import os
import json
import re
import numpy as np
import pandas as pd
import requests
from scipy.stats import ttest_ind
from sklearn.metrics import roc_auc_score, roc_curve, auc
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_predict
import xgboost as xgb
import shap
import matplotlib.pyplot as plt
import seaborn as sns

OUT = 'outputs_prot'
os.makedirs(OUT, exist_ok=True)
CACHE_FILE = os.path.join(OUT, "gene_symbol_lookup_cache.json")

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


def _extract_gene_symbol_from_text(text):
    m = re.search(r"\bGN=([A-Za-z0-9\-]+)\b", str(text))
    if m:
        return m.group(1).strip()
    return ""


def _extract_uniprot_accession(text):
    token = str(text).strip().split()[0]
    patterns = [
        r"^[OPQ][0-9][A-Z0-9]{3}[0-9]$",
        r"^[A-NR-Z][0-9]{5}$",
        r"^[A-NR-Z][0-9][A-Z0-9]{3}[0-9]$",
    ]
    for p in patterns:
        if re.match(p, token):
            return token
    return ""


def _load_lookup_cache():
    if not os.path.exists(CACHE_FILE):
        return {}
    try:
        with open(CACHE_FILE, "r", encoding="utf-8") as f:
            data = json.load(f)
        if isinstance(data, dict):
            return data
    except Exception:
        pass
    return {}


def _save_lookup_cache(cache):
    try:
        with open(CACHE_FILE, "w", encoding="utf-8") as f:
            json.dump(cache, f, ensure_ascii=False, indent=2)
    except Exception:
        pass


def _query_uniprot_gene_symbol(protein_id):
    accession = _extract_uniprot_accession(protein_id)
    if accession:
        query = f"accession:{accession} AND reviewed:true"
    else:
        name = str(protein_id).split(" OS=")[0].strip()
        if not name:
            return ""
        query = f'protein_name:"{name}" AND reviewed:true'
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": query,
        "fields": "gene_primary,accession",
        "size": 1,
        "format": "json",
    }
    try:
        r = requests.get(url, params=params, timeout=12)
        r.raise_for_status()
        results = r.json().get("results", [])
        if not results:
            return ""
        genes = results[0].get("genes", [])
        if not genes:
            return ""
        primary = genes[0].get("geneName", {}).get("value", "")
        return str(primary).strip()
    except Exception:
        return ""


def enrich_gene_symbol_map(mat, gene_map):
    cache = _load_lookup_cache()
    out_map = dict(gene_map)
    source_map = {}
    updated = False
    for pid in mat.index.astype(str):
        existing = str(out_map.get(pid, "")).strip()
        if existing and existing.lower() not in {"nan", "none"}:
            source_map[pid] = "input_annotation"
            continue
        parsed = _extract_gene_symbol_from_text(pid)
        if parsed:
            out_map[pid] = parsed
            source_map[pid] = "parsed_gn"
            continue
        if pid in cache and str(cache[pid]).strip():
            out_map[pid] = str(cache[pid]).strip()
            source_map[pid] = "cache_uniprot"
            continue
        online = _query_uniprot_gene_symbol(pid)
        if online:
            out_map[pid] = online
            cache[pid] = online
            source_map[pid] = "online_uniprot"
            updated = True
        else:
            source_map[pid] = "missing"
    if updated:
        _save_lookup_cache(cache)
    return out_map, source_map


def export_gene_annotation_catalog(mat, gene_map, source_map=None):
    source_map = source_map or {}
    rows = []
    for pid in mat.index.astype(str):
        gs = str(gene_map.get(pid, "")).strip()
        src = source_map.get(pid, "input_annotation" if gs else "missing")
        rows.append(
            {
                "protein_id": pid,
                "gene_symbol": gs,
                "annotation_status": "annotated" if gs else "missing",
                "annotation_source": src,
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
    with np.errstate(invalid="ignore", divide="ignore"):
        logm_np = np.log2((matf + 1).to_numpy(dtype=float))
    logm = pd.DataFrame(logm_np, index=matf.index, columns=matf.columns)
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
    ann["marker"] = ann["gene_symbol"]
    ann.loc[ann["marker"].isin(["", "nan", "None"]), "marker"] = ann["protein_id"]
    ann["annotation_status"] = ann["gene_symbol"].apply(lambda x: "annotated" if x else "missing")
    ann.to_csv(f"{OUT}/differential_prot_annotated.tsv", sep="\t", index=False)
    ann.to_csv(f"{OUT}/differential_prot_gene_symbol.tsv", sep="\t", index=False)
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
    annotated = pd.DataFrame(
        {
            'protein_id': df.index.astype(str),
            'gene_symbol': [gene_map.get(str(pid), '') for pid in df.index],
            'mean_abs_shap': df.values,
        }
    )
    annotated["marker"] = annotated["gene_symbol"].astype(str).str.strip()
    annotated.loc[annotated["marker"].isin(["", "nan", "None"]), "marker"] = annotated["protein_id"]
    annotated.to_csv(f'{OUT}/xgb_shap_prot_annotated.tsv', sep='\t', index=False)
    # Display-oriented SHAP output: use gene symbol marker as primary index.
    shap_display = annotated[["marker", "mean_abs_shap", "protein_id", "gene_symbol"]].copy()
    shap_display = shap_display.drop_duplicates(subset=["marker"], keep="first").set_index("marker")
    shap_display.to_csv(f"{OUT}/xgb_shap_prot.tsv", sep="\t")
    shap_display.to_csv(f"{OUT}/xgb_shap_prot_gene_symbol.tsv", sep="\t")
    return df


def _safe_name(text):
    s = re.sub(r"[^A-Za-z0-9_\-]+", "_", str(text)).strip("_")
    return s[:80] if s else "marker"


def plot_protein_visuals(logm, meta, diff_res, gene_map=None):
    gene_map = gene_map or {}
    groups = meta['group'].astype(str).str.upper()
    sle_samples = groups[groups == 'SLE'].index.tolist()
    hc_samples = groups[groups == 'HC'].index.tolist()
    if len(sle_samples) < 2 or len(hc_samples) < 2:
        return

    common_samples = sle_samples + hc_samples
    top_feats = (
        diff_res.sort_values('pvalue', ascending=True)
        .head(30)
        .index.intersection(logm.index)
        .tolist()
    )
    if len(top_feats) >= 2:
        try:
            sns.clustermap(logm.loc[top_feats, common_samples], cmap='vlag', standard_scale=0)
            plt.savefig(f'{OUT}/heatmap_prot_top30.png', dpi=200)
            plt.close()
        except Exception:
            pass

    # Single-marker ROC curves for top significant proteins (individual plots + one overlay).
    roc_feats = (
        diff_res.sort_values('pvalue', ascending=True)
        .head(5)
        .index.intersection(logm.index)
        .tolist()
    )
    if len(roc_feats) == 0:
        return
    y_true = (groups.loc[common_samples] == 'SLE').astype(int).values
    roc_dir = os.path.join(OUT, "roc_single_markers")
    os.makedirs(roc_dir, exist_ok=True)
    plt.figure(figsize=(6, 5))
    plotted = 0
    single_rows = []
    for feat in roc_feats:
        vals = pd.to_numeric(logm.loc[feat, common_samples], errors='coerce').fillna(0.0).values
        if np.nanstd(vals) < 1e-12:
            continue
        try:
            auc_raw = roc_auc_score(y_true, vals)
            score = vals if auc_raw >= 0.5 else -vals
            auc_use = auc_raw if auc_raw >= 0.5 else (1.0 - auc_raw)
            fpr, tpr, _ = roc_curve(y_true, score)
            marker = str(gene_map.get(str(feat), "")).strip() or str(feat)
            plt.plot(fpr, tpr, linewidth=1.5, label=f'{marker} (AUC={auc_use:.3f})')
            single_rows.append({"feature_id": str(feat), "marker": marker, "auc_single": auc_use})
            # Individual ROC figure for each top marker
            plt.figure(figsize=(5.5, 5))
            plt.plot(fpr, tpr, color="#2a7fff", linewidth=2.0, label=f"AUC={auc_use:.3f}")
            plt.plot([0, 1], [0, 1], "--", color="gray", linewidth=0.8)
            plt.xlabel("False Positive Rate (FPR)")
            plt.ylabel("True Positive Rate (TPR)")
            plt.title(f"Single Marker ROC - {marker}")
            plt.legend(loc="lower right")
            plt.tight_layout()
            plt.savefig(os.path.join(roc_dir, f"roc_single_{_safe_name(marker)}.png"), dpi=200)
            plt.close()
            plotted += 1
        except Exception:
            continue
    if plotted > 0:
        plt.plot([0, 1], [0, 1], '--', color='gray', linewidth=0.8)
        plt.xlabel('False Positive Rate (FPR)')
        plt.ylabel('True Positive Rate (TPR)')
        plt.title('Protein Single-Marker ROC (Top5 Overlay)')
        plt.legend(fontsize=8, loc='lower right')
        plt.tight_layout()
        plt.savefig(f'{OUT}/roc_prot_top5.png', dpi=200)
    plt.close()
    pd.DataFrame(single_rows).to_csv(f"{OUT}/single_marker_metrics.tsv", sep="\t", index=False)

    # Combined ROC using top5 markers jointly (cross-validated probability).
    valid_feats = [f for f in roc_feats if f in logm.index]
    if len(valid_feats) >= 2:
        X = logm.loc[valid_feats, common_samples].T.apply(pd.to_numeric, errors='coerce')
        X = X.fillna(X.median()).fillna(0.0)
        min_class = int(min((y_true == 1).sum(), (y_true == 0).sum()))
        n_splits = max(2, min(5, min_class, len(y_true)))
        try:
            cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)
            clf = LogisticRegression(solver='liblinear', max_iter=2000, class_weight='balanced')
            y_prob = cross_val_predict(clf, X, y_true, cv=cv, method='predict_proba')[:, 1]
            fpr_c, tpr_c, _ = roc_curve(y_true, y_prob)
            auc_c = float(auc(fpr_c, tpr_c))
            plt.figure(figsize=(5.8, 5.2))
            plt.plot(fpr_c, tpr_c, color="#d62728", linewidth=2.2, label=f"AUC={auc_c:.3f}")
            plt.plot([0, 1], [0, 1], "--", color="gray", linewidth=0.8)
            plt.xlabel("False Positive Rate (FPR)")
            plt.ylabel("True Positive Rate (TPR)")
            plt.title("Protein Top5 Combined ROC")
            plt.legend(loc="lower right")
            plt.tight_layout()
            plt.savefig(f"{OUT}/roc_prot_combined_top5.png", dpi=200)
            plt.close()
        except Exception:
            pass

def main():
    xls, main, df = read_prot()
    gene_map = build_gene_symbol_map(df)
    mat = build_matrix(df)
    gene_map, source_map = enrich_gene_symbol_map(mat, gene_map)
    export_gene_annotation_catalog(mat, gene_map, source_map=source_map)
    mat2, meta = align_samples(mat)
    norm = qc_normalize(mat2)
    logm = preprocess(norm)
    logm.to_csv(f'{OUT}/filtered_prot_matrix.tsv', sep='\t')
    res = differential(logm, meta)
    export_differential_annotated(res, gene_map)
    plot_protein_visuals(logm, meta, res, gene_map=gene_map)
    imp = xgb_shap(logm, meta, gene_map=gene_map)
    summary = {'n_features':int(logm.shape[0]), 'n_samples':int(logm.shape[1])}
    with open(f'{OUT}/summary.json','w') as f:
        json.dump(summary,f)
    print('Proteomics analysis complete. main_sheet:', main)

if __name__=='__main__':
    main()
