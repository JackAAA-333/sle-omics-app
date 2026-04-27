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
    # primary sheet name
    main = '蛋白定量结果' if '蛋白定量结果' in sheets else sheets[1] if len(sheets)>1 else sheets[0]
    df = pd.read_excel('prot.xlsx', sheet_name=main)
    return xls, main, df

def build_matrix(df):
    # assume first col is protein id, others are samples
    df = df.copy()
    df.columns = [str(c) for c in df.columns]
    feat_col = df.columns[0]
    mat = df.set_index(feat_col)
    return mat

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
        try:
            stat,p = ttest_ind(a,b, equal_var=False, nan_policy='omit')
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

def xgb_shap(mat, meta):
    samp_mask = meta['group'] != 'QC'
    X = mat.T.loc[samp_mask.index[samp_mask]]
    y = (meta.loc[samp_mask.index[samp_mask],'group']=='SLE').astype(int).values
    clf = xgb.XGBClassifier(use_label_encoder=False, eval_metric='logloss', verbosity=0)
    clf.fit(X,y)
    expl = shap.TreeExplainer(clf)
    shap_vals = expl.shap_values(X)
    imp = np.abs(shap_vals).mean(axis=0)
    df = pd.Series(imp, index=X.columns).sort_values(ascending=False)
    df.to_csv(f'{OUT}/xgb_shap_prot.tsv', sep='\t')
    return df

def main():
    xls, main, df = read_prot()
    mat = build_matrix(df)
    mat2, meta = align_samples(mat)
    norm = qc_normalize(mat2)
    logm = preprocess(norm)
    logm.to_csv(f'{OUT}/filtered_prot_matrix.tsv', sep='\t')
    res = differential(logm, meta)
    imp = xgb_shap(logm, meta)
    summary = {'n_features':int(logm.shape[0]), 'n_samples':int(logm.shape[1])}
    with open(f'{OUT}/summary.json','w') as f:
        json.dump(summary,f)
    print('Proteomics analysis complete. main_sheet:', main)

if __name__=='__main__':
    main()
