import os
import pandas as pd
import numpy as np
import xgboost as xgb
import shap

OUT = 'outputs_multiomics'
os.makedirs(OUT, exist_ok=True)

def load():
    met = pd.read_csv('outputs/filtered_matrix.tsv', sep='\t', index_col=0)
    prot = pd.read_csv('outputs_prot/filtered_prot_matrix.tsv', sep='\t', index_col=0)
    # transpose to samples x features
    metT = met.T
    protT = prot.T
    # align samples
    common = metT.index.intersection(protT.index)
    metT = metT.loc[common]
    protT = protT.loc[common]
    return metT, protT

def concat_and_model(metT, protT):
    X = pd.concat([metT, protT], axis=1)
    # load meta
    meta = pd.read_csv('outputs/sample_metadata.csv', index_col=0)
    meta = meta.loc[X.index]
    y = (meta['group']=='SLE').astype(int).values
    clf = xgb.XGBClassifier(use_label_encoder=False, eval_metric='logloss', verbosity=0)
    clf.fit(X, y)
    expl = shap.TreeExplainer(clf)
    shap_vals = expl.shap_values(X)
    imp = pd.Series(np.abs(shap_vals).mean(axis=0), index=X.columns).sort_values(ascending=False)
    imp.to_csv(f'{OUT}/xgb_shap_multi.tsv', sep='\t')
    # correlations between top features
    top_met = imp[imp.index.str.startswith('M')].head(30).index.intersection(X.columns)
    top_prot = imp[~imp.index.str.startswith('M')].head(30).index.intersection(X.columns)
    corr = X[top_met].corrwith(X[top_prot[0]], axis=0) if len(top_prot)>0 else pd.Series()
    # compute full pairwise Spearman
    import scipy.stats as stats
    rows = []
    for m in top_met:
        for p in top_prot:
            rho,pv = stats.spearmanr(X[m], X[p], nan_policy='omit')
            rows.append({'met':m,'prot':p,'rho':rho,'p':pv})
    pd.DataFrame(rows).to_csv(f'{OUT}/met_prot_spearman.tsv', sep='\t', index=False)
    return imp

def main():
    metT, protT = load()
    imp = concat_and_model(metT, protT)
    print('Multi-omics integration complete. Top features saved to', OUT)

if __name__=='__main__':
    main()
