import os
import re
import json
import time
import requests
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import hypergeom
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import RepeatedStratifiedKFold, cross_val_score
from sklearn.metrics import roc_auc_score
import xgboost as xgb
import shap

from preprocess_and_analyze import (
    load_metab, detect_samples, build_matrix, build_sample_metadata,
    qc_normalize, transform_and_scale
)

OUT = 'outputs_advanced'
os.makedirs(OUT, exist_ok=True)


def multiple_testing(df):
    # assume df has 'pvalue'
    from statsmodels.stats.multitest import multipletests
    p = df['pvalue'].fillna(1).values
    _, bh, _, _ = multipletests(p, method='fdr_bh')
    _, bonf, _, _ = multipletests(p, method='bonferroni')
    df['fdr_bh'] = bh
    df['p_bonf'] = bonf
    return df


def build_kegg_pathway_map(kegg_ids):
    # kegg_ids: iterable of KEGG compound IDs like C00031 or may include 'cpd:C00031'
    mapping = {}
    session = requests.Session()
    for kid in set([k for k in kegg_ids if pd.notna(k) and k!='']):
        kid_clean = str(kid)
        kid_clean = kid_clean.replace('cpd:','').strip()
        try:
            url = f'https://rest.kegg.jp/link/pathway/cpd:{kid_clean}'
            r = session.get(url, timeout=10)
            if r.status_code == 200 and r.text.strip():
                lines = r.text.strip().split('\n')
                pathways = [ln.split('\t')[1].split(':')[1] for ln in lines if '\t' in ln]
                mapping[kid_clean] = pathways
            else:
                mapping[kid_clean] = []
        except Exception:
            mapping[kid_clean] = []
        time.sleep(0.1)
    return mapping


def pathway_enrichment(sig_kegg, background_kegg, mapping):
    # mapping: kegg_id -> [pathway1, pathway2]
    # build pathway -> set(compounds)
    pw2cpd = {}
    for c,pws in mapping.items():
        for pw in pws:
            pw2cpd.setdefault(pw,set()).add(c)
    M = len(background_kegg)
    N = len(sig_kegg)
    results = []
    bg_set = set([b.replace('cpd:','') for b in background_kegg if pd.notna(b) and b!=''])
    sig_set = set([s.replace('cpd:','') for s in sig_kegg if pd.notna(s) and s!=''])
    for pw, cpds in pw2cpd.items():
        k = len(cpds & sig_set)
        K = len(cpds & bg_set)
        if K==0:
            continue
        # hypergeometric: P(X>=k)
        pval = hypergeom.sf(k-1, M, K, N)
        results.append((pw, k, K, pval))
    resdf = pd.DataFrame(results, columns=['pathway','k_sig','K_bg','pvalue']).sort_values('pvalue')
    if not resdf.empty:
        from statsmodels.stats.multitest import multipletests
        resdf['fdr_bh'] = multipletests(resdf['pvalue'].fillna(1), method='fdr_bh')[1]
    return resdf


def run_xgboost_shap(X, y, topk=20):
    rkf = RepeatedStratifiedKFold(n_splits=5, n_repeats=3, random_state=42)
    model = xgb.XGBClassifier(n_estimators=200, max_depth=4, use_label_encoder=False, eval_metric='auc', random_state=42)
    aucs = cross_val_score(model, X, y, cv=rkf, scoring='roc_auc', n_jobs=1)
    model.fit(X, y)
    explainer = shap.Explainer(model)
    shap_values = explainer(X)
    # mean absolute shap per feature
    mean_abs = np.abs(shap_values.values).mean(axis=0)
    feat_imp = pd.Series(mean_abs, index=X.columns).sort_values(ascending=False)
    top = feat_imp.head(topk)
    # save shap summary plot
    try:
        shap.plots.beeswarm(shap_values[:, :topk], show=False)
        import matplotlib.pyplot as plt
        plt.tight_layout()
        plt.savefig(os.path.join(OUT,'shap_beeswarm_top{}.png'.format(topk)), dpi=200)
        plt.close()
    except Exception:
        pass
    return aucs, feat_imp


def pubmed_search(keywords, retmax=5):
    # use NCBI E-utilities esearch to get PMIDs
    base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi'
    params = {'db':'pubmed','retmode':'json','retmax':retmax,'term':keywords}
    try:
        r = requests.get(base, params=params, timeout=10)
        j = r.json()
        ids = j.get('esearchresult',{}).get('idlist',[])
        return ids
    except Exception:
        return []


def main():
    df = load_metab('metab.xlsx')
    meta_cols, sample_cols = detect_samples(df)
    sample_meta = build_sample_metadata(sample_cols)
    mat = build_matrix(df, sample_cols)
    mat_norm = qc_normalize(mat, sample_meta)
    mat_log, _ = transform_and_scale(mat_norm)

    # strict filter >20% already applied earlier; reapply here
    na_frac = mat_log.isna().mean(axis=1)
    zero_frac = (mat_log==0).mean(axis=1)
    keep = (na_frac <= 0.2) & (zero_frac <= 0.2)
    mat_filt = mat_log.loc[keep]

    # differential test
    sle = sample_meta[sample_meta['group']=='SLE'].index.tolist()
    hc = sample_meta[sample_meta['group']=='HC'].index.tolist()
    res = []
    for feat,row in mat_filt.iterrows():
        x = row[sle].dropna().values
        y = row[hc].dropna().values
        if len(x)<3 or len(y)<3:
            p=np.nan
        else:
            _,p = stats.ttest_ind(x,y,equal_var=False)
        res.append((feat,p))
    resdf = pd.DataFrame(res, columns=['feature','pvalue']).set_index('feature')
    resdf = multiple_testing(resdf)
    resdf.to_csv(os.path.join(OUT,'differential_results_advanced.tsv'), sep='\t')

    # get KEGG IDs mapping from original df
    id_map = {}
    if 'KEGG ID' in df.columns:
        for i,row in df.iterrows():
            key = str(row.get('Met ID', i))
            id_map[key] = row.get('KEGG ID', '')

    # select significant features
    sig = resdf[resdf['fdr_bh']<0.05].index.tolist()
    sig_kegg = [id_map.get(s,'') for s in sig]
    bg_kegg = [id_map.get(s,'') for s in mat_filt.index.tolist()]

    mapping = build_kegg_pathway_map(bg_kegg)
    enrich = pathway_enrichment(sig_kegg, bg_kegg, mapping)
    enrich.to_csv(os.path.join(OUT,'kegg_enrichment.tsv'), sep='\t')

    # non-linear XGBoost + SHAP
    common = sample_meta[sample_meta['group'].isin(['SLE','HC'])].index.tolist()
    X = mat_filt[common].T.fillna(0)
    y = sample_meta.loc[common,'group'].map({'SLE':1,'HC':0}).astype(int)
    aucs, feat_imp = run_xgboost_shap(X, y, topk=30)
    pd.Series(aucs, name='xgb_cv_auc').to_csv(os.path.join(OUT,'xgb_cv_aucs.tsv'), sep='\t')
    feat_imp.to_csv(os.path.join(OUT,'xgb_shap_feature_importance.tsv'), sep='\t')

    # PubMed search for top 10 features (by XGB SHAP)
    top10 = feat_imp.head(10).index.tolist()
    pub_summary = {}
    for feat in top10:
        # try metabolite name from original df
        name = None
        row = df[df['Met ID']==feat]
        if not row.empty and 'Metabolite Name' in row.columns:
            name = row.iloc[0]['Metabolite Name']
        if not name:
            name = feat
        term = f"{name} AND systemic lupus erythematosus"
        ids = pubmed_search(term, retmax=5)
        pub_summary[feat] = {'name': name, 'pmids': ids, 'count': len(ids)}
        time.sleep(0.34)

    with open(os.path.join(OUT,'pubmed_summary.json'),'w',encoding='utf-8') as fh:
        json.dump(pub_summary, fh, ensure_ascii=False, indent=2)

    # write conclusions
    conclusions = {
        'n_sig': int((resdf['fdr_bh']<0.05).sum()),
        'top_xgb_auc_mean': float(np.mean(aucs)),
        'notes': 'See kegg_enrichment.tsv, xgb_shap_feature_importance.tsv, pubmed_summary.json for details.'
    }
    with open(os.path.join(OUT,'conclusions.json'),'w',encoding='utf-8') as fh:
        json.dump(conclusions, fh, ensure_ascii=False, indent=2)

    print('Advanced analysis complete. Outputs in', OUT)


if __name__ == '__main__':
    main()
