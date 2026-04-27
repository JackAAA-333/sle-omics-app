import os
import re
import json
import warnings
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegressionCV, LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import RepeatedStratifiedKFold, cross_val_score
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.exceptions import ConvergenceWarning
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import seaborn as sns


OUTDIR = 'outputs'
os.makedirs(OUTDIR, exist_ok=True)
warnings.filterwarnings(
    "ignore",
    message="invalid value encountered in cast",
    category=RuntimeWarning,
)


def load_metab(path='metab.xlsx'):
    df = pd.read_excel(path)
    return df


def load_prot(path='prot.xlsx', sheet_name='蛋白定量结果'):
    # list sheets and try to read the main quantification sheet
    try:
        xl = pd.ExcelFile(path)
        sheets = xl.sheet_names
    except Exception:
        return None, []
    if sheet_name in sheets:
        df = pd.read_excel(path, sheet_name=sheet_name)
    else:
        # fallback to first sheet
        df = pd.read_excel(path, sheet_name=0)
    return df, sheets


def detect_samples(df):
    pattern = re.compile(r'^(SLE|HC|QC)[-_]?\d+', re.IGNORECASE)
    sample_cols = [c for c in df.columns if pattern.match(str(c))]
    meta_cols = [c for c in df.columns if c not in sample_cols]
    return meta_cols, sample_cols


def build_matrix(df, sample_cols, id_col='Met ID'):
    if id_col in df.columns:
        features = df[id_col].fillna('').astype(str)
        # replace empty IDs with row index
        features = features.where(features != '', df.index.astype(str))
    else:
        features = df.index.astype(str)
    mat = df[sample_cols].copy()
    mat.index = features
    mat = mat.apply(pd.to_numeric, errors='coerce')
    mat = mat.replace([np.inf, -np.inf], np.nan)
    return mat


def sanitize_numeric_matrix(mat, name='matrix'):
    mat2 = mat.copy()
    mat2 = mat2.apply(pd.to_numeric, errors='coerce')
    inf_count = int(np.isinf(mat2.values).sum())
    nan_count_before = int(np.isnan(mat2.values).sum())
    if inf_count > 0:
        mat2 = mat2.replace([np.inf, -np.inf], np.nan)
    nan_count_after = int(np.isnan(mat2.values).sum())
    print(f"[QC] {name}: inf={inf_count}, nan_before={nan_count_before}, nan_after={nan_count_after}")
    return mat2


def build_sample_metadata(sample_cols):
    rows = []
    for s in sample_cols:
        m = re.match(r'^(SLE|HC|QC)[-_]?(\d+)', str(s), re.IGNORECASE)
        if m:
            grp = m.group(1).upper()
        else:
            grp = 'UNK'
        rows.append({'sample': s, 'group': grp})
    return pd.DataFrame(rows).set_index('sample')


def qc_normalize(mat, sample_meta):
    mat = sanitize_numeric_matrix(mat, name='raw_for_qc')
    qc_samples = sample_meta[sample_meta['group']=='QC'].index.tolist()
    eps = 1e-9
    if len(qc_samples) == 0:
        # fallback to median normalization across all samples
        ref = mat.median(axis=1)
    else:
        ref = mat[qc_samples].median(axis=1)
    ref = ref.replace(0, eps)
    mat_norm = mat.div(ref, axis=0)
    mat_norm = mat_norm.replace([np.inf, -np.inf], np.nan)
    return mat_norm


def transform_and_scale(mat):
    mat = sanitize_numeric_matrix(mat, name='normalized_matrix')
    non_positive = int((mat <= 0).sum().sum())
    if non_positive > 0:
        print(f"[QC] normalized_matrix: detected {non_positive} non-positive values, set to NaN before log2.")
    mat = mat.mask(mat <= 0, np.nan)
    # Silence expected invalid log warnings from NaN/filtered values.
    with np.errstate(invalid="ignore", divide="ignore"):
        mat_log_np = np.log2(mat.to_numpy(dtype=float))
    mat_log = pd.DataFrame(mat_log_np, index=mat.index, columns=mat.columns)
    mat_log = mat_log.replace([np.inf, -np.inf], np.nan)
    row_median = mat_log.median(axis=1)
    mat_log_filled = mat_log.T.fillna(row_median).T.fillna(0.0)
    scaler = StandardScaler()
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="invalid value encountered in cast", category=RuntimeWarning)
        mat_z = pd.DataFrame(
            scaler.fit_transform(mat_log_filled.T).T,
            index=mat_log.index,
            columns=mat_log.columns,
        )
    return mat_log, mat_z


def filter_features(mat, threshold=0.5):
    # drop if more than threshold fraction of samples are NA or zero
    na_frac = mat.isna().mean(axis=1)
    zero_frac = (mat==0).mean(axis=1)
    keep = (na_frac <= threshold) & (zero_frac <= threshold)
    out = mat.loc[keep].copy()
    out = out.replace([np.inf, -np.inf], np.nan)
    return out


def impute_for_modeling(mat):
    out = mat.copy().replace([np.inf, -np.inf], np.nan)
    row_median = out.median(axis=1)
    out = out.T.fillna(row_median).T.fillna(0.0)
    return out


def differential_test(mat, sample_meta):
    sle = sample_meta[sample_meta['group']=='SLE'].index.tolist()
    hc = sample_meta[sample_meta['group']=='HC'].index.tolist()
    res = []
    for feat, row in mat.iterrows():
        x = row[sle].dropna().values
        y = row[hc].dropna().values
        if len(x) < 3 or len(y) < 3:
            p = np.nan
            stat = np.nan
        else:
            stat, p = stats.ttest_ind(x, y, equal_var=False)
        mean_sle = np.nanmean(row[sle])
        mean_hc = np.nanmean(row[hc])
        fc = (mean_sle + 1e-9) / (mean_hc + 1e-9)
        res.append((feat, stat, p, mean_sle, mean_hc, fc))
    res_df = pd.DataFrame(res, columns=['feature','t_stat','pvalue','mean_sle','mean_hc','fold_change']).set_index('feature')
    res_df['neg_log10_p'] = -np.log10(res_df['pvalue'].replace(0, np.nan))
    res_df['fdr'] = multipletests(res_df['pvalue'].fillna(1), method='fdr_bh')[1]
    return res_df


def compute_auc(mat, sample_meta):
    sle = sample_meta[sample_meta['group']=='SLE'].index.tolist()
    hc = sample_meta[sample_meta['group']=='HC'].index.tolist()
    y = sample_meta['group'].replace({'SLE':1,'HC':0}).dropna()
    common = y.index.tolist()
    aucs = {}
    for feat, row in mat.iterrows():
        vals = row[common].values
        try:
            auc = roc_auc_score(y.loc[common].map({'SLE':1,'HC':0}), vals)
        except Exception:
            auc = np.nan
        aucs[feat] = auc
    return pd.Series(aucs, name='AUC')


def multivariate_models(X, y):
    X = X.replace([np.inf, -np.inf], np.nan).fillna(0.0)
    # LASSO with repeated CV
    rkf = RepeatedStratifiedKFold(n_splits=5, n_repeats=3, random_state=42)
    lasso = LogisticRegressionCV(
        Cs=10,
        penalty='l1',
        solver='saga',
        cv=rkf,
        scoring='roc_auc',
        max_iter=5000,
        n_jobs=1,
        use_legacy_attributes=True,
    )
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="'.*penalty.*deprecated.*'", category=FutureWarning)
        warnings.filterwarnings("ignore", message=".*default value for l1_ratios.*", category=FutureWarning)
        warnings.filterwarnings("ignore", category=ConvergenceWarning)
        lasso.fit(X, y)
    coef = pd.Series(lasso.coef_.ravel(), index=X.columns)

    rf = RandomForestClassifier(n_estimators=500, random_state=42)
    rf_auc = cross_val_score(rf, X, y, cv=rkf, scoring='roc_auc', n_jobs=1)

    return coef, rf_auc


def main():
    print('Loading metabolomics...')
    df = load_metab('metab.xlsx')
    meta_cols, sample_cols = detect_samples(df)
    print(f'{len(sample_cols)} samples detected; {len(meta_cols)} metadata cols')
    sample_meta = build_sample_metadata(sample_cols)
    sample_meta.to_csv(os.path.join(OUTDIR,'sample_metadata.csv'))

    mat = build_matrix(df, sample_cols)
    mat = sanitize_numeric_matrix(mat, name='raw_matrix')
    mat.to_csv(os.path.join(OUTDIR,'raw_matrix.tsv'), sep='\t')

    # QC normalization
    mat_norm = qc_normalize(mat, sample_meta)
    mat_norm.to_csv(os.path.join(OUTDIR,'normalized_matrix.tsv'), sep='\t')

    # transform & scale
    mat_log, mat_z = transform_and_scale(mat_norm)
    mat_log.to_csv(os.path.join(OUTDIR,'log_matrix.tsv'), sep='\t')
    mat_z.to_csv(os.path.join(OUTDIR,'z_matrix.tsv'), sep='\t')

    # filter
    mat_filt = filter_features(mat_log, threshold=0.5)
    mat_filt.to_csv(os.path.join(OUTDIR,'filtered_matrix.tsv'), sep='\t')
    mat_model = impute_for_modeling(mat_filt)
    mat_model.to_csv(os.path.join(OUTDIR,'filtered_matrix_model_ready.tsv'), sep='\t')

    # differential
    diff = differential_test(mat_model, sample_meta)
    diff = diff.sort_values('fdr')
    diff.to_csv(os.path.join(OUTDIR,'differential_results.tsv'), sep='\t')

    # AUC per feature
    aucs = compute_auc(mat_model, sample_meta)
    aucs.to_csv(os.path.join(OUTDIR,'feature_auc.tsv'), sep='\t')

    # merge
    merged = diff.join(aucs)
    merged.to_csv(os.path.join(OUTDIR,'candidates_single_marker.tsv'), sep='\t')

    # multivariate modeling
    # prepare X,y
    common_samples = sample_meta[sample_meta['group'].isin(['SLE','HC'])].index.tolist()
    X = mat_model[common_samples].T
    y = sample_meta.loc[common_samples,'group'].map({'SLE':1,'HC':0}).astype(int)

    scaler = StandardScaler()
    X = X.replace([np.inf, -np.inf], np.nan).fillna(0.0)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="invalid value encountered in cast", category=RuntimeWarning)
        Xs = pd.DataFrame(scaler.fit_transform(X), index=X.index, columns=X.columns)

    coef, rf_auc = multivariate_models(Xs, y)
    coef.to_csv(os.path.join(OUTDIR,'lasso_coefficients.tsv'), sep='\t')
    pd.Series(rf_auc, name='rf_auc').to_csv(os.path.join(OUTDIR,'rf_cv_aucs.tsv'), sep='\t')

    # save top features
    selected = coef[coef.abs()>1e-6].sort_values(key=lambda x: x.abs(), ascending=False)
    selected.to_csv(os.path.join(OUTDIR,'multivariate_selected_features.tsv'), sep='\t')

    # basic plots
    # volcano
    plt.figure(figsize=(6,5))
    fc = pd.to_numeric(merged['fold_change'], errors='coerce')
    log2_fc = pd.Series(np.nan, index=fc.index, dtype=float)
    pos_mask = fc > 0
    with np.errstate(invalid="ignore", divide="ignore"):
        log2_fc.loc[pos_mask] = np.log2(fc.loc[pos_mask].to_numpy(dtype=float))
    sns.scatterplot(x=log2_fc, y=merged['neg_log10_p'], hue=(merged['fdr']<0.05), legend=False, s=10)
    plt.xlabel('log2(Fold Change)')
    plt.ylabel('-log10(p)')
    plt.title('Volcano')
    plt.tight_layout()
    plt.savefig(os.path.join(OUTDIR,'volcano.png'), dpi=200)
    plt.close()

    # heatmap of top 30
    top = merged.head(30).index.tolist()
    sns.clustermap(mat_log.loc[top, common_samples], cmap='vlag', standard_scale=0)
    plt.savefig(os.path.join(OUTDIR,'heatmap_top30.png'), dpi=200)

    # report summary
    summary = {
        'n_samples': len(sample_cols),
        'n_features': mat.shape[0],
        'n_filtered_features': mat_filt.shape[0],
        'n_significant_fdr05': int((merged['fdr']<0.05).sum()),
        'rf_cv_auc_mean': float(np.mean(rf_auc)),
    }
    with open(os.path.join(OUTDIR,'summary.json'),'w',encoding='utf-8') as fh:
        json.dump(summary, fh, ensure_ascii=False, indent=2)

    print('Analysis complete. Outputs in', OUTDIR)


if __name__ == '__main__':
    main()
