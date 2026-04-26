import os
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from preprocess_and_analyze import (
    load_metab, detect_samples, build_matrix, build_sample_metadata,
    qc_normalize, transform_and_scale, differential_test, compute_auc,
    multivariate_models
)


OUT = 'outputs'
OUT_STRICT = 'outputs_strict'
os.makedirs(OUT_STRICT, exist_ok=True)


def make_thumbnail(src, dst, scale=0.25):
    try:
        img = mpimg.imread(src)
        h, w = img.shape[0], img.shape[1]
        fig = plt.figure(frameon=False)
        fig.set_size_inches(w*scale/100, h*scale/100)
        ax = plt.Axes(fig, [0., 0., 1., 1.])
        ax.set_axis_off()
        fig.add_axes(ax)
        ax.imshow(img)
        plt.savefig(dst, dpi=100)
        plt.close(fig)
        return True
    except Exception:
        return False


def export_top5_auc():
    auc_path = os.path.join(OUT, 'feature_auc.tsv')
    candidates = pd.read_csv(os.path.join(OUT,'top100_candidates.tsv'), sep='\t', index_col=0)
    if os.path.exists(auc_path):
        aucs = pd.read_csv(auc_path, sep='\t', index_col=0)
        # if single-column dataframe, squeeze
        if aucs.shape[1] == 1:
            aucs = aucs.iloc[:,0]
        top5 = candidates.head(5).index.tolist()
        out = aucs.reindex(top5)
        out.to_csv(os.path.join(OUT,'top5_auc.tsv'), sep='\t', header=True)
        return True
    return False


def rerun_strict(threshold=0.2):
    # load original raw data and sample meta
    df = load_metab('metab.xlsx')
    meta_cols, sample_cols = detect_samples(df)
    sample_meta = build_sample_metadata(sample_cols)
    mat = build_matrix(df, sample_cols)

    # QC normalize
    mat_norm = qc_normalize(mat, sample_meta)
    # log transform
    mat_log, _ = transform_and_scale(mat_norm)

    # stricter filter: remove features with > threshold fraction NA or zero
    na_frac = mat_log.isna().mean(axis=1)
    zero_frac = (mat_log==0).mean(axis=1)
    keep = (na_frac <= threshold) & (zero_frac <= threshold)
    mat_filt = mat_log.loc[keep]

    # differential
    diff = differential_test(mat_filt, sample_meta)
    diff = diff.sort_values('fdr')
    diff.to_csv(os.path.join(OUT_STRICT,'differential_results_strict.tsv'), sep='\t')

    # AUC
    aucs = compute_auc(mat_filt, sample_meta)
    aucs.to_csv(os.path.join(OUT_STRICT,'feature_auc_strict.tsv'), sep='\t')

    # multivariate
    common_samples = sample_meta[sample_meta['group'].isin(['SLE','HC'])].index.tolist()
    X = mat_filt[common_samples].T.fillna(0)
    y = sample_meta.loc[common_samples,'group'].map({'SLE':1,'HC':0}).astype(int)

    coef, rf_auc = multivariate_models(X, y)
    coef.to_csv(os.path.join(OUT_STRICT,'lasso_coefficients_strict.tsv'), sep='\t')
    pd.Series(rf_auc, name='rf_auc').to_csv(os.path.join(OUT_STRICT,'rf_cv_aucs_strict.tsv'), sep='\t')

    # save top candidates
    merged = diff.join(aucs)
    merged.to_csv(os.path.join(OUT_STRICT,'candidates_single_marker_strict.tsv'), sep='\t')

    # write a brief summary
    summary = {
        'n_samples': int(sample_meta.shape[0]),
        'n_features': int(mat.shape[0]),
        'n_filtered_features': int(mat_filt.shape[0]),
        'n_significant_fdr05': int((merged['fdr']<0.05).sum())
    }
    with open(os.path.join(OUT_STRICT,'summary.json'),'w',encoding='utf-8') as fh:
        json.dump(summary, fh, ensure_ascii=False, indent=2)

    return True


def main():
    # make thumbnails
    make_thumbnail(os.path.join(OUT,'roc_top5.png'), os.path.join(OUT,'roc_top5_thumb.png'))
    make_thumbnail(os.path.join(OUT,'heatmap_top100.png'), os.path.join(OUT,'heatmap_top100_thumb.png'))
    # export top5 auc
    export_top5_auc()
    # rerun strict filtering
    rerun_strict(threshold=0.2)
    print('Thumbnails and strict rerun complete. Outputs in', OUT, 'and', OUT_STRICT)


if __name__ == '__main__':
    main()
