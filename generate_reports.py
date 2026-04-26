import pandas as pd
import matplotlib.pyplot as plt
import json

OUT_MD = 'outputs/report.md'
FIG_DIR = 'outputs/report_figs'

def main():
    import os
    os.makedirs(FIG_DIR, exist_ok=True)
    # load summaries
    prot_sum = None
    try:
        prot_sum = pd.read_json('outputs_prot/summary.json', typ='series')
    except Exception:
        prot_sum = None
    adv = None
    try:
        with open('outputs_advanced/nested_cv_results.json') as f:
            adv = json.load(f)
    except Exception:
        adv = None

    # make simple plots: top 10 metab SHAP and prot SHAP if available
    try:
        shapm = pd.read_csv('outputs_advanced/xgb_shap_feature_importance.tsv', sep='\t', index_col=0)
        shapm.head(10).plot.barh(figsize=(6,4), legend=False)
        plt.title('Top metabolite SHAP')
        plt.tight_layout()
        plt.savefig(FIG_DIR + '/top_metab_shap.png', dpi=150)
        plt.close()
    except Exception:
        pass

    try:
        shapp = pd.read_csv('outputs_prot/xgb_shap_prot.tsv', sep='\t', index_col=0)
        shapp.head(10).plot.barh(figsize=(6,4), legend=False)
        plt.title('Top protein SHAP')
        plt.tight_layout()
        plt.savefig(FIG_DIR + '/top_prot_shap.png', dpi=150)
        plt.close()
    except Exception:
        pass

    # write markdown report
    with open(OUT_MD,'w') as f:
        f.write('# SLE 多组学整合报告\n\n')
        f.write('## 方法概述\n')
        f.write('- QC 基准归一化，log2(1+x) 变换，特征过滤（默认 50% 缺失）\n')
        f.write('- 单变量 Welch t 检验 + BH FDR，多重检验校正\n')
        f.write('- 多变量建模：LASSO logistic, RandomForest, XGBoost + SHAP\n')
        f.write('- 面板稳定性：嵌套交叉验证与重复外层划分（见 nested_cv_strict）\n')
        f.write('\n## 关键结果\n')
        if prot_sum is not None:
            f.write(f'- 蛋白组学：{prot_sum.get("n_features","?")} features, {prot_sum.get("n_samples","?")} samples\n')
        if adv is not None:
            f.write(f'- 嵌套 CV 简要：AUC mean = {adv.get("lr_cv_auc_mean","?")}, xgb mean = {adv.get("xgb_cv_auc_mean","?")}\n')
        f.write('\n## 文件与图\n')
        f.write('- Top 代谢物注释： [outputs_advanced/metabolite_annotations.tsv](outputs_advanced/metabolite_annotations.tsv)\n')
        f.write('- PubMed 汇总： [outputs_advanced/pubmed_expanded.tsv](outputs_advanced/pubmed_expanded.tsv)\n')
        f.write('- 最终面板预测与 ROC： [outputs_advanced/final_panel.csv](outputs_advanced/final_panel.csv), [outputs_advanced/final_panel_roc.png](outputs_advanced/final_panel_roc.png)\n')
        f.write('- 多组学 SHAP： [outputs_multiomics/xgb_shap_multi.tsv](outputs_multiomics/xgb_shap_multi.tsv)\n')
        f.write('\n## 局限与建议\n')
        f.write('- 当前模型需独立验证；建议外部队列或留出集进行验证。\n')
    print('Report generated:', OUT_MD)

if __name__=='__main__':
    main()
import os
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import roc_curve, auc, roc_auc_score

OUT = 'outputs'
os.makedirs(OUT, exist_ok=True)

def safe_read_tsv(path):
    return pd.read_csv(path, sep='\t', index_col=0)

def main():
    # read available files
    sample_meta = pd.read_csv(os.path.join(OUT,'sample_metadata.csv'), index_col=0)
    raw = safe_read_tsv(os.path.join(OUT,'raw_matrix.tsv'))
    filtered = safe_read_tsv(os.path.join(OUT,'filtered_matrix.tsv'))
    diff = safe_read_tsv(os.path.join(OUT,'differential_results.tsv'))
    candidates = safe_read_tsv(os.path.join(OUT,'candidates_single_marker.tsv'))
    # try log matrix if exists
    log_path = os.path.join(OUT,'log_matrix.tsv')
    if os.path.exists(log_path):
        mat_log = safe_read_tsv(log_path)
    else:
        mat_log = raw.copy()

    # summary
    n_samples = sample_meta.shape[0]
    n_features = raw.shape[0]
    n_filtered = filtered.shape[0]
    n_signif = int((diff['fdr']<0.05).sum()) if 'fdr' in diff.columns else None
    rf_mean = None
    rf_path = os.path.join(OUT,'rf_cv_aucs.tsv')
    if os.path.exists(rf_path):
        try:
            rf = pd.read_csv(rf_path, sep='\t', index_col=0, squeeze=True)
            rf_mean = float(rf.mean())
        except Exception:
            rf_mean = None

    summary = {
        'n_samples': int(n_samples),
        'n_features': int(n_features),
        'n_filtered_features': int(n_filtered),
        'n_significant_fdr05': int(n_signif) if n_signif is not None else None,
        'rf_cv_auc_mean': rf_mean
    }
    with open(os.path.join(OUT,'summary.json'),'w',encoding='utf-8') as fh:
        json.dump(summary, fh, ensure_ascii=False, indent=2)

    # export top100 candidates
    top100 = candidates.head(100)
    top100.to_csv(os.path.join(OUT,'top100_candidates.tsv'), sep='\t')

    # prepare labels
    labels = sample_meta['group'].replace({'SLE':'SLE','HC':'HC'})
    common = labels[labels.isin(['SLE','HC'])].index.tolist()

    # ensure mat_log rows are features and columns are samples
    # if mat_log has features as rows (index) then keep
    # top features list from top100 index
    top_feats = top100.index.tolist()

    # ROC for top5
    top5 = top_feats[:5]
    plt.figure(figsize=(6,6))
    for feat in top5:
        if feat not in mat_log.index:
            continue
        vals = mat_log.loc[feat, common].astype(float)
        y_true = labels.loc[common].map({'SLE':1,'HC':0}).values
        try:
            fpr, tpr, _ = roc_curve(y_true, vals)
            roc_auc = auc(fpr, tpr)
            plt.plot(fpr, tpr, label=f"{feat} (AUC={roc_auc:.2f})")
        except Exception:
            continue
    plt.plot([0,1],[0,1],'k--',linewidth=0.6)
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC - top5 candidates')
    plt.legend(loc='lower right', fontsize='small')
    plt.tight_layout()
    plt.savefig(os.path.join(OUT,'roc_top5.png'), dpi=200)
    plt.close()

    # heatmap top100
    heat_feats = [f for f in top_feats if f in mat_log.index][:100]
    if len(heat_feats) > 0:
        data = mat_log.loc[heat_feats, common]
        # z-score rows
        data_z = data.sub(data.mean(axis=1), axis=0).div(data.std(axis=1).replace(0,1), axis=0)
        sns.clustermap(data_z, cmap='vlag', col_cluster=True, row_cluster=True, figsize=(8,10))
        plt.savefig(os.path.join(OUT,'heatmap_top100.png'), dpi=200)
        plt.close()

    print('Reports generated: summary.json, top100_candidates.tsv, roc_top5.png, heatmap_top100.png')


if __name__ == '__main__':
    main()
