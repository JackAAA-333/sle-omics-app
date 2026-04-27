import json
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

OUT_ROOT = Path("outputs")
OUT_MD = OUT_ROOT / "report.md"
FIG_DIR = OUT_ROOT / "report_figs"


def _read_json(path: Path):
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None


def _read_table(path: Path, sep: str = "\t", index_col=None):
    if not path.exists():
        return None
    try:
        return pd.read_csv(path, sep=sep, index_col=index_col)
    except Exception:
        return None


def _plot_top_shap(path: Path, title: str, out_png: Path, topk: int = 12):
    series = None
    try:
        df = pd.read_csv(path, sep="\t", index_col=0)
        if df.shape[1] >= 1:
            series = pd.to_numeric(df.iloc[:, 0], errors="coerce").dropna()
        else:
            series = pd.to_numeric(df.squeeze(), errors="coerce").dropna()
    except Exception:
        try:
            series = pd.read_csv(path, sep="\t", index_col=0, header=None).iloc[:, 0]
            series = pd.to_numeric(series, errors="coerce").dropna()
        except Exception:
            series = None
    if series is None or len(series) == 0:
        return False

    top = series.sort_values(ascending=False).head(topk).sort_values(ascending=True)
    plt.figure(figsize=(8, 5.5))
    plt.barh(top.index.astype(str), top.values, color="#2a7fff")
    plt.xlabel("Mean |SHAP|")
    plt.title(title)
    plt.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_png, dpi=220)
    plt.close()
    return True


def _gene_annotation_summary():
    p = _read_table(Path("outputs_prot/protein_gene_annotations.tsv"), sep="\t")
    if p is None or p.empty:
        return None
    n_total = int(p.shape[0])
    n_ann = int((p["annotation_status"] == "annotated").sum()) if "annotation_status" in p.columns else 0
    return {
        "n_total": n_total,
        "n_annotated": n_ann,
        "ratio": (n_ann / n_total) if n_total > 0 else 0.0,
    }


def main():
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    summary_met = _read_json(Path("outputs/summary.json"))
    summary_prot = _read_json(Path("outputs_prot/summary.json"))
    summary_cand = _read_json(Path("outputs_candidate/panel_summary.json"))
    summary_lit = _read_json(Path("outputs_literature/literature_meta_recent3y_summary.json"))
    summary_causal = _read_json(Path("outputs_causal/causal_summary.json"))
    gene_stat = _gene_annotation_summary()

    img_met_shap = FIG_DIR / "top_metabolite_shap.png"
    img_prot_shap = FIG_DIR / "top_protein_shap.png"
    has_met_img = _plot_top_shap(
        Path("outputs_advanced/xgb_shap_feature_importance.tsv"),
        "Top Metabolite Predictive Contributors",
        img_met_shap,
    )
    has_prot_img = _plot_top_shap(
        Path("outputs_prot/xgb_shap_prot.tsv"),
        "Top Protein Predictive Contributors",
        img_prot_shap,
    )

    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    with OUT_MD.open("w", encoding="utf-8") as f:
        f.write("# 全流程非线性自动组学分析报告\n\n")
        f.write("## 执行概览\n")
        if summary_met:
            f.write(
                f"- 代谢组样本/特征：{summary_met.get('n_samples', 'N/A')} / {summary_met.get('n_features', 'N/A')}。\n"
            )
        if summary_prot:
            f.write(
                f"- 蛋白组样本/特征：{summary_prot.get('n_samples', 'N/A')} / {summary_prot.get('n_features', 'N/A')}。\n"
            )
        if summary_lit:
            f.write(
                f"- 近三年文献命中：{summary_lit.get('total_recent_3y_hits', 'N/A')} 篇，"
                f"覆盖 marker 数：{summary_lit.get('n_markers_with_recent_evidence', 'N/A')}。\n"
            )
        if summary_cand:
            f.write(
                f"- 候选面板数：{summary_cand.get('n_panel_candidates', 'N/A')}，"
                f"严格清单数：{summary_cand.get('n_strict_candidates', 'N/A')}。\n"
            )
        if summary_causal and isinstance(summary_causal, dict):
            cfg = summary_causal.get("config", {})
            if cfg:
                f.write(
                    f"- 因果推断配置：DML {cfg.get('n_splits_dml', 'N/A')}折，"
                    f"Bootstrap {cfg.get('n_bootstrap', 'N/A')}次。\n"
                )
        f.write("\n## 基因注释完整性\n")
        if gene_stat:
            f.write(
                f"- 蛋白 ID 注释到基因符号：{gene_stat['n_annotated']}/{gene_stat['n_total']} "
                f"({gene_stat['ratio']*100:.1f}%)。\n"
            )
            f.write("- 注释明细：`outputs_prot/protein_gene_annotations.tsv`\n")
            f.write("- 差异结果注释版：`outputs_prot/differential_prot_annotated.tsv`\n")
            f.write("- SHAP 注释版：`outputs_prot/xgb_shap_prot_annotated.tsv`\n")
        else:
            f.write("- 本次未检测到蛋白注释统计文件（可能未上传蛋白组数据）。\n")

        f.write("\n## 图形结果总览\n")
        if has_met_img:
            f.write("![Top Metabolite SHAP](report_figs/top_metabolite_shap.png)\n\n")
        if has_prot_img:
            f.write("![Top Protein SHAP](report_figs/top_protein_shap.png)\n\n")
        f.write("![增强火山图](../outputs_candidate/figs/enhanced_volcano_causal.png)\n\n")
        f.write("![推荐面板热图](../outputs_candidate/figs/heatmap_recommended_panel.png)\n\n")
        f.write("![严格面板热图](../outputs_candidate/figs/heatmap_strict_panel.png)\n\n")
        f.write("![推荐面板 ROC](../outputs_candidate/figs/roc_recommended_panel.png)\n\n")
        f.write("![严格面板 ROC](../outputs_candidate/figs/roc_strict_panel.png)\n\n")

        f.write("## 核心结果文件\n")
        f.write("- 候选开发报告：`outputs_candidate/联合试剂盒候选开发报告.md`\n")
        f.write("- 严格临床转化清单：`outputs_candidate/strict_translational_candidates.tsv`\n")
        f.write("- 候选性能指标（AUC/F1）：`outputs_candidate/panel_model_metrics.tsv`\n")
        f.write("- 因果结果：`outputs_causal/causal_marker_effects.tsv`\n")
        f.write("- 文献荟萃：`outputs_literature/literature_meta_recent3y.tsv`\n")
        f.write("- 多组学 SHAP：`outputs_multiomics/xgb_shap_multi.tsv`\n")

        f.write("\n## 解释与建议\n")
        f.write("- 本报告以“统计 + 机器学习 + 因果 + 文献”四证据汇总排序。\n")
        f.write("- 严格清单优先用于临床转化候选，但仍建议外部队列验证后再进入试剂开发阶段。\n")

    print("Report generated:", OUT_MD)


if __name__ == "__main__":
    main()
