import json
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import yaml
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import auc, f1_score, roc_curve
from sklearn.model_selection import StratifiedKFold, cross_val_predict

OUT_DIR = Path("outputs_candidate")
FIG_DIR = OUT_DIR / "figs"
DEFAULT_THRESHOLD_PRESETS = {
    "宽松": (0.05, 1.5),
    "标准": (0.08, 2.0),
    "严格": (0.12, 3.0),
    "auto": (None, None),
}
DEFAULT_PANEL_EVAL_CFG = {
    "cv_folds": 5,
    "rf_n_estimators": 500,
    "random_state": 42,
}


def _load_config():
    cfg_path = Path("config.yaml")
    if not cfg_path.exists():
        return DEFAULT_THRESHOLD_PRESETS.copy(), DEFAULT_PANEL_EVAL_CFG.copy()
    try:
        raw = yaml.safe_load(cfg_path.read_text(encoding="utf-8")) or {}
    except Exception:
        return DEFAULT_THRESHOLD_PRESETS.copy(), DEFAULT_PANEL_EVAL_CFG.copy()

    presets = DEFAULT_THRESHOLD_PRESETS.copy()
    preset_raw = raw.get("strict_threshold_presets", {}) or {}
    for key, zh in [("loose", "宽松"), ("standard", "标准"), ("strict", "严格")]:
        block = preset_raw.get(key, {}) or {}
        d = block.get("min_abs_dml", presets[zh][0])
        e = block.get("min_e_value", presets[zh][1])
        presets[zh] = (float(d), float(e))

    eval_cfg = DEFAULT_PANEL_EVAL_CFG.copy()
    eval_raw = raw.get("panel_evaluation", {}) or {}
    for k in eval_cfg:
        if k in eval_raw:
            eval_cfg[k] = eval_raw[k]
    return presets, eval_cfg


THRESHOLD_PRESETS, PANEL_EVAL_CFG = _load_config()


def _safe_read_table(path: Path, sep: str = "\t", index_col=None):
    if not path.exists():
        return None
    try:
        return pd.read_csv(path, sep=sep, index_col=index_col)
    except Exception:
        return None


def _read_series(path: Path, sep: str = "\t"):
    if not path.exists():
        return None
    try:
        s = pd.read_csv(path, sep=sep, index_col=0, header=None).iloc[:, 0]
        s = pd.to_numeric(s, errors="coerce").dropna()
        return s
    except Exception:
        try:
            s = pd.read_csv(path, sep=sep, index_col=0).iloc[:, 0]
            s = pd.to_numeric(s, errors="coerce").dropna()
            return s
        except Exception:
            return None


def _safe_read_matrix(path: Path):
    if not path.exists():
        return None
    try:
        return pd.read_csv(path, sep="\t", index_col=0)
    except Exception:
        return None


def _norm01(series: pd.Series) -> pd.Series:
    if series is None or len(series) == 0:
        return pd.Series(dtype=float)
    s = pd.to_numeric(series, errors="coerce")
    s = s.fillna(0.0)
    vmax = float(s.max())
    vmin = float(s.min())
    if vmax <= vmin:
        return pd.Series(1.0, index=s.index, dtype=float)
    return (s - vmin) / (vmax - vmin)


def _prepare_candidates():
    rows = []

    prot_shap_annot = _safe_read_table(Path("outputs_prot/xgb_shap_prot_annotated.tsv"), sep="\t")
    prot_diff = _safe_read_table(Path("outputs_prot/differential_prot.tsv"), sep="\t", index_col=0)
    if prot_shap_annot is not None and not prot_shap_annot.empty:
        p = prot_shap_annot.copy()
        if "protein_id" not in p.columns:
            p["protein_id"] = p.iloc[:, 0].astype(str)
        p["feature_id"] = p["protein_id"].astype(str)
        p["marker"] = p.get("gene_symbol", "").astype(str).str.strip()
        p.loc[p["marker"].isin(["", "nan", "None"]), "marker"] = p["feature_id"]
        p["modality"] = "protein"
        p["shap_value"] = pd.to_numeric(p.get("mean_abs_shap", np.nan), errors="coerce")
        if prot_diff is not None and not prot_diff.empty:
            p["pvalue"] = p["feature_id"].map(prot_diff.get("pvalue", pd.Series(dtype=float)))
            p["fdr"] = p["feature_id"].map(prot_diff.get("fdr_bh", pd.Series(dtype=float)))
            p["auc"] = p["feature_id"].map(prot_diff.get("auc", pd.Series(dtype=float)))
        rows.append(p[["modality", "feature_id", "marker", "shap_value", "pvalue", "fdr", "auc"]])

    met_shap = _read_series(Path("outputs_advanced/xgb_shap_feature_importance.tsv"), sep="\t")
    met_diff_basic = _safe_read_table(Path("outputs/differential_results.tsv"), sep="\t", index_col=0)
    met_diff_adv = _safe_read_table(Path("outputs_advanced/differential_results_advanced.tsv"), sep="\t", index_col=0)
    met_name_map_df = None
    try:
        met_name_map_df = pd.read_excel("metab.xlsx")
    except Exception:
        met_name_map_df = None
    met_name_map = {}
    if met_name_map_df is not None and "Met ID" in met_name_map_df.columns and "Metabolite Name" in met_name_map_df.columns:
        for _, r in met_name_map_df.iterrows():
            k = str(r.get("Met ID", "")).strip()
            v = str(r.get("Metabolite Name", "")).strip()
            if k and v and v.lower() != "nan":
                met_name_map[k] = v

    if met_shap is not None and len(met_shap) > 0:
        m = pd.DataFrame({"feature_id": met_shap.index.astype(str), "shap_value": met_shap.values})
        m["marker"] = m["feature_id"].map(met_name_map).fillna(m["feature_id"])
        m["modality"] = "metabolite"
        if met_diff_basic is not None and not met_diff_basic.empty:
            m["pvalue"] = m["feature_id"].map(met_diff_basic.get("pvalue", pd.Series(dtype=float)))
            m["fdr"] = m["feature_id"].map(met_diff_basic.get("fdr", pd.Series(dtype=float)))
            m["auc"] = m["feature_id"].map(met_diff_basic.get("AUC", pd.Series(dtype=float)))
            m["fold_change"] = m["feature_id"].map(met_diff_basic.get("fold_change", pd.Series(dtype=float)))
        if met_diff_adv is not None and not met_diff_adv.empty:
            adv_fdr = m["feature_id"].map(met_diff_adv.get("fdr_bh", pd.Series(dtype=float)))
            m["fdr"] = pd.to_numeric(m.get("fdr"), errors="coerce").combine_first(pd.to_numeric(adv_fdr, errors="coerce"))
        rows.append(m[["modality", "feature_id", "marker", "shap_value", "pvalue", "fdr", "auc", "fold_change"]])

    if len(rows) == 0:
        return pd.DataFrame()
    all_df = pd.concat(rows, axis=0, ignore_index=True)
    all_df["marker"] = all_df["marker"].astype(str).str.strip()
    all_df = all_df[(all_df["marker"] != "") & (all_df["marker"].str.lower() != "nan")]
    return all_df


def _attach_literature(df: pd.DataFrame):
    lit = _safe_read_table(Path("outputs_literature/literature_meta_recent3y.tsv"), sep="\t")
    if lit is None or lit.empty:
        df["recent_3y_count"] = 0
        df["pmids"] = ""
        return df
    lit2 = lit.copy()
    lit2["feature_id"] = lit2.get("feature_id", "").astype(str)
    lit2["marker"] = lit2.get("marker", "").astype(str)
    lit2["modality"] = lit2.get("modality", "").astype(str)
    lit2["recent_3y_count"] = pd.to_numeric(lit2.get("recent_3y_count", 0), errors="coerce").fillna(0).astype(int)
    lit2 = lit2.sort_values("recent_3y_count", ascending=False).drop_duplicates(["modality", "feature_id"], keep="first")
    merge_cols = ["modality", "feature_id", "recent_3y_count", "pmids"]
    out = df.merge(lit2[merge_cols], on=["modality", "feature_id"], how="left")
    out["recent_3y_count"] = pd.to_numeric(out["recent_3y_count"], errors="coerce").fillna(0).astype(int)
    out["pmids"] = out["pmids"].fillna("")
    return out


def _attach_causal(df: pd.DataFrame):
    causal = _safe_read_table(Path("outputs_causal/causal_marker_effects.tsv"), sep="\t")
    if causal is None or causal.empty:
        df["dml_effect_theta"] = np.nan
        df["ci95_low"] = np.nan
        df["ci95_high"] = np.nan
        df["or_adjusted"] = np.nan
        df["e_value"] = np.nan
        return df
    c = causal.copy()
    c["modality"] = c["modality"].astype(str)
    c["feature_id"] = c["feature_id"].astype(str)
    merge_cols = ["modality", "feature_id", "dml_effect_theta", "ci95_low", "ci95_high", "or_adjusted", "e_value"]
    out = df.merge(c[merge_cols], on=["modality", "feature_id"], how="left")
    return out


def _add_strict_filter(df: pd.DataFrame, min_abs_dml: float | None = None, min_e_value: float | None = None):
    if df.empty:
        return df
    out = df.copy()
    abs_dml = pd.to_numeric(out.get("dml_effect_theta", np.nan), errors="coerce").abs().fillna(0.0)
    e_val = pd.to_numeric(out.get("e_value", np.nan), errors="coerce").fillna(1.0)
    if min_abs_dml is None:
        min_abs_dml = float(abs_dml.quantile(0.6))
    if min_e_value is None:
        min_e_value = float(max(1.25, e_val.quantile(0.6)))
    out["abs_dml"] = abs_dml
    out["strict_pass"] = (out["abs_dml"] >= min_abs_dml) & (e_val >= min_e_value)
    out["strict_threshold_abs_dml"] = min_abs_dml
    out["strict_threshold_e_value"] = min_e_value
    return out


def _resolve_thresholds(profile: str, min_abs_dml: float | None, min_e_value: float | None):
    if profile in THRESHOLD_PRESETS:
        p_dml, p_e = THRESHOLD_PRESETS[profile]
        if min_abs_dml is None:
            min_abs_dml = p_dml
        if min_e_value is None:
            min_e_value = p_e
    return min_abs_dml, min_e_value


def _score_and_rank(df: pd.DataFrame):
    if df.empty:
        return df
    df = df.copy()
    df["shap_value"] = pd.to_numeric(df["shap_value"], errors="coerce").fillna(0.0)
    df["fdr"] = pd.to_numeric(df.get("fdr", np.nan), errors="coerce")
    df["auc"] = pd.to_numeric(df.get("auc", np.nan), errors="coerce")

    # 证据组件：模型贡献、统计显著性、判别能力、近三年文献支持
    df["pred_component"] = _norm01(df["shap_value"])
    fdr_proxy = -np.log10(df["fdr"].fillna(1.0).clip(lower=1e-12))
    df["stat_component"] = _norm01(fdr_proxy)
    auc_dev = (df["auc"].fillna(0.5) - 0.5).abs() * 2.0
    df["auc_component"] = auc_dev.clip(lower=0, upper=1)
    df["lit_component"] = _norm01(df["recent_3y_count"])
    df["causal_component"] = _norm01(pd.to_numeric(df.get("dml_effect_theta", np.nan), errors="coerce").abs().fillna(0.0))
    e_val_norm = _norm01(pd.to_numeric(df.get("e_value", np.nan), errors="coerce").fillna(1.0))
    df["causal_component"] = 0.7 * df["causal_component"] + 0.3 * e_val_norm

    df["effect_component"] = 0.6 * df["stat_component"] + 0.4 * df["auc_component"]
    df["total_score"] = (
        0.45 * df["pred_component"]
        + 0.20 * df["effect_component"]
        + 0.15 * df["lit_component"]
        + 0.20 * df["causal_component"]
    )
    df = df.sort_values("total_score", ascending=False).reset_index(drop=True)
    df["rank"] = np.arange(1, len(df) + 1)
    return df


def _pick_panel(df: pd.DataFrame, panel_size: int = 10):
    if df.empty:
        return df
    df = df.copy()
    prot = df[df["modality"] == "protein"].head(5)
    met = df[df["modality"] == "metabolite"].head(5)
    picked = pd.concat([prot, met], axis=0).drop_duplicates(["modality", "feature_id"])
    if len(picked) < panel_size:
        remain = df[~df["feature_id"].isin(picked["feature_id"])]
        picked = pd.concat([picked, remain.head(panel_size - len(picked))], axis=0)
    picked = picked.head(panel_size).sort_values("total_score", ascending=False).reset_index(drop=True)
    picked["panel_rank"] = np.arange(1, len(picked) + 1)
    return picked


def _pick_strict_panel(df: pd.DataFrame, panel_size: int = 8):
    if df.empty:
        return df
    strict = df[df.get("strict_pass", False)].copy().sort_values("total_score", ascending=False)
    if strict.empty:
        strict = df.copy().sort_values("causal_component", ascending=False).head(panel_size)
    else:
        prot = strict[strict["modality"] == "protein"].head(panel_size // 2)
        met = strict[strict["modality"] == "metabolite"].head(panel_size // 2)
        strict = pd.concat([prot, met], axis=0).drop_duplicates(["modality", "feature_id"])
        if len(strict) < panel_size:
            remain = df[~df["feature_id"].isin(strict["feature_id"])].sort_values("total_score", ascending=False)
            strict = pd.concat([strict, remain.head(panel_size - len(strict))], axis=0)
    strict = strict.head(panel_size).sort_values("total_score", ascending=False).reset_index(drop=True)
    strict["strict_rank"] = np.arange(1, len(strict) + 1)
    return strict


def _reason_text(row: pd.Series):
    comps = []
    comps.append(f"模型贡献权重={row['pred_component']:.2f}")
    comps.append(f"统计/判别证据={row['effect_component']:.2f}")
    comps.append(f"近三年文献命中={int(row['recent_3y_count'])}篇")
    comps.append(f"因果证据={row.get('causal_component', 0):.2f}")
    auc = row.get("auc", np.nan)
    if pd.notna(auc):
        comps.append(f"AUC={float(auc):.3f}")
    fdr = row.get("fdr", np.nan)
    if pd.notna(fdr):
        comps.append(f"FDR={float(fdr):.2e}")
    theta = row.get("dml_effect_theta", np.nan)
    if pd.notna(theta):
        comps.append(f"DML效应={float(theta):.3f}")
    or_adj = row.get("or_adjusted", np.nan)
    if pd.notna(or_adj):
        comps.append(f"调整OR={float(or_adj):.3f}")
    e_val = row.get("e_value", np.nan)
    if pd.notna(e_val):
        comps.append(f"E-value={float(e_val):.3f}")
    return "；".join(comps)


def _plot_panel(panel_df: pd.DataFrame):
    if panel_df.empty:
        return
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    names = panel_df["marker"].tolist()
    y = np.arange(len(names))

    # 综合评分条形图
    plt.figure(figsize=(10, 6))
    plt.barh(y, panel_df["total_score"], color="#2a7fff")
    plt.yticks(y, names)
    plt.gca().invert_yaxis()
    plt.xlabel("综合候选评分")
    plt.title("联合试剂盒候选分子综合评分")
    plt.tight_layout()
    plt.savefig(FIG_DIR / "panel_total_score.png", dpi=200)
    plt.close()

    # 证据拆解堆叠图
    plt.figure(figsize=(10, 6))
    plt.barh(y, panel_df["pred_component"], label="模型贡献(45%)", color="#1f77b4")
    plt.barh(y, panel_df["effect_component"], left=panel_df["pred_component"], label="统计/判别(20%)", color="#ff7f0e")
    left2 = panel_df["pred_component"] + panel_df["effect_component"]
    plt.barh(y, panel_df["lit_component"], left=left2, label="文献支持(15%)", color="#2ca02c")
    left3 = left2 + panel_df["lit_component"]
    plt.barh(y, panel_df["causal_component"], left=left3, label="因果证据(20%)", color="#d62728")
    plt.yticks(y, names)
    plt.gca().invert_yaxis()
    plt.xlabel("分项证据得分(0-1)")
    plt.title("候选分子证据分项构成")
    plt.legend(loc="lower right")
    plt.tight_layout()
    plt.savefig(FIG_DIR / "panel_evidence_breakdown.png", dpi=200)
    plt.close()


def _plot_enhanced_volcano(all_df: pd.DataFrame):
    if all_df.empty or "pvalue" not in all_df.columns:
        return
    v = all_df.copy()
    v["pvalue"] = pd.to_numeric(v["pvalue"], errors="coerce")
    v = v[v["pvalue"].notna()].copy()
    if v.empty:
        return
    x = pd.to_numeric(v.get("dml_effect_theta", np.nan), errors="coerce").fillna(0.0)
    y = -np.log10(v["pvalue"].clip(lower=1e-300))
    hue = v.get("strict_pass", False)
    plt.figure(figsize=(10, 6))
    sns.scatterplot(x=x, y=y, hue=hue, palette={True: "#d62728", False: "#1f77b4"}, s=45, edgecolor="white", linewidth=0.3)
    top = v.sort_values("total_score", ascending=False).head(12)
    for _, r in top.iterrows():
        plt.text(float(r.get("dml_effect_theta", 0)), float(-np.log10(max(float(r.get("pvalue", 1)), 1e-300))), str(r["marker"]), fontsize=8)
    plt.axvline(0, color="grey", linestyle="--", linewidth=0.8)
    plt.xlabel("DML 因果效应 (theta)")
    plt.ylabel("-log10(p-value)")
    plt.title("增强火山图：统计显著性 × 因果效应")
    plt.tight_layout()
    plt.savefig(FIG_DIR / "enhanced_volcano_causal.png", dpi=220)
    plt.close()


def _build_feature_matrix():
    meta_path = Path("outputs/sample_metadata.csv")
    if not meta_path.exists():
        return None, None
    meta = pd.read_csv(meta_path, index_col=0)
    meta = meta[meta["group"].isin(["SLE", "HC"])].copy()
    if meta.empty:
        return None, None
    y = (meta["group"] == "SLE").astype(int)
    met = _safe_read_matrix(Path("outputs/filtered_matrix.tsv"))
    prot = _safe_read_matrix(Path("outputs_prot/filtered_prot_matrix.tsv"))
    return (met, prot, y)


def _evaluate_panel_metrics(panel_df: pd.DataFrame, prefix: str):
    data = _build_feature_matrix()
    if data is None or panel_df.empty:
        return None
    met, prot, y = data
    samples = y.index.tolist()
    X_parts = []
    for _, r in panel_df.iterrows():
        feat = str(r["feature_id"])
        mod = str(r["modality"])
        mat = prot if mod == "protein" else met
        if mat is None or feat not in mat.index:
            continue
        row = pd.to_numeric(mat.loc[feat], errors="coerce")
        row = row.reindex(samples).fillna(float(row.median() if row.notna().any() else 0))
        X_parts.append(row.rename(f"{mod}:{feat}"))
    if len(X_parts) < 2:
        return None
    X = pd.concat(X_parts, axis=1).astype(float)
    cv = StratifiedKFold(
        n_splits=int(PANEL_EVAL_CFG.get("cv_folds", 5)),
        shuffle=True,
        random_state=int(PANEL_EVAL_CFG.get("random_state", 42)),
    )
    clf = RandomForestClassifier(
        n_estimators=int(PANEL_EVAL_CFG.get("rf_n_estimators", 500)),
        random_state=int(PANEL_EVAL_CFG.get("random_state", 42)),
        class_weight="balanced",
    )
    y_prob = cross_val_predict(clf, X, y.values, cv=cv, method="predict_proba")[:, 1]
    y_pred = (y_prob >= 0.5).astype(int)
    fpr, tpr, _ = roc_curve(y.values, y_prob)
    auc_val = float(auc(fpr, tpr))
    f1_val = float(f1_score(y.values, y_pred))

    plt.figure(figsize=(6, 6))
    plt.plot(fpr, tpr, color="#2a7fff", label=f"AUC={auc_val:.3f}")
    plt.plot([0, 1], [0, 1], "k--", linewidth=0.8)
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title(f"ROC - {prefix} panel")
    plt.legend(loc="lower right")
    plt.tight_layout()
    plt.savefig(FIG_DIR / f"roc_{prefix}_panel.png", dpi=220)
    plt.close()

    return {"panel": prefix, "n_features": int(X.shape[1]), "auc_cv": auc_val, "f1_cv": f1_val}


def _plot_panel_heatmap(panel_df: pd.DataFrame, prefix: str):
    data = _build_feature_matrix()
    if data is None or panel_df.empty:
        return
    met, prot, y = data
    samples = y.index.tolist()
    X_parts = []
    for _, r in panel_df.iterrows():
        feat = str(r["feature_id"])
        mod = str(r["modality"])
        mat = prot if mod == "protein" else met
        if mat is None or feat not in mat.index:
            continue
        row = pd.to_numeric(mat.loc[feat], errors="coerce")
        row = row.reindex(samples).fillna(float(row.median() if row.notna().any() else 0))
        name = str(r["marker"])
        X_parts.append(row.rename(name))
    if len(X_parts) < 2:
        return
    M = pd.concat(X_parts, axis=1).T
    Mz = M.sub(M.mean(axis=1), axis=0).div(M.std(axis=1).replace(0, 1), axis=0)
    col_colors = y.map({1: "#d62728", 0: "#1f77b4"})
    cg = sns.clustermap(
        Mz,
        cmap="vlag",
        col_colors=col_colors,
        col_cluster=True,
        row_cluster=True,
        figsize=(10, 8),
    )
    cg.fig.suptitle(f"Heatmap - {prefix} panel", y=1.02)
    cg.fig.savefig(FIG_DIR / f"heatmap_{prefix}_panel.png", dpi=220, bbox_inches="tight")
    plt.close(cg.fig)


def _write_markdown(panel_df: pd.DataFrame, strict_df: pd.DataFrame, all_df: pd.DataFrame, metrics_df: pd.DataFrame):
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    md = OUT_DIR / "联合试剂盒候选开发报告.md"
    with md.open("w", encoding="utf-8") as f:
        f.write("# 联合试剂盒候选 Biomarker 开发报告\n\n")
        f.write("## 评分策略（含最新因果推断，可解释）\n")
        f.write("- 模型贡献权重（45%）：来自 XGBoost SHAP 绝对值均值，表示分子对预测输出的平均贡献。\n")
        f.write("- 统计/判别证据（20%）：综合 FDR 显著性与 AUC 区分能力。\n")
        f.write("- 近三年文献支持（15%）：PubMed 中 SLE 主题下最近三年的文献命中数量。\n")
        f.write("- 因果证据（20%）：Double Machine Learning (PLR) 估计因果效应 + E-value 稳健性评估。\n\n")
        f.write("## 候选联合面板（推荐）\n")
        if panel_df.empty:
            f.write("暂无可用候选分子（输入数据或上游分析输出不足）。\n")
        else:
            f.write("![综合评分图](figs/panel_total_score.png)\n\n")
            f.write("![证据构成图](figs/panel_evidence_breakdown.png)\n\n")
            f.write("| 排名 | 分子 | 类型 | 综合评分 | SHAP贡献 | 统计/判别 | 文献 | 因果证据 |\n")
            f.write("|---:|---|---|---:|---:|---:|---:|---:|\n")
            for _, r in panel_df.iterrows():
                f.write(
                    f"| {int(r['panel_rank'])} | {r['marker']} | {r['modality']} | "
                    f"{r['total_score']:.3f} | {r['pred_component']:.3f} | {r['effect_component']:.3f} | "
                    f"{int(r['recent_3y_count'])} | {r.get('causal_component', 0):.3f} |\n"
                )
            f.write("\n## 每个候选分子的“为什么”\n")
            for _, r in panel_df.iterrows():
                f.write(f"- **{r['marker']}**（{r['modality']}）: {_reason_text(r)}。\n")

        f.write("\n## 更严格临床转化候选清单（因果阈值过滤）\n")
        if strict_df.empty:
            f.write("暂无通过严格阈值的分子（可在样本量扩大后重跑）。\n")
        else:
            dml_t = float(strict_df["strict_threshold_abs_dml"].iloc[0])
            e_t = float(strict_df["strict_threshold_e_value"].iloc[0])
            f.write(f"- 阈值：`|DML| >= {dml_t:.3f}` 且 `E-value >= {e_t:.3f}`\n\n")
            f.write("| 排名 | 分子 | 类型 | 综合评分 | |DML| | E-value |\n")
            f.write("|---:|---|---|---:|---:|---:|\n")
            for _, r in strict_df.iterrows():
                f.write(
                    f"| {int(r['strict_rank'])} | {r['marker']} | {r['modality']} | "
                    f"{r['total_score']:.3f} | {float(abs(r.get('dml_effect_theta', 0))):.3f} | {float(r.get('e_value', np.nan)):.3f} |\n"
                )

        f.write("\n## 必要性能与可视化结果\n")
        f.write("![增强火山图](figs/enhanced_volcano_causal.png)\n\n")
        f.write("![推荐面板热图](figs/heatmap_recommended_panel.png)\n\n")
        f.write("![严格面板热图](figs/heatmap_strict_panel.png)\n\n")
        f.write("![推荐面板ROC](figs/roc_recommended_panel.png)\n\n")
        f.write("![严格面板ROC](figs/roc_strict_panel.png)\n\n")
        if metrics_df is not None and not metrics_df.empty:
            f.write("| 面板 | 特征数 | CV AUC | CV F1 |\n")
            f.write("|---|---:|---:|---:|\n")
            for _, r in metrics_df.iterrows():
                f.write(f"| {r['panel']} | {int(r['n_features'])} | {r['auc_cv']:.3f} | {r['f1_cv']:.3f} |\n")

        f.write("\n## 产出文件\n")
        f.write("- `outputs_candidate/panel_candidates.tsv`：推荐联合面板清单。\n")
        f.write("- `outputs_candidate/strict_translational_candidates.tsv`：严格临床转化清单（因果阈值过滤）。\n")
        f.write("- `outputs_candidate/all_candidates_scored.tsv`：全部候选分子评分明细。\n")
        f.write("- `outputs_causal/causal_marker_effects.tsv`：DML 因果效应、95%CI、调整 OR 与 E-value。\n")
        f.write("- `outputs_candidate/figs/panel_total_score.png`：综合评分图。\n")
        f.write("- `outputs_candidate/figs/panel_evidence_breakdown.png`：证据分项构成图。\n")
        f.write("- `outputs_candidate/figs/enhanced_volcano_causal.png`：增强火山图。\n")
        f.write("- `outputs_candidate/figs/heatmap_recommended_panel.png` / `heatmap_strict_panel.png`：面板热图。\n")
        f.write("- `outputs_candidate/figs/roc_recommended_panel.png` / `roc_strict_panel.png`：ROC 曲线。\n")
        f.write("- `outputs_candidate/panel_model_metrics.tsv`：AUC/F1 指标汇总。\n")
        f.write("\n## 说明\n")
        f.write("- 本报告用于候选分子优先级排序，不替代外部队列和临床验证。\n")
        f.write("- 若为单组学输入，结果仅基于该组学候选自动生成。\n")


def main(min_abs_dml: float | None = None, min_e_value: float | None = None, strict_profile: str = "标准"):
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    min_abs_dml, min_e_value = _resolve_thresholds(strict_profile, min_abs_dml, min_e_value)
    all_df = _prepare_candidates()
    all_df = _attach_literature(all_df) if not all_df.empty else all_df
    all_df = _attach_causal(all_df) if not all_df.empty else all_df
    all_df = _score_and_rank(all_df) if not all_df.empty else all_df
    all_df = _add_strict_filter(all_df, min_abs_dml=min_abs_dml, min_e_value=min_e_value) if not all_df.empty else all_df
    panel_df = _pick_panel(all_df, panel_size=10) if not all_df.empty else pd.DataFrame()
    strict_df = _pick_strict_panel(all_df, panel_size=8) if not all_df.empty else pd.DataFrame()

    if not all_df.empty:
        all_df.to_csv(OUT_DIR / "all_candidates_scored.tsv", sep="\t", index=False)
    if not panel_df.empty:
        panel_df.to_csv(OUT_DIR / "panel_candidates.tsv", sep="\t", index=False)
        _plot_panel(panel_df)
        _plot_panel_heatmap(panel_df, "recommended")
    if not strict_df.empty:
        strict_df.to_csv(OUT_DIR / "strict_translational_candidates.tsv", sep="\t", index=False)
        _plot_panel_heatmap(strict_df, "strict")
    else:
        pd.DataFrame(columns=["panel_rank", "marker", "modality"]).to_csv(
            OUT_DIR / "panel_candidates.tsv", sep="\t", index=False
        )
        pd.DataFrame(columns=["strict_rank", "marker", "modality"]).to_csv(
            OUT_DIR / "strict_translational_candidates.tsv", sep="\t", index=False
        )

    _plot_enhanced_volcano(all_df if not all_df.empty else pd.DataFrame())
    metrics_rows = []
    m1 = _evaluate_panel_metrics(panel_df, "recommended")
    m2 = _evaluate_panel_metrics(strict_df, "strict")
    if m1:
        metrics_rows.append(m1)
    if m2:
        metrics_rows.append(m2)
    metrics_df = pd.DataFrame(metrics_rows)
    if not metrics_df.empty:
        metrics_df.to_csv(OUT_DIR / "panel_model_metrics.tsv", sep="\t", index=False)
    else:
        pd.DataFrame(columns=["panel", "n_features", "auc_cv", "f1_cv"]).to_csv(
            OUT_DIR / "panel_model_metrics.tsv", sep="\t", index=False
        )

    _write_markdown(panel_df, strict_df, all_df, metrics_df)

    summary = {
        "strict_profile": strict_profile,
        "strict_threshold_abs_dml": None if all_df.empty else float(all_df["strict_threshold_abs_dml"].iloc[0]),
        "strict_threshold_e_value": None if all_df.empty else float(all_df["strict_threshold_e_value"].iloc[0]),
        "n_all_candidates": int(0 if all_df.empty else all_df.shape[0]),
        "n_panel_candidates": int(0 if panel_df.empty else panel_df.shape[0]),
        "n_strict_candidates": int(0 if strict_df.empty else strict_df.shape[0]),
        "panel_top5": [] if panel_df.empty else panel_df.head(5)[["marker", "modality", "total_score"]].to_dict(orient="records"),
    }
    with open(OUT_DIR / "panel_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, ensure_ascii=False, indent=2)
    print("Biomarker candidate report generated in", OUT_DIR)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate biomarker candidate report with strict causal thresholds.")
    parser.add_argument("--min-abs-dml", type=float, default=None, help="Strict filter threshold for |DML|.")
    parser.add_argument("--min-e-value", type=float, default=None, help="Strict filter threshold for E-value.")
    parser.add_argument(
        "--strict-profile",
        type=str,
        default="标准",
        choices=["宽松", "标准", "严格", "auto"],
        help="Threshold preset profile.",
    )
    args = parser.parse_args()
    main(min_abs_dml=args.min_abs_dml, min_e_value=args.min_e_value, strict_profile=args.strict_profile)
