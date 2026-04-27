import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegressionCV
from sklearn.model_selection import KFold
from sklearn.ensemble import RandomForestRegressor

OUT = Path("outputs_causal")
OUT.mkdir(parents=True, exist_ok=True)
DEFAULT_CAUSAL_CFG = {
    # Commonly accepted robust defaults in biomedical ML workflows.
    "n_top_markers_per_modality": 30,
    "n_pc": 10,
    "min_samples_per_marker": 30,
    "n_splits_dml": 5,
    "n_bootstrap": 1000,
    "rf_n_estimators": 500,
    "seed": 42,
}


def _load_causal_cfg():
    cfg_path = Path("config.yaml")
    if not cfg_path.exists():
        return DEFAULT_CAUSAL_CFG.copy()
    try:
        raw = yaml.safe_load(cfg_path.read_text(encoding="utf-8")) or {}
        c = raw.get("causal_inference", {}) or {}
        out = DEFAULT_CAUSAL_CFG.copy()
        out.update({k: c[k] for k in out.keys() if k in c})
        return out
    except Exception:
        return DEFAULT_CAUSAL_CFG.copy()


CAUSAL_CFG = _load_causal_cfg()


def _load_meta():
    p = Path("outputs/sample_metadata.csv")
    if not p.exists():
        return None
    meta = pd.read_csv(p, index_col=0)
    meta = meta[meta["group"].isin(["SLE", "HC"])].copy()
    meta["y"] = (meta["group"] == "SLE").astype(float)
    return meta


def _load_matrix(path: Path):
    if not path.exists():
        return None
    try:
        return pd.read_csv(path, sep="\t", index_col=0)
    except Exception:
        return None


def _load_top_candidates():
    parts = []
    prot = Path("outputs_prot/xgb_shap_prot_annotated.tsv")
    if prot.exists():
        dfp = pd.read_csv(prot, sep="\t")
        if "protein_id" in dfp.columns:
            tmp = pd.DataFrame(
                {
                    "modality": "protein",
                    "feature_id": dfp["protein_id"].astype(str),
                    "marker": dfp.get("gene_symbol", dfp["protein_id"]).astype(str),
                    "shap": pd.to_numeric(dfp.get("mean_abs_shap", 0), errors="coerce").fillna(0),
                }
            )
            parts.append(tmp.head(CAUSAL_CFG["n_top_markers_per_modality"]))
    met = Path("outputs_advanced/xgb_shap_feature_importance.tsv")
    if met.exists():
        s = pd.read_csv(met, sep="\t", index_col=0).iloc[:, 0]
        tmp = pd.DataFrame(
            {"modality": "metabolite", "feature_id": s.index.astype(str), "marker": s.index.astype(str), "shap": s.values}
        ).head(CAUSAL_CFG["n_top_markers_per_modality"])
        parts.append(tmp)
    if len(parts) == 0:
        return pd.DataFrame(columns=["modality", "feature_id", "marker", "shap"])
    return pd.concat(parts, ignore_index=True).drop_duplicates(["modality", "feature_id"])


def _build_confounders(meta, met_mat, prot_mat, n_pc=5):
    samples = meta.index.tolist()
    blocks = []
    for mat in [met_mat, prot_mat]:
        if mat is None or mat.empty:
            continue
        common = [s for s in samples if s in mat.columns]
        if len(common) < 10:
            continue
        X = mat[common].T.fillna(0.0).astype(float)
        k = max(1, min(n_pc, X.shape[1], X.shape[0] - 1))
        pcs = PCA(n_components=k, random_state=42).fit_transform(X)
        pc_df = pd.DataFrame(pcs, index=common, columns=[f"PC_{len(blocks)}_{i+1}" for i in range(k)])
        blocks.append(pc_df.reindex(samples).fillna(0.0))
    if len(blocks) == 0:
        return pd.DataFrame(index=samples)
    Z = pd.concat(blocks, axis=1)
    return Z


def _dml_plr_binary_y(y, t, z, n_splits=5, seed=42, rf_n_estimators=500):
    # DML partially linear model:
    # y = theta*t + g(z) + u
    # t = m(z) + v
    # theta = E[(y-g(z))*(t-m(z))]/E[(t-m(z))^2]
    y = np.asarray(y, dtype=float)
    t = np.asarray(t, dtype=float)
    if z is None or z.shape[1] == 0:
        z = np.zeros((len(y), 1))
    z = np.asarray(z, dtype=float)
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=seed)
    y_res = np.zeros_like(y, dtype=float)
    t_res = np.zeros_like(t, dtype=float)
    for tr, te in kf.split(z):
        gy = RandomForestRegressor(n_estimators=rf_n_estimators, random_state=seed)
        gt = RandomForestRegressor(n_estimators=rf_n_estimators, random_state=seed + 1)
        gy.fit(z[tr], y[tr])
        gt.fit(z[tr], t[tr])
        y_res[te] = y[te] - gy.predict(z[te])
        t_res[te] = t[te] - gt.predict(z[te])
    denom = np.dot(t_res, t_res)
    if abs(denom) < 1e-12:
        return np.nan
    theta = float(np.dot(t_res, y_res) / denom)
    return theta


def _bootstrap_ci(y, t, z, n_boot=1000, seed=42, n_splits=5, rf_n_estimators=500):
    rng = np.random.default_rng(seed)
    n = len(y)
    vals = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        v = _dml_plr_binary_y(
            y[idx],
            t[idx],
            z[idx],
            n_splits=n_splits,
            seed=seed,
            rf_n_estimators=rf_n_estimators,
        )
        if np.isfinite(v):
            vals.append(v)
    if len(vals) < 10:
        return np.nan, np.nan
    lo, hi = np.percentile(vals, [2.5, 97.5])
    return float(lo), float(hi)


def _e_value_from_or(or_val):
    # VanderWeele E-value for point estimate
    # For OR >=1: E = OR + sqrt(OR*(OR-1))
    if not np.isfinite(or_val) or or_val <= 0:
        return np.nan
    or_use = or_val if or_val >= 1 else 1.0 / or_val
    return float(or_use + np.sqrt(max(or_use * (or_use - 1.0), 0.0)))


def main():
    meta = _load_meta()
    if meta is None or meta.empty:
        pd.DataFrame().to_csv(OUT / "causal_marker_effects.tsv", sep="\t", index=False)
        with open(OUT / "causal_summary.json", "w", encoding="utf-8") as f:
            json.dump({"note": "No sample metadata available for causal inference."}, f, ensure_ascii=False, indent=2)
        print("No metadata; causal output empty.")
        return

    met_mat = _load_matrix(Path("outputs/filtered_matrix.tsv"))
    prot_mat = _load_matrix(Path("outputs_prot/filtered_prot_matrix.tsv"))
    cand = _load_top_candidates()
    Z = _build_confounders(meta, met_mat, prot_mat, n_pc=CAUSAL_CFG["n_pc"])
    y = meta["y"].values.astype(float)
    samples = meta.index.tolist()

    records = []
    for _, r in cand.iterrows():
        mod = r["modality"]
        feat = str(r["feature_id"])
        marker = str(r["marker"])
        mat = prot_mat if mod == "protein" else met_mat
        if mat is None or feat not in mat.index:
            continue
        common = [s for s in samples if s in mat.columns]
        if len(common) < CAUSAL_CFG["min_samples_per_marker"]:
            continue
        t = pd.to_numeric(mat.loc[feat, common], errors="coerce").fillna(0.0).values.astype(float)
        yy = meta.loc[common, "y"].values.astype(float)
        zz = Z.loc[common].values.astype(float) if not Z.empty else np.zeros((len(common), 1))

        theta = _dml_plr_binary_y(
            yy,
            t,
            zz,
            n_splits=CAUSAL_CFG["n_splits_dml"],
            seed=CAUSAL_CFG["seed"],
            rf_n_estimators=CAUSAL_CFG["rf_n_estimators"],
        )
        lo, hi = _bootstrap_ci(
            yy,
            t,
            zz,
            n_boot=CAUSAL_CFG["n_bootstrap"],
            seed=CAUSAL_CFG["seed"],
            n_splits=CAUSAL_CFG["n_splits_dml"],
            rf_n_estimators=CAUSAL_CFG["rf_n_estimators"],
        )

        # Logistic OR for intuitive biomedical interpretation
        try:
            X_lr = pd.DataFrame(zz)
            X_lr["t"] = t
            lr = LogisticRegressionCV(
                cv=5,
                l1_ratios=(0.0,),
                max_iter=5000,
                scoring="roc_auc",
                class_weight="balanced",
                use_legacy_attributes=True,
            )
            lr.fit(X_lr, yy.astype(int))
            beta_t = float(lr.coef_.ravel()[-1])
            or_t = float(np.exp(beta_t))
        except Exception:
            or_t = np.nan
        e_val = _e_value_from_or(or_t)

        records.append(
            {
                "modality": mod,
                "feature_id": feat,
                "marker": marker,
                "dml_effect_theta": theta,
                "ci95_low": lo,
                "ci95_high": hi,
                "or_adjusted": or_t,
                "e_value": e_val,
                "n_samples": len(common),
            }
        )

    out = pd.DataFrame(records)
    if out.empty:
        out.to_csv(OUT / "causal_marker_effects.tsv", sep="\t", index=False)
        with open(OUT / "causal_summary.json", "w", encoding="utf-8") as f:
            json.dump({"note": "No candidate marker had enough data for DML."}, f, ensure_ascii=False, indent=2)
        print("No valid causal effects estimated.")
        return

    out["abs_dml"] = out["dml_effect_theta"].abs()
    out = out.sort_values("abs_dml", ascending=False).reset_index(drop=True)
    out.to_csv(OUT / "causal_marker_effects.tsv", sep="\t", index=False)

    summary = {
        "method": "Double Machine Learning (PLR) + bootstrap CI + adjusted OR + E-value",
        "config": CAUSAL_CFG,
        "n_markers": int(out.shape[0]),
        "top_markers": out.head(10)[["modality", "marker", "dml_effect_theta", "or_adjusted", "e_value"]].to_dict(orient="records"),
    }
    with open(OUT / "causal_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)
    print("Causal inference outputs written to", OUT)


if __name__ == "__main__":
    main()
