import json
import os
from datetime import datetime
from pathlib import Path

import pandas as pd
import requests

OUT = Path("outputs_literature")
OUT.mkdir(parents=True, exist_ok=True)


def _pubmed_recent3y(marker: str, retmax: int = 10):
    end_year = datetime.now().year
    start_year = end_year - 2
    term = f'("{marker}") AND ("systemic lupus erythematosus" OR SLE)'
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {
        "db": "pubmed",
        "retmode": "json",
        "retmax": retmax,
        "term": term,
        "mindate": f"{start_year}/01/01",
        "maxdate": f"{end_year}/12/31",
        "datetype": "pdat",
    }
    try:
        r = requests.get(url, params=params, timeout=12)
        r.raise_for_status()
        j = r.json().get("esearchresult", {})
        count = int(j.get("count", 0))
        pmids = j.get("idlist", [])
        return {
            "query": term,
            "recent_3y_count": count,
            "pmids": pmids,
        }
    except Exception:
        return {
            "query": term,
            "recent_3y_count": 0,
            "pmids": [],
        }


def _collect_protein_markers():
    p = Path("outputs_prot/xgb_shap_prot_annotated.tsv")
    if not p.exists():
        return []
    try:
        df = pd.read_csv(p, sep="\t")
        rows = []
        for _, row in df.head(20).iterrows():
            protein_id = str(row.get("protein_id", "")).strip()
            gene_symbol = str(row.get("gene_symbol", "")).strip()
            marker = gene_symbol if gene_symbol and gene_symbol.lower() != "nan" else protein_id
            if marker:
                rows.append(("protein", marker, protein_id, gene_symbol))
        return rows
    except Exception:
        return []


def _collect_metabolite_markers():
    p = Path("outputs_advanced/xgb_shap_feature_importance.tsv")
    if not p.exists():
        return []
    try:
        imp = pd.read_csv(p, sep="\t", index_col=0)
        markers = [str(x) for x in imp.head(20).index]
    except Exception:
        return []

    name_map = {}
    metab_path = Path("metab.xlsx")
    if metab_path.exists():
        try:
            mdf = pd.read_excel(metab_path)
            if "Met ID" in mdf.columns and "Metabolite Name" in mdf.columns:
                for _, row in mdf.iterrows():
                    key = str(row.get("Met ID", "")).strip()
                    val = str(row.get("Metabolite Name", "")).strip()
                    if key and val and val.lower() != "nan":
                        name_map[key] = val
        except Exception:
            pass

    rows = []
    for met_id in markers:
        marker = name_map.get(met_id, met_id)
        rows.append(("metabolite", marker, met_id, ""))
    return rows


def main():
    marker_rows = _collect_protein_markers() + _collect_metabolite_markers()
    if len(marker_rows) == 0:
        empty = pd.DataFrame(
            columns=[
                "modality",
                "marker",
                "feature_id",
                "gene_symbol",
                "query",
                "recent_3y_count",
                "pmids",
            ]
        )
        empty.to_csv(OUT / "literature_meta_recent3y.tsv", sep="\t", index=False)
        with open(OUT / "literature_meta_recent3y_summary.json", "w", encoding="utf-8") as f:
            json.dump({"note": "No markers available for literature meta summary."}, f, ensure_ascii=False, indent=2)
        print("No markers found; literature output created with empty content.")
        return

    records = []
    for modality, marker, feat_id, gene_symbol in marker_rows:
        rs = _pubmed_recent3y(marker)
        records.append(
            {
                "modality": modality,
                "marker": marker,
                "feature_id": feat_id,
                "gene_symbol": gene_symbol,
                "query": rs["query"],
                "recent_3y_count": rs["recent_3y_count"],
                "pmids": ",".join(rs["pmids"]),
            }
        )

    out_df = pd.DataFrame(records).sort_values(["recent_3y_count", "modality"], ascending=[False, True])
    out_df.to_csv(OUT / "literature_meta_recent3y.tsv", sep="\t", index=False)

    summary = {
        "n_markers": int(out_df.shape[0]),
        "n_markers_with_recent_evidence": int((out_df["recent_3y_count"] > 0).sum()),
        "total_recent_3y_hits": int(out_df["recent_3y_count"].sum()),
        "top_markers": out_df.head(10)[["modality", "marker", "recent_3y_count"]].to_dict(orient="records"),
    }
    with open(OUT / "literature_meta_recent3y_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)

    print("Literature meta summary generated in", OUT)


if __name__ == "__main__":
    main()
