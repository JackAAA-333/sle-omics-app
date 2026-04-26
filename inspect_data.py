import sys
import json
import pandas as pd

FILES = ["prot.xlsx", "metab.xlsx"]

def summarize(path):
    info = {"file": path}
    try:
        df = pd.read_excel(path)
    except Exception as e:
        info["error"] = str(e)
        return info

    info["shape"] = df.shape
    info["columns"] = list(df.columns)
    info["dtypes"] = {c: str(df[c].dtype) for c in df.columns}
    info["missing_counts"] = {c: int(df[c].isna().sum()) for c in df.columns}

    # look for group columns containing sle/hc
    groups = {}
    lowered = None
    for c in df.columns:
        try:
            vals = df[c].dropna().astype(str).str.lower().unique().tolist()
        except Exception:
            vals = []
        if any(v in ("sle","hc","control","case","case","healthy","healthy control","healthy_control") for v in vals):
            groups[c] = vals[:20]
    info["group_candidates"] = groups

    # if samples are columns (features in first column), try transpose detection
    if df.shape[0] < df.shape[1]:
        info["note"] = "rows < columns — possibly features in rows and samples in columns"

    return info

def main():
    results = []
    for f in FILES:
        results.append(summarize(f))

    out = {"results": results}
    print(json.dumps(out, indent=2, ensure_ascii=False))
    with open("inspect_summary.json", "w", encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=2)

if __name__ == '__main__':
    main()
