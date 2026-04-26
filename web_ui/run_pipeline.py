import json
import os
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd

OUTPUT_DIRS = ["outputs", "outputs_advanced", "outputs_prot", "outputs_multiomics"]
PIPELINE_SCRIPTS = [
    "preprocess_and_analyze.py",
    "proteomics_analysis.py",
    "multiomics_integration.py",
    "generate_reports.py",
]


def _read_tabular(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix in {".xlsx", ".xls"}:
        return pd.read_excel(path)
    if suffix == ".tsv":
        return pd.read_csv(path, sep="\t")
    return pd.read_csv(path)


def _normalize_to_xlsx(src: str, dst: Path) -> Path:
    src_path = Path(src)
    dst.parent.mkdir(parents=True, exist_ok=True)
    if src_path.suffix.lower() in {".xlsx", ".xls"}:
        shutil.copy(src_path, dst)
        return dst
    df = _read_tabular(src_path)
    df.to_excel(dst, index=False)
    return dst


def prepare_run_dirs(base_outdir: str) -> Path:
    root = Path(base_outdir).resolve()
    run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = root / run_id
    (run_dir / "inputs").mkdir(parents=True, exist_ok=True)
    (run_dir / "artifacts").mkdir(parents=True, exist_ok=True)
    return run_dir


def run_pipeline(protein_path: str, metabolite_path: str, outdir: str = "outputs_web"):
    repo_root = Path(__file__).resolve().parent.parent
    run_dir = prepare_run_dirs(outdir)
    inputs_dir = run_dir / "inputs"
    artifacts_dir = run_dir / "artifacts"

    yield f"[INFO] 本次任务目录：{run_dir}\n"
    prot_norm = _normalize_to_xlsx(protein_path, inputs_dir / "protein_input.xlsx")
    met_norm = _normalize_to_xlsx(metabolite_path, inputs_dir / "metabolite_input.xlsx")
    yield "[INFO] 输入文件已规范化为 xlsx。\n"

    target_prot = repo_root / "prot.xlsx"
    target_met = repo_root / "metab.xlsx"
    shutil.copy(prot_norm, target_prot)
    shutil.copy(met_norm, target_met)
    yield "[INFO] 已写入管线标准输入：prot.xlsx / metab.xlsx\n"

    for script in PIPELINE_SCRIPTS:
        cmd = [sys.executable, script]
        yield f"[RUN ] {' '.join(cmd)}\n"
        proc = subprocess.Popen(
            cmd,
            cwd=repo_root,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        assert proc.stdout is not None
        for line in proc.stdout:
            yield line
        proc.wait()
        if proc.returncode != 0:
            raise RuntimeError(f"{script} 运行失败，退出码={proc.returncode}")
        yield f"[DONE] {script}\n"

    summary = {}
    for output_name in OUTPUT_DIRS:
        src = repo_root / output_name
        if src.exists():
            dst = artifacts_dir / output_name
            if dst.exists():
                shutil.rmtree(dst)
            shutil.copytree(src, dst)
            yield f"[INFO] 已归档：{output_name}\n"

    summary_path = artifacts_dir / "outputs" / "summary.json"
    if summary_path.exists():
        try:
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
        except Exception:
            summary = {}

    result = {
        "run_dir": str(run_dir),
        "artifacts_dir": str(artifacts_dir),
        "summary": summary,
    }
    result_file = run_dir / "run_result.json"
    result_file.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    yield f"[INFO] 结果索引已写入：{result_file}\n"
    return result
