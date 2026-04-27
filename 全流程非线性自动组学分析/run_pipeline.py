import json
import re
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd

OUTPUT_DIRS = ["outputs", "outputs_advanced", "outputs_prot", "outputs_multiomics", "outputs_literature"]
LOGO_PATH = Path(
    "/Users/jacka/.cursor/projects/Users-jacka-Desktop-SLE/assets/__logo__-___-40948cae-d905-4025-b378-50bdd1b3110b.png"
)
ORG_NAME = "SLE 多组学联合分析项目组"
PROJECT_NAME = "全流程非线性自动组学分析"


def _read_tabular(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix in {".xlsx", ".xls"}:
        return pd.read_excel(path)
    if suffix == ".tsv":
        return pd.read_csv(path, sep="\t")
    return pd.read_csv(path)


def _normalize_sample_headers(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    renamed = {}
    seen = set()
    auto_counter = {"SLE": 1, "HC": 1, "QC": 1}
    pattern = re.compile(r"(SLE|HC|QC)", re.IGNORECASE)

    for col in df.columns:
        col_str = str(col).strip()
        m = pattern.search(col_str.upper())
        if not m:
            renamed[col] = col
            continue
        group = m.group(1).upper()
        digits = re.findall(r"\d+", col_str)
        idx = digits[-1] if digits else str(auto_counter[group])
        if not digits:
            auto_counter[group] += 1
        candidate = f"{group}-{idx}"
        if candidate in seen:
            suffix = 2
            while f"{candidate}_{suffix}" in seen:
                suffix += 1
            candidate = f"{candidate}_{suffix}"
        seen.add(candidate)
        renamed[col] = candidate

    df.rename(columns=renamed, inplace=True)
    return df


def _normalize_to_xlsx(src: str, dst: Path) -> Path:
    src_path = Path(src)
    dst.parent.mkdir(parents=True, exist_ok=True)
    if src_path.suffix.lower() in {".xlsx", ".xls"}:
        df = pd.read_excel(src_path)
        df = _normalize_sample_headers(df)
        df.to_excel(dst, index=False)
        return dst
    df = _read_tabular(src_path)
    df = _normalize_sample_headers(df)
    df.to_excel(dst, index=False)
    return dst


def prepare_run_dirs(base_outdir: str, run_name: str | None = None) -> Path:
    root = Path(base_outdir).resolve()
    run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    if run_name:
        safe_name = re.sub(r"[^\w\-\u4e00-\u9fff]", "_", run_name.strip())
        run_id = safe_name or run_id
        target = root / run_id
        if target.exists():
            run_id = f"{run_id}_{datetime.now().strftime('%H%M%S')}"
    run_dir = root / run_id
    (run_dir / "inputs").mkdir(parents=True, exist_ok=True)
    (run_dir / "artifacts").mkdir(parents=True, exist_ok=True)
    return run_dir


def _build_report_cover(artifacts_dir: Path) -> Path:
    branding_dir = artifacts_dir / "branding"
    branding_dir.mkdir(parents=True, exist_ok=True)
    logo_dst = branding_dir / "logo.png"
    if LOGO_PATH.exists():
        shutil.copy(LOGO_PATH, logo_dst)

    cover_path = artifacts_dir / "report_cover.md"
    cover_text = (
        f"# {PROJECT_NAME} 报告封面\n\n"
        f"![Logo](branding/logo.png)\n\n"
        f"- 项目名称：{PROJECT_NAME}\n"
        f"- 署名单位：{ORG_NAME}\n"
        f"- 生成时间：{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
    )
    cover_path.write_text(cover_text, encoding="utf-8")
    return cover_path


def _clean_previous_outputs(repo_root: Path):
    for output_name in OUTPUT_DIRS:
        target = repo_root / output_name
        if target.exists():
            shutil.rmtree(target)


def _select_pipeline_scripts(has_protein: bool, has_metabolite: bool):
    scripts = []
    if has_metabolite:
        scripts.append("preprocess_and_analyze.py")
    if has_protein:
        scripts.append("proteomics_analysis.py")
    if has_metabolite and has_protein:
        scripts.append("multiomics_integration.py")
    scripts.append("literature_meta_recent3y.py")
    if has_metabolite:
        scripts.append("generate_reports.py")
    return scripts


def run_pipeline(
    protein_path: str | None = None,
    metabolite_path: str | None = None,
    outdir: str = "outputs_web",
    run_name: str | None = None,
):
    repo_root = Path(__file__).resolve().parent.parent
    has_protein = bool(protein_path)
    has_metabolite = bool(metabolite_path)
    if not has_protein and not has_metabolite:
        raise ValueError("至少需要上传一类组学数据（蛋白组或代谢组）。")

    run_dir = prepare_run_dirs(outdir, run_name=run_name)
    inputs_dir = run_dir / "inputs"
    artifacts_dir = run_dir / "artifacts"

    yield f"[INFO] 本次任务目录：{run_dir}\n"
    if has_protein:
        prot_norm = _normalize_to_xlsx(protein_path, inputs_dir / "protein_input.xlsx")
        target_prot = repo_root / "prot.xlsx"
        shutil.copy(prot_norm, target_prot)
        yield "[INFO] 蛋白组输入已规范化并写入 prot.xlsx。\n"
    if has_metabolite:
        met_norm = _normalize_to_xlsx(metabolite_path, inputs_dir / "metabolite_input.xlsx")
        target_met = repo_root / "metab.xlsx"
        shutil.copy(met_norm, target_met)
        yield "[INFO] 代谢组输入已规范化并写入 metab.xlsx。\n"

    _clean_previous_outputs(repo_root)
    scripts = _select_pipeline_scripts(has_protein=has_protein, has_metabolite=has_metabolite)
    yield f"[INFO] 本次将执行脚本：{', '.join(scripts)}\n"

    for script in scripts:
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

    cover_file = _build_report_cover(artifacts_dir)
    yield f"[INFO] 报告封面已生成：{cover_file.name}\n"

    summary_candidates = [
        artifacts_dir / "outputs" / "summary.json",
        artifacts_dir / "outputs_prot" / "summary.json",
    ]
    for summary_path in summary_candidates:
        if summary_path.exists():
            try:
                summary = json.loads(summary_path.read_text(encoding="utf-8"))
                break
            except Exception:
                summary = {}

    result = {
        "run_dir": str(run_dir),
        "artifacts_dir": str(artifacts_dir),
        "summary": summary,
        "has_protein": has_protein,
        "has_metabolite": has_metabolite,
    }
    result_file = run_dir / "run_result.json"
    result_file.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    yield f"[INFO] 结果索引已写入：{result_file}\n"
    return result
