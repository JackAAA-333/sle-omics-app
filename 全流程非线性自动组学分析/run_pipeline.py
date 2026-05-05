import json
import re
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd

OUTPUT_DIRS = [
    "outputs",
    "outputs_advanced",
    "outputs_prot",
    "outputs_multiomics",
    "outputs_literature",
    "outputs_causal",
    "outputs_candidate",
]
LOGO_PATH = Path(
    "/Users/jacka/.cursor/projects/Users-jacka-Desktop-SLE/assets/__logo__-___-40948cae-d905-4025-b378-50bdd1b3110b.png"
)
ORG_NAME = "SLE 多组学联合分析项目组"
PROJECT_NAME = "全流程非线性自动组学分析"
HEADER_KB_PATH = Path(__file__).resolve().parent.parent / "header_alias_registry.json"


def _read_tabular(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix in {".xlsx", ".xls"}:
        return pd.read_excel(path)
    if suffix == ".tsv":
        return pd.read_csv(path, sep="\t")
    return pd.read_csv(path)


def _norm_col_key(text: str) -> str:
    return re.sub(r"[^a-z0-9\u4e00-\u9fff]+", "", str(text).strip().lower())


def _find_column_by_alias(df: pd.DataFrame, aliases: list[str]) -> str | None:
    alias_set = {_norm_col_key(a) for a in aliases}
    for col in df.columns:
        if _norm_col_key(col) in alias_set:
            return str(col)
    return None


def _load_header_kb() -> dict:
    if not HEADER_KB_PATH.exists():
        return {}
    try:
        data = json.loads(HEADER_KB_PATH.read_text(encoding="utf-8"))
        return data if isinstance(data, dict) else {}
    except Exception:
        return {}


def _save_header_kb(kb: dict):
    try:
        HEADER_KB_PATH.write_text(json.dumps(kb, ensure_ascii=False, indent=2), encoding="utf-8")
    except Exception:
        pass


def _kb_aliases(kb: dict, modality: str, standard_col: str) -> list[str]:
    items = (((kb.get(modality, {}) or {}).get("column_aliases", {}) or {}).get(standard_col, []) or [])
    return [str(x) for x in items if str(x).strip()]


def _kb_sample_aliases(kb: dict, group: str) -> list[str]:
    items = ((kb.get("sample_group_aliases", {}) or {}).get(group, []) or [])
    return [str(x) for x in items if str(x).strip()]


def _record_kb_alias(kb: dict, modality: str, standard_col: str, raw_col: str):
    kb.setdefault(modality, {}).setdefault("column_aliases", {}).setdefault(standard_col, [])
    arr = kb[modality]["column_aliases"][standard_col]
    if raw_col not in arr:
        arr.append(raw_col)


def _record_kb_sample_alias(kb: dict, group: str, raw_col: str):
    kb.setdefault("sample_group_aliases", {}).setdefault(group, [])
    arr = kb["sample_group_aliases"][group]
    if raw_col not in arr:
        arr.append(raw_col)


def _classify_sample_group(col_name: str, kb: dict | None = None) -> str | None:
    key = _norm_col_key(col_name)
    sle_tokens = ["sle", "lupus", "case", "patient", "disease", "病例", "患者", "实验组", "病组"] + _kb_sample_aliases(kb or {}, "SLE")
    hc_tokens = ["hc", "ctrl", "control", "healthy", "normal", "对照", "健康", "正常", "阴性"] + _kb_sample_aliases(kb or {}, "HC")
    qc_tokens = ["qc", "pool", "pooled", "reference", "质控", "质控样"] + _kb_sample_aliases(kb or {}, "QC")
    if any(t in key for t in qc_tokens):
        return "QC"
    if any(t in key for t in sle_tokens):
        return "SLE"
    if any(t in key for t in hc_tokens):
        return "HC"
    m = re.search(r"(sle|hc|qc)", key, flags=re.IGNORECASE)
    if m:
        return m.group(1).upper()
    return None


def _normalize_sample_headers(df: pd.DataFrame, reserved_cols: set[str] | None = None, kb: dict | None = None) -> tuple[pd.DataFrame, list[dict]]:
    df = df.copy()
    reserved_cols = reserved_cols or set()
    renamed = {}
    mapping = []
    seen = set()
    auto_counter = {"SLE": 1, "HC": 1, "QC": 1}

    for col in df.columns:
        col_str = str(col).strip()
        if str(col) in reserved_cols:
            renamed[col] = col
            continue
        group = _classify_sample_group(col_str, kb=kb)
        if not group:
            renamed[col] = col
            continue
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
        mapping.append({"raw": col_str, "standard": candidate, "type": "sample_group", "group": group})
        if kb is not None:
            _record_kb_sample_alias(kb, group, col_str)

    df.rename(columns=renamed, inplace=True)
    return df, mapping


def _normalize_metabolomics_headers(df: pd.DataFrame, kb: dict | None = None) -> tuple[pd.DataFrame, list[dict]]:
    df = df.copy()
    mapping = []
    id_col = _find_column_by_alias(
        df,
        [
            "Met ID",
            "met_id",
            "metabolite id",
            "metabolite_id",
            "feature id",
            "feature_id",
            "compound id",
            "compound_id",
            "id",
            "代谢物ID",
            "代谢物编号",
        ] + _kb_aliases(kb or {}, "metabolite", "Met ID"),
    )
    name_col = _find_column_by_alias(
        df,
        [
            "Metabolite Name",
            "metabolite_name",
            "metabolite",
            "compound name",
            "compound_name",
            "name",
            "代谢物名称",
            "代谢物",
            "名称",
        ] + _kb_aliases(kb or {}, "metabolite", "Metabolite Name"),
    )
    if id_col is None:
        id_col = str(df.columns[0])
    rename_map = {}
    if id_col != "Met ID":
        rename_map[id_col] = "Met ID"
        mapping.append({"raw": str(id_col), "standard": "Met ID", "type": "feature_id"})
        if kb is not None:
            _record_kb_alias(kb, "metabolite", "Met ID", str(id_col))
    if name_col and name_col != "Metabolite Name":
        rename_map[name_col] = "Metabolite Name"
        mapping.append({"raw": str(name_col), "standard": "Metabolite Name", "type": "feature_name"})
        if kb is not None:
            _record_kb_alias(kb, "metabolite", "Metabolite Name", str(name_col))
    if rename_map:
        df.rename(columns=rename_map, inplace=True)
    reserved = {"Met ID", "Metabolite Name"}
    df, sample_mapping = _normalize_sample_headers(df, reserved_cols=reserved, kb=kb)
    mapping.extend(sample_mapping)
    return df, mapping


def _normalize_proteomics_headers(df: pd.DataFrame, kb: dict | None = None) -> tuple[pd.DataFrame, list[dict]]:
    df = df.copy()
    mapping = []
    id_col = _find_column_by_alias(
        df,
        [
            "Protein ID",
            "protein_id",
            "Protein IDs",
            "Protein.Group",
            "Accession",
            "UniProt",
            "Uniprot ID",
            "Entry",
            "蛋白ID",
            "蛋白编号",
        ] + _kb_aliases(kb or {}, "protein", "Protein ID"),
    )
    if id_col is None:
        id_col = str(df.columns[0])
    gene_col = _find_column_by_alias(
        df,
        [
            "Gene Symbol",
            "Gene",
            "Gene Name",
            "Gene Names",
            "symbol",
            "基因名",
            "基因符号",
        ] + _kb_aliases(kb or {}, "protein", "Gene Symbol"),
    )
    rename_map = {}
    if id_col != "Protein ID":
        rename_map[id_col] = "Protein ID"
        mapping.append({"raw": str(id_col), "standard": "Protein ID", "type": "feature_id"})
        if kb is not None:
            _record_kb_alias(kb, "protein", "Protein ID", str(id_col))
    if gene_col and gene_col != "Gene Symbol":
        rename_map[gene_col] = "Gene Symbol"
        mapping.append({"raw": str(gene_col), "standard": "Gene Symbol", "type": "feature_name"})
        if kb is not None:
            _record_kb_alias(kb, "protein", "Gene Symbol", str(gene_col))
    if rename_map:
        df.rename(columns=rename_map, inplace=True)
    reserved = {"Protein ID", "Gene Symbol"}
    df, sample_mapping = _normalize_sample_headers(df, reserved_cols=reserved, kb=kb)
    mapping.extend(sample_mapping)
    return df, mapping


def _normalize_to_xlsx(src: str, dst: Path, modality: str, kb: dict | None = None) -> tuple[Path, list[dict]]:
    src_path = Path(src)
    dst.parent.mkdir(parents=True, exist_ok=True)
    if src_path.suffix.lower() in {".xlsx", ".xls"}:
        df = pd.read_excel(src_path)
        if modality == "protein":
            df, mapping = _normalize_proteomics_headers(df, kb=kb)
        elif modality == "metabolite":
            df, mapping = _normalize_metabolomics_headers(df, kb=kb)
        else:
            df, mapping = _normalize_sample_headers(df, kb=kb)
        df.to_excel(dst, index=False)
        return dst, mapping
    df = _read_tabular(src_path)
    if modality == "protein":
        df, mapping = _normalize_proteomics_headers(df, kb=kb)
    elif modality == "metabolite":
        df, mapping = _normalize_metabolomics_headers(df, kb=kb)
    else:
        df, mapping = _normalize_sample_headers(df, kb=kb)
    df.to_excel(dst, index=False)
    return dst, mapping


def _write_header_mapping_report(run_dir: Path, mappings: dict):
    report_data = {
        "generated_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "mappings": mappings,
    }
    json_path = run_dir / "header_mapping_audit.json"
    json_path.write_text(json.dumps(report_data, ensure_ascii=False, indent=2), encoding="utf-8")
    md_path = run_dir / "header_mapping_audit.md"
    lines = [
        "# 表头识别与标准化审计报告",
        "",
        f"- 生成时间：{report_data['generated_at']}",
        "",
    ]
    for modality in ["protein", "metabolite"]:
        rows = mappings.get(modality, [])
        lines.append(f"## {modality} 映射")
        if not rows:
            lines.append("无映射记录。")
            lines.append("")
            continue
        lines.append("| 原始表头 | 标准表头 | 类型 |")
        lines.append("|---|---|---|")
        for r in rows:
            lines.append(f"| {r.get('raw','')} | {r.get('standard','')} | {r.get('type','')} |")
        lines.append("")
    md_path.write_text("\n".join(lines), encoding="utf-8")
    return json_path, md_path


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
    scripts.append("causal_inference_latest.py")
    scripts.append("biomarker_candidate_report.py")
    scripts.append("generate_reports.py")
    return scripts


def run_pipeline(
    protein_path: str | None = None,
    metabolite_path: str | None = None,
    outdir: str = "outputs_web",
    run_name: str | None = None,
    strict_profile: str = "标准",
    strict_min_abs_dml: float | None = None,
    strict_min_e_value: float | None = None,
):
    repo_root = Path(__file__).resolve().parent.parent
    has_protein = bool(protein_path)
    has_metabolite = bool(metabolite_path)
    if not has_protein and not has_metabolite:
        raise ValueError("至少需要上传一类组学数据（蛋白组或代谢组）。")

    run_dir = prepare_run_dirs(outdir, run_name=run_name)
    inputs_dir = run_dir / "inputs"
    artifacts_dir = run_dir / "artifacts"
    kb = _load_header_kb()
    header_mappings = {"protein": [], "metabolite": []}

    yield f"[INFO] 本次任务目录：{run_dir}\n"
    if has_protein:
        prot_norm, prot_map = _normalize_to_xlsx(protein_path, inputs_dir / "protein_input.xlsx", modality="protein", kb=kb)
        header_mappings["protein"] = prot_map
        target_prot = repo_root / "prot.xlsx"
        shutil.copy(prot_norm, target_prot)
        yield "[INFO] 蛋白组输入已规范化并写入 prot.xlsx。\n"
        if prot_map:
            for m in prot_map:
                yield f"[MAP ] protein: {m.get('raw')} -> {m.get('standard')} ({m.get('type')})\n"
    if has_metabolite:
        met_norm, met_map = _normalize_to_xlsx(metabolite_path, inputs_dir / "metabolite_input.xlsx", modality="metabolite", kb=kb)
        header_mappings["metabolite"] = met_map
        target_met = repo_root / "metab.xlsx"
        shutil.copy(met_norm, target_met)
        yield "[INFO] 代谢组输入已规范化并写入 metab.xlsx。\n"
        if met_map:
            for m in met_map:
                yield f"[MAP ] metabolite: {m.get('raw')} -> {m.get('standard')} ({m.get('type')})\n"
    _save_header_kb(kb)
    _, audit_md = _write_header_mapping_report(run_dir, header_mappings)
    yield f"[INFO] 表头映射审计报告已生成：{audit_md.name}\n"

    _clean_previous_outputs(repo_root)
    scripts = _select_pipeline_scripts(has_protein=has_protein, has_metabolite=has_metabolite)
    yield f"[INFO] 本次将执行脚本：{', '.join(scripts)}\n"

    for script in scripts:
        cmd = [sys.executable, script]
        if script == "biomarker_candidate_report.py":
            cmd.extend(["--strict-profile", str(strict_profile or "标准")])
            if strict_min_abs_dml is not None:
                cmd.extend(["--min-abs-dml", str(float(strict_min_abs_dml))])
            if strict_min_e_value is not None:
                cmd.extend(["--min-e-value", str(float(strict_min_e_value))])
        yield f"[RUN ] {' '.join(cmd)}\n"
        proc = subprocess.Popen(
            cmd,
            cwd=repo_root,
            stdin=subprocess.DEVNULL,
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
    for p in [run_dir / "header_mapping_audit.json", run_dir / "header_mapping_audit.md"]:
        if p.exists():
            shutil.copy(p, artifacts_dir / p.name)
            yield f"[INFO] 已归档：{p.name}\n"

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
