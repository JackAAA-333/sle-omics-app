import zipfile
from pathlib import Path

import streamlit as st

import run_pipeline

LOGO_PATH = Path(
    "/Users/jacka/.cursor/projects/Users-jacka-Desktop-SLE/assets/__logo__-___-40948cae-d905-4025-b378-50bdd1b3110b.png"
)
ORG_NAME = "SLE 多组学联合分析项目组"
PROJECT_NAME = "全流程非线性自动组学分析"

st.set_page_config(page_title=PROJECT_NAME, page_icon=str(LOGO_PATH) if LOGO_PATH.exists() else "🧬", layout="wide")


def save_upload(upload, target_dir: Path, name_hint: str) -> str:
    target_dir.mkdir(parents=True, exist_ok=True)
    dest = target_dir / name_hint
    with open(dest, "wb") as f:
        f.write(upload.getbuffer())
    return str(dest)


def collect_files(root: Path):
    return [p for p in sorted(root.rglob("*")) if p.is_file()]


def make_zip(source_dir: Path, zip_path: Path):
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for file in collect_files(source_dir):
            zf.write(file, arcname=str(file.relative_to(source_dir)))


header_left, header_right = st.columns([1, 7])
with header_left:
    if LOGO_PATH.exists():
        st.image(str(LOGO_PATH), width=120)
with header_right:
    st.title(PROJECT_NAME)
st.caption("上传蛋白组学与代谢组学数据，自动执行全流程非线性整合分析并输出专业结果。")

left, right = st.columns(2)
with left:
    protein_file = st.file_uploader("蛋白组数据", type=["csv", "tsv", "xlsx"], key="prot")
with right:
    metabolite_file = st.file_uploader("代谢组数据", type=["csv", "tsv", "xlsx"], key="met")

outdir = st.text_input("结果保存目录（父目录）", value="outputs_web")
run_name = st.text_input("结果子文件夹名（可选，自定义本次任务文件夹）", value="")
run_button = st.button("开始自动分析", type="primary", use_container_width=True)

mode_title = ""
mode_desc = ""
planned_modules = []
mode_color = "#6b7280"
if protein_file and metabolite_file:
    mode_title = "双组学模式"
    mode_desc = "已上传蛋白组 + 代谢组，将执行单组学分析与多组学整合分析。"
    mode_color = "#7c3aed"
    planned_modules = [
        "preprocess_and_analyze.py（代谢组分析）",
        "proteomics_analysis.py（蛋白组分析）",
        "multiomics_integration.py（多组学整合）",
        "generate_reports.py（汇总报告）",
    ]
elif protein_file:
    mode_title = "单组学模式（蛋白组）"
    mode_desc = "仅上传蛋白组，将执行蛋白组分析。"
    mode_color = "#2563eb"
    planned_modules = [
        "proteomics_analysis.py（蛋白组分析）",
    ]
elif metabolite_file:
    mode_title = "单组学模式（代谢组）"
    mode_desc = "仅上传代谢组，将执行代谢组分析与报告生成。"
    mode_color = "#0ea5e9"
    planned_modules = [
        "preprocess_and_analyze.py（代谢组分析）",
        "generate_reports.py（汇总报告）",
    ]
else:
    mode_title = "待上传数据"
    mode_desc = "请至少上传一类组学数据（蛋白组或代谢组）。"
    mode_color = "#6b7280"

st.markdown("### 当前模式")
st.markdown(
    f"<span style='display:inline-block;padding:4px 10px;border-radius:9999px;"
    f"background:{mode_color};color:white;font-size:0.9rem;font-weight:600;'>{mode_title}</span>",
    unsafe_allow_html=True,
)
if protein_file or metabolite_file:
    st.success(mode_desc)
else:
    st.info(mode_desc)
if planned_modules:
    st.markdown("**本次将执行模块：**")
    for module in planned_modules:
        st.markdown(f"- {module}")

st.markdown("### 流程说明")
st.markdown(
    "- 自动数据规范化与样本表头识别\n"
    "- 支持只上传单组学（蛋白或代谢）也可分析\n"
    "- 双组学时自动执行非线性整合建模（XGBoost + SHAP）\n"
    "- 自动汇总报告、图像预览与结果打包"
)

if run_button:
    if not protein_file and not metabolite_file:
        st.error("请至少上传一类组学文件（蛋白组或代谢组）。")
    else:
        staging = Path(outdir) / "_uploads"
        saved_prot = None
        saved_met = None
        if protein_file:
            saved_prot = save_upload(protein_file, staging, f"protein{Path(protein_file.name).suffix.lower()}")
        if metabolite_file:
            saved_met = save_upload(metabolite_file, staging, f"metabolite{Path(metabolite_file.name).suffix.lower()}")

        st.info("任务已启动，正在执行完整分析流程。")
        log_area = st.empty()
        logs = []
        result = {}

        generator = run_pipeline.run_pipeline(saved_prot, saved_met, outdir=outdir, run_name=run_name.strip() or None)
        try:
            while True:
                line = next(generator)
                logs.append(line)
                log_area.text_area("运行日志", value="".join(logs), height=360)
        except StopIteration as stop:
            result = stop.value or {}
        except Exception as exc:
            st.error(f"分析失败：{exc}")

        if result:
            st.success("分析完成，可查看关键指标和下载完整结果。")
            run_dir = Path(result["run_dir"])
            artifacts_dir = Path(result["artifacts_dir"])
            summary = result.get("summary", {})

            m1, m2, m3 = st.columns(3)
            m1.metric("样本数", summary.get("n_samples", "N/A"))
            m2.metric("初始特征数", summary.get("n_features", "N/A"))
            m3.metric("FDR<0.05 特征", summary.get("n_significant_fdr05", "N/A"))

            tabs = st.tabs(["图像预览", "结果下载"])
            with tabs[0]:
                image_candidates = [
                    artifacts_dir / "outputs" / "volcano.png",
                    artifacts_dir / "outputs" / "heatmap_top30.png",
                    artifacts_dir / "outputs" / "roc_top5.png",
                ]
                shown = 0
                for img in image_candidates:
                    if img.exists():
                        st.image(str(img), caption=str(img.relative_to(artifacts_dir)), use_container_width=True)
                        shown += 1
                if shown == 0:
                    st.warning("本次运行未检测到可预览图像，可在下载区获取全部文件。")

            with tabs[1]:
                zip_path = run_dir / "analysis_results.zip"
                cover_path = artifacts_dir / "report_cover.md"
                make_zip(artifacts_dir, zip_path)
                st.download_button(
                    "下载全部结果 (ZIP)",
                    data=zip_path.read_bytes(),
                    file_name=zip_path.name,
                    mime="application/zip",
                    use_container_width=True,
                )
                if cover_path.exists():
                    st.download_button(
                        "下载报告封面",
                        data=cover_path.read_bytes(),
                        file_name=cover_path.name,
                        mime="text/markdown",
                        use_container_width=True,
                    )
                for fp in collect_files(artifacts_dir):
                    rel = fp.relative_to(artifacts_dir)
                    st.download_button(
                        f"下载 {rel}",
                        data=fp.read_bytes(),
                        file_name=fp.name,
                    )

st.markdown("---")
st.caption(f"Copyright © 2026 {ORG_NAME} | 项目：{PROJECT_NAME}")
