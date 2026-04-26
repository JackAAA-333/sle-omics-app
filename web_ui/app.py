import tempfile
import zipfile
from pathlib import Path

import streamlit as st

try:
    import run_pipeline
except ImportError:
    from web_ui import run_pipeline


st.set_page_config(page_title="SLE 多组学自动分析平台", page_icon="🧬", layout="wide")


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


st.title("🧬 SLE 多组学自动对齐分析平台")
st.caption("上传蛋白组学与代谢组学数据，自动执行全流程分析并生成可下载结果。")

left, right = st.columns(2)
with left:
    protein_file = st.file_uploader("蛋白组数据", type=["csv", "tsv", "xlsx"], key="prot")
with right:
    metabolite_file = st.file_uploader("代谢组数据", type=["csv", "tsv", "xlsx"], key="met")

outdir = st.text_input("结果保存目录", value="outputs_web")
run_button = st.button("开始自动分析", type="primary", use_container_width=True)

st.markdown("### 流程说明")
st.markdown(
    "- 自动数据规范化与样本对齐\n"
    "- 代谢组/蛋白组单组学分析\n"
    "- 多组学整合建模（XGBoost + SHAP）\n"
    "- 自动汇总报告与图表归档"
)

if run_button:
    if not protein_file or not metabolite_file:
        st.error("请先上传蛋白组和代谢组文件。")
    else:
        staging = Path(outdir) / "_uploads"
        saved_prot = save_upload(protein_file, staging, f"protein{Path(protein_file.name).suffix.lower()}")
        saved_met = save_upload(metabolite_file, staging, f"metabolite{Path(metabolite_file.name).suffix.lower()}")

        st.info("任务已启动，正在执行完整分析流程。")
        log_area = st.empty()
        logs = []
        result = {}

        generator = run_pipeline.run_pipeline(saved_prot, saved_met, outdir=outdir)
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
                make_zip(artifacts_dir, zip_path)
                st.download_button(
                    "下载全部结果 (ZIP)",
                    data=zip_path.read_bytes(),
                    file_name=zip_path.name,
                    mime="application/zip",
                    use_container_width=True,
                )
                for fp in collect_files(artifacts_dir):
                    rel = fp.relative_to(artifacts_dir)
                    st.download_button(
                        f"下载 {rel}",
                        data=fp.read_bytes(),
                        file_name=fp.name,
                    )
