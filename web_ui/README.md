# SLE 多组学自动分析 Web UI

该页面提供一个可视化入口：上传蛋白组学和代谢组学文件后，自动执行全流程分析并输出结果包。

## 主要能力

- 自动将上传文件规范化为管线可读格式（xlsx）
- 自动对齐并依次运行：
  - `preprocess_and_analyze.py`
  - `proteomics_analysis.py`
  - `multiomics_integration.py`
  - `generate_reports.py`
- 自动归档 `outputs*` 目录，生成下载 ZIP
- 页面展示关键指标和图像预览

## 快速开始

1. 安装依赖

```bash
source .venv/bin/activate
pip install -r web_ui/requirements.txt
```

2. 启动页面

```bash
cd web_ui
streamlit run app.py
```

3. 上传蛋白组学和代谢组学数据，点击“开始自动分析”。
