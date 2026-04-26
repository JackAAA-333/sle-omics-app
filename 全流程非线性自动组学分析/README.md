# 全流程非线性自动组学分析

该页面提供一个可视化入口：上传蛋白组学和代谢组学文件后，自动执行全流程非线性分析并输出结果包。

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
pip install -r "全流程非线性自动组学分析/requirements.txt"
```

2. 启动页面

```bash
cd "全流程非线性自动组学分析"
streamlit run app.py
```

3. 上传蛋白组学和代谢组学数据，点击“开始自动分析”。

## 双击即用（macOS App）

已提供启动脚本与 AppleScript：

- `start_app.sh`
- `app_launcher.applescript`

可在项目根目录执行：

```bash
osacompile -o "/Users/jacka/Desktop/SLE分析/全流程非线性自动组学分析.app" "/Users/jacka/Desktop/SLE分析/全流程非线性自动组学分析/app_launcher.applescript"
```

之后双击 `全流程非线性自动组学分析.app` 即可自动启动并打开页面。
