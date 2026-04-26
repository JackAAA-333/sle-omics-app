from pathlib import Path
import runpy


APP_PATH = Path(__file__).resolve().parent / "全流程非线性自动组学分析" / "app.py"

if not APP_PATH.exists():
    raise FileNotFoundError(f"未找到应用入口文件: {APP_PATH}")

runpy.run_path(str(APP_PATH), run_name="__main__")
