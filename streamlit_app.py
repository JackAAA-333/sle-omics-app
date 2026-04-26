from pathlib import Path
import runpy
import sys


APP_PATH = Path(__file__).resolve().parent / "全流程非线性自动组学分析" / "app.py"
APP_DIR = APP_PATH.parent

if not APP_PATH.exists():
    raise FileNotFoundError(f"未找到应用入口文件: {APP_PATH}")

if str(APP_DIR) not in sys.path:
    sys.path.insert(0, str(APP_DIR))

runpy.run_path(str(APP_PATH), run_name="__main__")
