import socket
import subprocess
import sys
import time
import traceback
from pathlib import Path


PROJECT_NAME = "全流程非线性自动组学分析"
APP_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = APP_DIR.parent
RUNTIME_DIR = APP_DIR / ".runtime"
LOG_FILE = RUNTIME_DIR / "streamlit.log"
PID_FILE = RUNTIME_DIR / "streamlit.pid"
PORT = 8501
URL = f"http://127.0.0.1:{PORT}"
LAUNCH_LOG = RUNTIME_DIR / "desktop_launcher.log"


def is_port_open(port: int) -> bool:
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
        sock.settimeout(0.5)
        return sock.connect_ex(("127.0.0.1", port)) == 0


def start_streamlit_if_needed() -> None:
    if is_port_open(PORT):
        return

    RUNTIME_DIR.mkdir(parents=True, exist_ok=True)
    with open(LOG_FILE, "a", encoding="utf-8") as logf:
        proc = subprocess.Popen(
            [
                sys.executable,
                "-m",
                "streamlit",
                "run",
                str(APP_DIR / "app.py"),
                "--server.headless",
                "true",
                "--server.port",
                str(PORT),
                "--server.address",
                "127.0.0.1",
            ],
            cwd=str(PROJECT_ROOT),
            stdout=logf,
            stderr=subprocess.STDOUT,
        )
    PID_FILE.write_text(str(proc.pid), encoding="utf-8")

    for _ in range(40):
        if is_port_open(PORT):
            return
        time.sleep(0.5)
    raise RuntimeError("Streamlit 服务启动超时，请检查 .runtime/streamlit.log")


def write_launch_log(message: str) -> None:
    RUNTIME_DIR.mkdir(parents=True, exist_ok=True)
    with open(LAUNCH_LOG, "a", encoding="utf-8") as f:
        f.write(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}\n")


def open_standalone_browser_window() -> None:
    # Fallback: launch a standalone app-like window (not regular browser tab).
    chrome_paths = [
        "/Applications/Google Chrome.app",
        "/Applications/Chromium.app",
        "/Applications/Microsoft Edge.app",
    ]
    for app_path in chrome_paths:
        if Path(app_path).exists():
            subprocess.Popen(["open", "-na", app_path, "--args", f"--app={URL}"])
            write_launch_log(f"Fallback launched with app mode: {app_path}")
            return
    subprocess.Popen(["open", URL])
    write_launch_log("Fallback launched with default browser.")


def open_desktop_window() -> None:
    # Use standalone app-mode window for reliability on macOS launch contexts.
    open_standalone_browser_window()


if __name__ == "__main__":
    try:
        start_streamlit_if_needed()
        open_desktop_window()
    except Exception as exc:
        write_launch_log(f"launcher fatal: {exc}")
        write_launch_log(traceback.format_exc())
        open_standalone_browser_window()
