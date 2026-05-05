#!/bin/zsh
set -euo pipefail

# 仓库根（含 preprocess_and_analyze.py 等）；可由 .app 内 launcher 注入
SOURCE_ROOT="${SLE_OMIX_SOURCE_ROOT:-$(cd "$(dirname "$0")/.." && pwd)}"
WORK_ROOT="$HOME/Library/Application Support/SLE-OmiX-App"
APP_RELATIVE_DIR="全流程非线性自动组学分析"
PROJECT_ROOT="$WORK_ROOT"
APP_DIR="$PROJECT_ROOT/$APP_RELATIVE_DIR"
VENV_DIR="$PROJECT_ROOT/.venv"
RUNTIME_DIR="$APP_DIR/.runtime"
PID_FILE="$RUNTIME_DIR/streamlit.pid"
LOG_FILE="$RUNTIME_DIR/streamlit.log"

mkdir -p "$PROJECT_ROOT" "$RUNTIME_DIR"

# 快速重进：8501 已有 Streamlit 在跑则跳过 rsync / pip，直接拉起窗口（改代码后需先在页面里停服务或杀进程再开，以触发同步）
if [[ -x "$VENV_DIR/bin/python" ]] && "$VENV_DIR/bin/python" -c "
import socket
s = socket.socket()
s.settimeout(0.25)
try:
    raise SystemExit(0 if s.connect_ex(('127.0.0.1', 8501)) == 0 else 1)
finally:
    s.close()
" 2>/dev/null; then
	source "$VENV_DIR/bin/activate"
	exec python "$APP_DIR/desktop_launcher.py" >>"$LOG_FILE" 2>&1
fi

# 可选：开发机跳过同步（export SLE_OMIX_SKIP_RSYNC=1）；首次安装勿用
if [[ "${SLE_OMIX_SKIP_RSYNC:-}" != "1" ]] || [[ ! -f "$APP_DIR/desktop_launcher.py" ]]; then
	rsync -a --delete \
		--exclude ".git/" \
		--exclude ".venv/" \
		--exclude ".history/" \
		--exclude "outputs/" \
		--exclude "outputs_*/" \
		--exclude "*.app/" \
		--exclude "__pycache__/" \
		--exclude "*.pyc" \
		--exclude ".pytest_cache/" \
		--exclude ".mypy_cache/" \
		--exclude ".cursor/" \
		"$SOURCE_ROOT/" "$PROJECT_ROOT/"
fi

if [[ ! -x "$VENV_DIR/bin/python" ]]; then
	if [[ -x "/opt/homebrew/bin/python3.14" ]]; then
		"/opt/homebrew/bin/python3.14" -m venv "$VENV_DIR"
	else
		python3 -m venv "$VENV_DIR"
	fi
fi

source "$VENV_DIR/bin/activate"

# 仅在新创建的 venv 里升级 pip；日常启动不再跑，加快冷启动
if [[ ! -f "$VENV_DIR/.sle_pip_bootstrapped" ]]; then
	python -m pip install --upgrade pip >>"$LOG_FILE" 2>&1
	touch "$VENV_DIR/.sle_pip_bootstrapped"
fi

if ! python -c "import streamlit, webview, yaml" >/dev/null 2>&1; then
	pip install -r "$PROJECT_ROOT/requirements.txt" >>"$LOG_FILE" 2>&1
fi

python "$APP_DIR/desktop_launcher.py" >>"$LOG_FILE" 2>&1
