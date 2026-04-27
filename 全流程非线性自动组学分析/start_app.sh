#!/bin/zsh
set -euo pipefail

SOURCE_ROOT="/Users/jacka/Desktop/SLE分析"
WORK_ROOT="$HOME/Library/Application Support/SLE-OmiX-App"
APP_RELATIVE_DIR="全流程非线性自动组学分析"
PROJECT_ROOT="$WORK_ROOT"
APP_DIR="$PROJECT_ROOT/$APP_RELATIVE_DIR"
VENV_DIR="$PROJECT_ROOT/.venv"
RUNTIME_DIR="$APP_DIR/.runtime"
PID_FILE="$RUNTIME_DIR/streamlit.pid"
LOG_FILE="$RUNTIME_DIR/streamlit.log"

mkdir -p "$PROJECT_ROOT"

# Run from Application Support to avoid desktop folder permission restrictions for app bundles.
rsync -a --delete \
  --exclude ".git/" \
  --exclude ".venv/" \
  --exclude ".history/" \
  --exclude "outputs/" \
  --exclude "outputs_*/" \
  --exclude "*.app/" \
  "$SOURCE_ROOT/" "$PROJECT_ROOT/"

mkdir -p "$RUNTIME_DIR"

if [[ ! -x "$VENV_DIR/bin/python" ]]; then
  if [[ -x "/opt/homebrew/bin/python3.14" ]]; then
    "/opt/homebrew/bin/python3.14" -m venv "$VENV_DIR"
  else
    python3 -m venv "$VENV_DIR"
  fi
fi

source "$VENV_DIR/bin/activate"
python -m pip install --upgrade pip >>"$LOG_FILE" 2>&1

if ! python -c "import streamlit, webview" >/dev/null 2>&1; then
  pip install -r "$PROJECT_ROOT/requirements.txt" >>"$LOG_FILE" 2>&1
fi

python "$APP_DIR/desktop_launcher.py" >>"$LOG_FILE" 2>&1
