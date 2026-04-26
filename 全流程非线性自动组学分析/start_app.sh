#!/bin/zsh
set -euo pipefail

PROJECT_ROOT="/Users/jacka/Desktop/SLE分析"
APP_DIR="$PROJECT_ROOT/全流程非线性自动组学分析"
VENV_DIR="$PROJECT_ROOT/.venv"
RUNTIME_DIR="$APP_DIR/.runtime"
PID_FILE="$RUNTIME_DIR/streamlit.pid"
LOG_FILE="$RUNTIME_DIR/streamlit.log"
PORT="8501"
URL="http://127.0.0.1:${PORT}"

mkdir -p "$RUNTIME_DIR"

if [[ -f "$PID_FILE" ]]; then
  EXISTING_PID="$(<"$PID_FILE")"
  if kill -0 "$EXISTING_PID" 2>/dev/null; then
    open "$URL"
    exit 0
  fi
fi

if lsof -iTCP:"$PORT" -sTCP:LISTEN >/dev/null 2>&1; then
  open "$URL"
  exit 0
fi

if [[ ! -x "$VENV_DIR/bin/python" ]]; then
  /usr/bin/python3 -m venv "$VENV_DIR"
fi

source "$VENV_DIR/bin/activate"

if ! python -c "import streamlit" >/dev/null 2>&1; then
  pip install -r "$APP_DIR/requirements.txt" >>"$LOG_FILE" 2>&1
fi

nohup python -m streamlit run "$APP_DIR/app.py" \
  --server.headless true \
  --server.port "$PORT" \
  --server.address 127.0.0.1 >>"$LOG_FILE" 2>&1 &

echo $! >"$PID_FILE"

for _ in {1..20}; do
  if lsof -iTCP:"$PORT" -sTCP:LISTEN >/dev/null 2>&1; then
    break
  fi
  sleep 0.5
done

open "$URL"
