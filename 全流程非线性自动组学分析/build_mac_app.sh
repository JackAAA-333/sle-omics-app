#!/bin/zsh
set -euo pipefail
# 生成可双击的 .app（可提交到 Git：启动逻辑用 path to me，无本机硬编码路径）：
# 1) osacompile app_launcher.applescript（要求 .app 与同名的「全流程非线性自动组学分析」目录同父目录）
# 2) 复制图标、删 Assets.car、改 Info.plist 后必须 xattr -cr + codesign
# 用法：./build_mac_app.sh
#      COPY_TO_DESKTOP=1 ./build_mac_app.sh

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
OUT_APP="$ROOT/全流程非线性自动组学分析.app"
ICON="$HERE/logo.icns"

sign_app() {
	local app="$1"
	xattr -cr "$app" 2>/dev/null || true
	codesign --force --deep --sign - "$app"
	codesign --verify "$app"
}

TMP_PARENT="$(mktemp -d /tmp/sle-omics-app.XXXXXX)"
TMP_APP="$TMP_PARENT/全流程非线性自动组学分析.app"

rm -rf "$OUT_APP"
osacompile -o "$TMP_APP" "$HERE/app_launcher.applescript"
cp "$ICON" "$TMP_APP/Contents/Resources/logo.icns"
# osacompile 会带上 Assets.car + CFBundleIconName=applet，Finder 会优先用资产目录而忽略 logo.icns
rm -f "$TMP_APP/Contents/Resources/Assets.car"
/usr/libexec/PlistBuddy -c "Delete :CFBundleIconName" "$TMP_APP/Contents/Info.plist" 2>/dev/null || true
# Apple 约定：CFBundleIconFile 为文件名且不含 .icns 后缀
if plutil -replace CFBundleIconFile -string "logo" "$TMP_APP/Contents/Info.plist" 2>/dev/null; then
	:
else
	/usr/libexec/PlistBuddy -c "Set :CFBundleIconFile logo" "$TMP_APP/Contents/Info.plist" 2>/dev/null || \
		/usr/libexec/PlistBuddy -c "Add :CFBundleIconFile string logo" "$TMP_APP/Contents/Info.plist"
fi
/usr/libexec/PlistBuddy -c "Set :CFBundleDisplayName 全流程非线性自动组学分析" "$TMP_APP/Contents/Info.plist" 2>/dev/null || true

sign_app "$TMP_APP"
ditto --norsrc "$TMP_APP" "$OUT_APP"
sign_app "$OUT_APP"
rm -rf "$TMP_PARENT"

echo "已生成: $OUT_APP"

if [[ "${COPY_TO_DESKTOP:-}" == "1" ]]; then
	rm -rf "$HOME/Desktop/全流程非线性自动组学分析.app"
	ditto --norsrc "$OUT_APP" "$HOME/Desktop/全流程非线性自动组学分析.app"
	sign_app "$HOME/Desktop/全流程非线性自动组学分析.app"
	echo "已复制并签名: $HOME/Desktop/全流程非线性自动组学分析.app"
fi
