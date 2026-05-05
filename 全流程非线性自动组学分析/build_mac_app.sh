#!/bin/zsh
set -euo pipefail
# 生成可双击的 .app：
# 1) AppleScript 内联绝对路径到 start_app.sh（不依赖 path to me / 额外 launcher）
# 2) 复制图标并更新 Info.plist 后必须 xattr -cr + codesign，否则 Finder 会拒绝（invalid plist / 签名损坏）
# 用法：./build_mac_app.sh
#      COPY_TO_DESKTOP=1 ./build_mac_app.sh

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
OUT_APP="$ROOT/全流程非线性自动组学分析.app"
ICON="$HERE/logo.icns"
START_SH="$HERE/start_app.sh"

sign_app() {
	local app="$1"
	xattr -cr "$app" 2>/dev/null || true
	codesign --force --deep --sign - "$app"
	codesign --verify "$app"
}

TMP_PARENT="$(mktemp -d /tmp/sle-omics-app.XXXXXX)"
TMP_APP="$TMP_PARENT/全流程非线性自动组学分析.app"
ASM="$TMP_PARENT/launcher.applescript"

# AppleScript 源由 shell 展开路径（路径中含 " 则不支持，一般无此情况）
cat >"$ASM" <<APPLESCRIPT
on run
	do shell script "zsh " & quoted form of "${START_SH}" & " >/dev/null 2>&1 &"
end run
APPLESCRIPT

rm -rf "$OUT_APP"
osacompile -o "$TMP_APP" "$ASM"
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
