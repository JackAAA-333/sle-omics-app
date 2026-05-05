-- 与 app_launcher.applescript 相同逻辑（静默启动占位）
on run
	set appPath to POSIX path of (path to me)
	if appPath ends with "/" then set appPath to text 1 thru -2 of appPath
	set inner to "APP=" & quoted form of appPath & "; export APP; zsh \"${APP%.app}/start_app.sh\" >/dev/null 2>&1 &"
	do shell script "/bin/zsh -c " & quoted form of inner
end run
