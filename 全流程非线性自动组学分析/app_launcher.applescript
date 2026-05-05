-- 与仓库根目录旁的「全流程非线性自动组学分析」文件夹配套使用：
-- …/全流程非线性自动组学分析.app 与 …/全流程非线性自动组学分析/start_app.sh 同父目录。
-- build_mac_app.sh 会据此源文件 osacompile 并处理图标与签名。
on run
	set appPath to POSIX path of (path to me)
	if appPath ends with "/" then set appPath to text 1 thru -2 of appPath
	set inner to "APP=" & quoted form of appPath & "; export APP; zsh \"${APP%.app}/start_app.sh\" >/dev/null 2>&1 &"
	do shell script "/bin/zsh -c " & quoted form of inner
end run
