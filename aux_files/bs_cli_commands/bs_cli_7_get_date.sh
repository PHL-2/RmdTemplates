# get date of app session
$HOME/bs.exe appsession get --id $1 -f json | grep DateCompleted | cut -d'"' -f4