# get the latest DRAGEN app version
$HOME/bs.exe application list -f csv | grep "$1" | cut -d, -f3