# launch app; CSV format doesn't return what we want nicely to using json format; the first Id should be the AppSession number
$HOME/bs.exe launch application -n "$1" --app-version "$2" \
     -o project-id:$3 \
     -o input-type:project \
     -o input-project-id:$3 \
     -o virus-primers:$4 \
     -o basespace-labs-disclaimer:Accepted \
     -f json | grep Id | cut -d'"' -f4 | head -n1