# run bs cli to create a new project and save the project id
$HOME/bs.exe projects create -n $1 -f csv | grep $1 | cut -d, -f2