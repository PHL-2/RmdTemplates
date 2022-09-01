#each command has a retry loop for 5 tries
#run bs cli to create a new project and save the project id
n=0
until [ "$n" -ge 5 ]
do
   PROJECT_ID=$(HTTPS_PROXY=proxy.phila.gov:8080 $HOME/bs.exe projects create -n $1 -f csv | grep $1 | cut -d, -f2) && break
   n=$((n+1))
   echo "Creating project tries: $n times"
   sleep 15
done

echo "BaseSpace Project was created with ID $PROJECT_ID"

sleep 5

#upload the data to the project ID
n=0
until [ "$n" -ge 5 ]
do
   HTTPS_PROXY=proxy.phila.gov:8080 $HOME/bs.exe upload dataset -p $PROJECT_ID --recursive "$2" && break
   n=$((n+1))
   echo "Uploading samples tries: $n times"
   sleep 15
done


sleep 5

#get the latest DRAGEN app version
n=0
until [ "$n" -ge 5 ]
do
   APP_VERSION=$(HTTPS_PROXY=proxy.phila.gov:8080 $HOME/bs.exe application list -f csv | grep "$3" | cut -d, -f3) && break
   n=$((n+1))
   echo "Getting app version tries: $n times"
   sleep 15
done

sleep 5

#CSV format doesn't return what we want nicely to using json format; the first Id should be the AppSession number
n=0
until [ "$n" -ge 5 ]
do
  APP_SESSION_ID=$(HTTPS_PROXY=proxy.phila.gov:8080 $HOME/bs.exe launch application -n "$3" --app-version $APP_VERSION \
    -o project-id:$PROJECT_ID \
    -o input-type:project \
    -o input-project-id:$PROJECT_ID \
    -o virus-primers:$4 \
    -o basespace-labs-disclaimer:Accepted \
    -f json | grep Id | cut -d'"' -f4 | head -n1) && break
   n=$((n+1))
   echo "Launching app tries: $n times"
   sleep 15
done

#wait until the app finishes running
n=0
until [ "$n" -ge 5 ]
do
   HTTPS_PROXY=proxy.phila.gov:8080 $HOME/bs.exe await appsession $APP_SESSION_ID && break
   n=$((n+1))
   echo "Awaiting app tries: $n times"
   sleep 15
done

sleep 5

n=0
until [ "$n" -ge 5 ]
do
   HTTPS_PROXY=proxy.phila.gov:8080 $HOME/bs.exe download project -i $PROJECT_ID \
    --exclude "*" \
    --include "*.zip" \
    --include "Undetermined_report_metrics.json" \
    --include "*detect.report.csv" \
    --include "*lineage_report.csv" \
    --include "*nextclade.tsv" \
    -o "$5" \
    --no-metadata \
    --overwrite && break
   n=$((n+1))
   echo "Download files tries: $n times"
   sleep 15
done

sleep 5

n=0
until [ "$n" -ge 5 ]
do
   APP_DATE=$(HTTPS_PROXY=proxy.phila.gov:8080 $HOME/bs.exe appsession get --id $APP_SESSION_ID -f json | grep DateCompleted | cut -d'"' -f4) && break
   n=$((n+1))
   echo "Getting date tries: $n times"
   sleep 15
done

echo $APP_VERSION
echo $APP_DATE
