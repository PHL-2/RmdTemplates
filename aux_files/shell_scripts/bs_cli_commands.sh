#run bs cli to create a new project and save the project id
PROJECT_ID=$(HTTPS_PROXY=proxy.phila.gov:8080 $HOME/bs.exe projects create -n $1 -f csv | grep $1 | cut -d, -f2)

echo "BaseSpace Project was created with ID $PROJECT_ID"

#upload the data to the project ID
HTTPS_PROXY=proxy.phila.gov:8080 $HOME/bs.exe upload dataset -p $PROJECT_ID --recursive $2

#get the latest DRAGEN app version
APP_VERSION=$(HTTPS_PROXY=proxy.phila.gov:8080 $HOME/bs.exe application list -f csv | grep "$3" | cut -d, -f3)

#CSV format doesn't return what we want nicely to using json format; the first Id should be the AppSession number
APP_SESSION_ID=$(HTTPS_PROXY=proxy.phila.gov:8080 $HOME/bs.exe launch application -n "$3" --app-version $APP_VERSION \
  -o project-id:$PROJECT_ID \
  -o input-type:project \
  -o input-project-id:$PROJECT_ID \
  -o virus-primers:$4 \
  -o basespace-labs-disclaimer:Accepted \
  -f json | grep Id | cut -d'"' -f4 | head -n1)

HTTPS_PROXY=proxy.phila.gov:8080 $HOME/bs.exe await $APP_SESSION_ID

HTTPS_PROXY=proxy.phila.gov:8080 $HOME/bs.exe download project -i $PROJECT_ID \
  --extension=detect.report.csv \
  --extension=lineage_report.csv \
  --extension=nextclade.tsv \
  --extension=json \
  --exclude "kmer_config*" \
  -o "$5" \
  --no-metadata \
  --overwrite

HTTPS_PROXY=proxy.phila.gov:8080 $HOME/bs.exe download project -i $PROJECT_ID \
  --exclude "*" \
  --include "ScrubbedFASTQsNode*.zip" \
  -o "$5" \
  --no-metadata \
  --overwrite

APP_DATE=$(HTTPS_PROXY=proxy.phila.gov:8080 $HOME/bs.exe appsession get --id $APP_SESSION_ID -f json | grep DateCompleted | cut -d'"' -f4)

echo $APP_VERSION
echo $APP_DATE
