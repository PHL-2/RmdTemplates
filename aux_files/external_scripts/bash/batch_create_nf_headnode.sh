#!/bin/bash

# Function to display help message
show_help() {
  echo -e "Usage: $0 [-h] [-r] [-d] [-f <Job definition>] [-i <Docker image>] [-a </path/to/aws/cli>] [-n <Job name>] [-q <Job queue>] [-c <Command>] [-z <Time zone>]"
  echo -e "Options:"
  echo -e "  -h                     Show help message"
  echo -e "  -r                     Set option to register the job definition. Only need to be performed once unless updating the job definition"
  echo -e "                           Use with the -f, -i, and -a options"
  echo -e "  -d                     Set option to describe the job definition"
  echo -e "                           Use with the -f option"
  echo -e "  -f <Job definition>    Name of the job definition (default: nf-headnode)"
  echo -e "  -i <Docker image>      URI of the Docker image hosted on AWS ECR (used only for registering the job definition)"
  echo -e "  -a </path/to/aws/cli>  Path to the AWS CLI binary installed on the AMI (used only for registering the job definition)"
  echo -e "                           Default will be /home/ec2-user/miniconda"
  echo -e "  -n <Job name>          Job name that will appear in the AWS Batch job queue"
  echo -e "  -q <Job queue>         Job queue of the submitted job"
  echo -e "  -c <Command>           Command to run on the AWS Batch head node"
  echo -e "  -z <Time zone>         Time zone to display log timestamps (default: America/New_York)"
  exit 1
}

# Initialize variables
register_job_def="false"
describe_job_def="false"
job_def="nf-headnode"
docker_uri=""
awsclipath="/home/ec2-user/miniconda"
job_name=""
job_queue=""
command=""
time_zone="America/New_York"

# Define options
while getopts "hrdf:i:a:n:q:c:" opt; do
  case $opt in
    h)
      show_help
      ;;
    r)
      register_job_def="true"
      echo "Registering the job definition with the arguments passed to the -f, -i, and -a options"
      ;;
    d)
      describe_job_def="true"
      echo "Describing the job definition with the argument passed to the -f option"
      ;;
    f)
      job_def="$OPTARG"
      echo "Job definition: ${job_def}"
      ;;
    i)
      docker_uri="$OPTARG"
      echo "Docker image: ${docker_uri}"
      ;;
    a)
      awsclipath="$OPTARG"
      # Remove / at the end
      if [[ "${awsclipath: -1}" == "/" ]]; then
        awsclipath="${awsclipath:0:-1}"
      fi
      ;;
    n)
      job_name="$OPTARG"
      echo "Job name: ${job_name}"
      ;;
    q)
      job_queue="$OPTARG"
      echo "Job queue: ${job_queue}"
      ;;
    c)
      command="$OPTARG"
      echo "Command to run: ${command}"
      ;;
    z)
      time_zone="$OPTARG"
      echo "Time zone: ${time_zone}"
      ;;
    \?)
      show_help
      ;;
    :)
      show_help
      ;;
  esac
done

# Remove parsed options from positional parameters
shift $((OPTIND -1))

if [[ -n $@ ]]; then
  echo "Remaining arguments will be ignored: $@"
fi

# Main process
if [[ "$register_job_def" == "true" ]]; then

  if [[ -z "$docker_uri" || -z "$awsclipath" ]]; then
    echo -e "\nError:\nMust pass a Docker image for the job registration\n"
    exit 1
  fi

  # Register the job
  aws batch register-job-definition \
    --job-definition-name "$job_def" \
    --type container \
    --container-properties '{
      "image": "'"$docker_uri"'",
      "vcpus": 2,
      "memory": 4096,
      "command": ["/bin/bash", "-c", "for i in $(seq 1 10); do echo Current time $(date +%T); sleep 60; done"],
      "volumes": [
        {
          "name": "aws-cli",
          "host": {
            "sourcePath": "'"$awsclipath"'"
          }
        }
      ],
      "mountPoints": [
        {
          "sourceVolume": "aws-cli",
          "containerPath": "'"$awsclipath"'",
          "readOnly": true
        }
      ]
    }'

  echo -e "\nJob definition ${job_def} registered\n"
fi

if [[ "$describe_job_def" == "true" ]]; then

  # Check status of job definition
  aws batch describe-job-definitions --job-definition-name "$job_def"
fi

if [[ "$register_job_def" == "false" && "$describe_job_def" == "false" ]]; then

  if [[ -z "$job_name" || -z "$job_queue" || -z "$command" ]]; then
    echo -e "\nError:\n-n, -q, -c options cannot be empty when submitting a new run\n";
    show_help;
    exit 1
  fi

  JOB_ID=$(aws batch submit-job \
    --job-name "$job_name" \
    --job-queue "$job_queue" \
    --job-definition "$job_def" \
    --container-overrides '{"command": ["unbuffer", "/bin/bash", "-c", "'"$command"'"]}' \
    --query jobId --output text)

  echo -e "\nJob ID: $JOB_ID"
  echo -e "If you need to terminate the headnode job, run the following command:\n"
  echo -e "\033[31maws batch terminate-job --job-id $JOB_ID --reason \"Reason for canceling\"\033[0m\n"
  sleep 5

  while true; do
    LOG_STREAM=$(aws batch describe-jobs --jobs "$JOB_ID" \
      --query 'jobs[0].container.logStreamName' --output text)
    if [[ "$LOG_STREAM" != "None" ]]; then
      echo -e "\nLog stream: $LOG_STREAM"
      sleep 20
      break
    fi
    echo -ne "\rWaiting for headnode to start up..."

    # Get the job status
    STATUS=$(aws batch describe-jobs --jobs "$JOB_ID" \
      --query 'jobs[0].status' --output text)

    if [[ "$STATUS" == "FAILED" ]]; then
      echo -e "\nJob finished with status: $STATUS"
      exit 1
    fi

    sleep 20
  done

  # Run aws logs tail until the job completes
  while true; do

    aws logs tail /aws/batch/job \
      --log-stream-names "$LOG_STREAM" \
      --since 60s \
      --format short 2>/dev/null \
      | while read -r line; do
          ts=$(echo "$line" | awk '{print $1}')
          msg=$(echo "$line" | cut -d' ' -f2-)
          echo "$(TZ="$time_zone" date -d "$ts UTC" '+%Y-%m-%d %r %Z') $msg"
        done
    #  | sed -E 's/^[^ ]* //'
    sleep 2

    # Get the job status
    STATUS=$(aws batch describe-jobs --jobs "$JOB_ID" \
      --query 'jobs[0].status' --output text)

    if [[ "$STATUS" == "SUCCEEDED" ]]; then
      echo -e "\nJob finished with status: $STATUS"
      exit 0
    elif [[ "$STATUS" == "FAILED" ]]; then
      echo -e "\nJob finished with status: $STATUS"
      exit 1
    fi
  
    sleep 55
  done
fi