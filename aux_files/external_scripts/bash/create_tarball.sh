#!/bin/bash

# Function to display help message
show_help() {
  echo -e "Usage: $0 [-h] [-t </tmp/dir/>] [-i <MiSeq/NextSeq2000>] [-d <YYYY-MM-DD>] [-b] [-n <hostname>] [-p </path/to/runs/>] [-s <s3://path/to/bucket/>] [-r <prefix>] [-o] [-c]"
  echo -e "Options:"
  echo -e "  -h                         Show help message"
  echo -e "  -t </tmp/dir/>             Temporary directory to hold sequencing run. Default will be ~/tmp_bs_dl/"
  echo -e "  -i <MiSeq/NextSeq2000>     Sequencing instrument"
  echo -e "  -d <YYYY-MM-DD>            Date of the sequencing run"
  echo -e "  -b                         Set option if sequencing run to archive was already uploaded to BaseSpace"
  echo -e "  -n <hostname>              Remote hostname of sequencing instrument in .ssh/config file"
  echo -e "  -p </path/to/runs/>        Path to the sequencing directory in remote hostname"
  echo -e "                               MiSeq default: '' (blank as it is assumed an sFTP server is set up on 'D:\Illumina\MiSeqOutput')"
  echo -e "                               NextSeq2000 default: '/usr/local/illumina/runs'"
  echo -e "  -s <s3://path/to/bucket/>  AWS S3 bucket path to upload sequencing run. Note that the full path must already exist"
  echo -e "  -r <prefix>                Provide a prefix to upload a record (no data files) of the sequencing run to BaseSpace"
  echo -e "                               This string will be added to the front of the Run Name"
  echo -e "  -o                         Set option to upload a dummy SampleSheet.csv to BaseSpace. To be used with the -r option"
  echo -e "                               Otherwise, a format version 1 sample sheet may cause the run record to have a 'Failed' status on BaseSpace"
  echo -e "  -c                         Set option to clean up all data files after running this script. Do not set if want to use local files"
  exit 1
}

# Function to create tar with progress
tar_w_progress() {
  local input_fp=$1
  local output_fp=$2
  local checkpoint_char="="
  local bars=100

  cd "$input_fp" || exit 1

  files_to_tar=$(find . \( -wholename "./*SampleSheet*.csv" -o -name "*.xml" -o -path "./Data" -o -path "./InterOp" \) \
                 ! -path "./Alignment_*" ! -path "./Analysis/*")

  checkpoint_N=$(du -akc --apparent-size $files_to_tar | tail -n 1 | cut -f1 | xargs printf "%d/$((bars-1))" | xargs echo | bc)

  echo -n "Estimate: "
  printf "[%${bars}s]\n" | tr ' ' "${checkpoint_char}"
  echo -n "Progress: [ "

  tar \
    --checkpoint="$checkpoint_N" \
    --checkpoint-action="ttyout=\b${checkpoint_char}>" \
    --record-size=1K \
    --create $files_to_tar | \
    gzip -n > ${output_fp}.tar.gz

  echo ']'
  echo 'Tar completed!'
}

# Initialize variables
tmp_dir=~/tmp_bs_dl/
instrument=""
date_string=""
basespace="false"
hostname=""
hostname_path=""
s3_path=""
record_prefix=""
overwrite="false"
clean_up_local="false"


# Define options
while getopts "ht:i:d:bn:p:s:r:oc" opt; do
  case $opt in
    h)
      show_help
      ;;
    t)
      tmp_dir="$OPTARG"
      # Append / to the end if not provided
      if [[ "${tmp_dir: -1}" != "/" ]]; then
        tmp_dir="${tmp_dir}/"
      fi
      ;;
    i)
      instrument="$OPTARG"
      echo "Instrument: ${instrument}"
      if [[ "$instrument" != "MiSeq" && "$instrument" != "NextSeq2000" ]]; then
        echo -e "\nError:\nInstrument option either needs to be MiSeq or NextSeq2000\n"
        show_help
      fi
      ;;
    d)
      date_string="$OPTARG"
      echo "Date: ${date_string}"
      if [[ ! "$date_string" =~ ^[2-9]{1}[0-9]{3}-[01]{1}[0-9]{1}-[0123]{1}[0-9]{1}$ ]]; then
        echo -e "\nError:\nInvalid date format. Use YYYY-MM-DD\n"
        show_help
      fi
      ;;
    b)
      basespace="true"
      echo "b flag set: Sequencing run will be downloaded from BaseSpace"
      ;;
    n)
      hostname="$OPTARG"
      echo "HostName: ${hostname}"
      ;;
    p)
      hostname_path="$OPTARG"
      echo "Path to sequencing directory on remote: ${hostname_path}"
      ;;
    s)
      s3_path="$OPTARG"
      # Append s3:// to the front if not provided
      if [[ "${s3_path:0:5}" != "s3://" ]]; then
        s3_path="s3://${s3_path}"
      fi
      # Append / to the end if not provided
      if [[ "${s3_path: -1}" != "/" ]]; then
        s3_path="${s3_path}/"
      fi
      echo "Existing S3 bucket path: ${s3_path}"
      ;;
    r)
      record_prefix="$OPTARG"
      echo "A record of the sequencing run (prepended with '${record_prefix}') will be created on BaseSpace"
      ;;
    o)
      overwrite="true"
      ;;
    c)
      clean_up_local="true"
      echo "c flag set: The files generated in the local temporary will be deleted"
      echo "            Stop the script now and remove this flag if you only want the local files"
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

# Check arguments
if [[ -z "$instrument" ]]; then
  echo -e "\nError:\nMust pass either MiSeq or NextSeq2000 to the -i option\n"
  exit 1
elif [[ -z "$date_string" ]]; then
  echo -e "\nError:\nMust pass a date in the format YYYY-MM-DD to the -d option\n"
  exit 1
elif [[ "$basespace" == "false" && -z "$hostname" ]]; then
  echo -e "\nError:\n-b flag is not set; therefore a HostName to the sequencer must be passed -n option\n";
  exit 1
fi

# Check bs and aws cli
if [[ "$basespace" == "true" || -n "$record_prefix" ]]; then
  echo -e "\nChecking if bs cli is available and working..."
  invisible_whoami=$(bs whoami) || bs_error_code=$?

  if [[ "$bs_error_code" -ne 0 ]];then
    echo -e "\nError:\n'bs whoami' command failed with exit status ${bs_error_code}\n"
    exit 1
  fi
fi

if [[ -n "$s3_path" ]]; then
  echo -e "\nChecking if aws cli is available and working..."
  invisible_bucket=$(aws s3 ls "${s3_path}") || s3_error_code=$?

  if [[ "$s3_error_code" -ne 0 ]];then
    echo -e "\nError:\n'aws s3 ls ${s3_path}' command failed with exit status ${s3_error_code}"
    echo -e "If the provided S3 bucket path does not exist yet, upload an empty file to create it first\n"
    echo -e "Example:\n$ touch file.txt\n$ aws s3 cp file.txt ${s3_path}\n"
    exit 1
  fi
fi

if [[ "$instrument" == "MiSeq" ]]; then
  instrument_init="M"
elif [[ "$instrument" == "NextSeq2000" ]]; then
  instrument_init="VH"
else
  exit 1
fi

# Define regex pattern for Illumina sequencing folder
date_reformat="${date_string:2:2}${date_string:5:2}${date_string:8:2}"
run_regex="${date_reformat}_${instrument_init}[0-9]*_[0-9]*_[0-9A-Z-]*"
lower_case_instru=$(echo "$instrument" | tr '[:upper:]' '[:lower:]')
extend_tmp_dir="${tmp_dir}${date_string}-${lower_case_instru}/"

# Main process
# Download sequencing run from BaseSpace
if [[ "$basespace" == "true" ]];then
  bs_entry=$(bs list runs -f csv | grep -E "$run_regex" | grep -vE "$record_prefix")
  num_entries=$(bs list runs -f csv | grep -E "$run_regex" | grep -vE "$record_prefix" | wc -l)

  if [[ -z "$bs_entry" ]]; then
    echo -e "\nError:\nSequencing run loaded on ${date_string} from the ${instrument} is not on BaseSpace\n"
    exit 1
  elif [[ "$num_entries" -gt 1 ]]; then
    echo -e "\nError:\nThis script cannot handle multiple runs with the same date and sequencer\n"
    echo "$bs_entry"
    exit 1
  fi

  run_id=$(echo "$bs_entry" | cut -f2 -d,)
  sequencing_run=$(echo "$bs_entry" | cut -f1 -d,)
  actual_temp_dir="${extend_tmp_dir}${sequencing_run}/"

  echo -e "\nDownloading BaseSpace run ${sequencing_run} to: ${extend_tmp_dir}"
  bs download run --id "$run_id" --output "$actual_temp_dir"

  echo -e "\nCreating tarball of ${sequencing_run} in ${extend_tmp_dir}:\n"
  tar_w_progress "${extend_tmp_dir}${sequencing_run}" "${extend_tmp_dir}${sequencing_run}"
else
# Download sequencing run from sequencers
  # If hostname path for the NextSeq was not provided
  if [[ -z "$hostname_path" && "$instrument" == "NextSeq2000" ]]; then
    hostname_path="/usr/local/illumina/runs"
  fi

  completed_run=$(echo "ls ${hostname_path}/*/RunCompletionStatus.xml" | sftp "$hostname" | grep "$run_regex") || sftp_error_code=$?
  sequencing_run=$(echo "$completed_run" | sed "s/ /\n/g" | rev | cut -f2 -d"/" | rev)
  num_entries=$(echo "$sequencing_run" | wc -l)

  if [[ "$sftp_error_code" -ne 0 ]]; then
    echo -e "\nError:\n'sftp' command to look for the ${date_string} run on ${hostname}:${hostname_path} failed with exit status ${sftp_error_code}\n"
    exit 1
  elif [[ -z "$sequencing_run" ]]; then
    echo -e "\nError:\nRun loaded on ${date_string} could not be found on the ${instrument}\n"
    echo -e "Was there a run loaded on this date?\n"
    echo -e "Otherwise, the sequencer may still be running...\n"
    exit 1
  elif [[ "$num_entries" -gt 1 ]]; then
    echo -e "\nError:\nThis script cannot handle multiple runs on the same sequencer\n"
    echo "$sequencing_run"
    exit 1
  fi

  full_hostname_path="${hostname}:${hostname_path}/${sequencing_run}"

  if [[ "$instrument" == "MiSeq" ]]; then

    actual_temp_dir="${extend_tmp_dir}${sequencing_run}/"
    echo -e "\nDownloading ${instrument} run ${full_hostname_path} to: ${extend_tmp_dir}"
    mkdir -p "$actual_temp_dir"

    scp -rB \
      "${full_hostname_path}/*SampleSheet*.csv" \
      "${full_hostname_path}/*.xml" \
      "${full_hostname_path}/InterOp/" \
      "${full_hostname_path}/Logs/" \
      "${full_hostname_path}/Data/" \
      "$actual_temp_dir"

    echo -e "\nCreating tarball of ${sequencing_run} in ${extend_tmp_dir}:\n"
    tar_w_progress "${extend_tmp_dir}${sequencing_run}" "${extend_tmp_dir}${sequencing_run}"

  elif [[ "$instrument" == "NextSeq2000" ]]; then

  echo -e "\nCreating the tarball in directory - ${full_hostname_path}/"
  ssh -t "$hostname" "$(typeset -f tar_w_progress); tar_w_progress '${hostname_path}/${sequencing_run}' '${hostname_path}/${sequencing_run}/${sequencing_run}'"
  sleep 2

  echo -e "\nDownloading ${instrument} run ${sequencing_run}.tar.gz to: ${extend_tmp_dir}"
  mkdir -p "$extend_tmp_dir"
  scp -rB "$full_hostname_path/${sequencing_run}.tar.gz" "$extend_tmp_dir"
  sleep 2
  fi
fi

# Generate md5sum
echo -e "\nGenerating md5sum file..."
cd "$extend_tmp_dir" && \
  md5sum "${sequencing_run}.tar.gz" > "${extend_tmp_dir}${sequencing_run}.md5"
echo -e "\nArchiving complete! Local files can be found at ${extend_tmp_dir}\n"

# Upload to S3 if -s option provided
if [[ -n "$s3_path" ]]; then
  echo -e "-s option provided. Proceeding to upload to the S3 bucket path"
  echo -e "\nChecking if the run has already been uploaded..."
  existing_run=$(aws s3 ls "${s3_path}${date_string}/${sequencing_run}.tar.gz") || existing_run_code=$?

  if [[ "$existing_run_code" -eq 0 ]];then
    echo -e "\nError:\nThe run already exists at ${s3_path}${date_string}/${sequencing_run}.tar.gz"
    echo -e "To reupload, manually delete this run first\n"
    exit 1
  fi

  aws s3 cp "${extend_tmp_dir}${sequencing_run}.tar.gz" "${s3_path}${date_string}/"
  aws s3 cp "${extend_tmp_dir}${sequencing_run}.md5" "${s3_path}${date_string}/"
fi

# Record the run to BaseSpace
if [[ -n "$record_prefix" ]]; then
  echo -e "\n-r option provided. Proceeding to record the sequencing run to BaseSpace"
  echo -e "\nCleaning any old record files..."
  bs_ul_fp="${extend_tmp_dir}bs_record/${sequencing_run}/"
  rm -rf "$bs_ul_fp"
  mkdir -p "$bs_ul_fp"

  if [[ "$instrument" == "MiSeq" || "$basespace" == "true" ]]; then
    cp -r \
      "${actual_temp_dir}"*"SampleSheet"*".csv" \
      "${actual_temp_dir}"*".xml" \
      "${actual_temp_dir}InterOp/" \
      "$bs_ul_fp"

  # If run is from NextSeq, extract the relevant files from tarball
  elif [[ "$instrument" == "NextSeq2000" ]]; then
    echo -e "\nExtracting the necessary files from ${extend_tmp_dir}${sequencing_run}.tar.gz..."

    tar -xzf "${extend_tmp_dir}${sequencing_run}.tar.gz" \
      --directory "$bs_ul_fp" \
      "./*SampleSheet*.csv" \
      "./InterOp/" \
      "./*.xml" \
      --exclude "./InterOp/C[0-9]*" \
      --exclude "./Config" \
      --exclude "./Logs" \
      --exclude "./Recipe" \
      --checkpoint=1000 \
      --checkpoint-action=ttyout='%{%Y-%m-%d %H:%M:%S}t (%d sec): #%u, %T%*\r'
  fi

  # If the -o option is used, create a new SampleSheet.csv to upload to BaseSpace
  # Uploaded runs will have a 'Failed' or 'Needs Attention' status without the BCLConvert and Cloud sections
  if [[ "$overwrite" == "true" ]]; then

    if [[ -f "${bs_ul_fp}SampleSheet.csv" ]]; then
      mv "${bs_ul_fp}SampleSheet.csv" "${bs_ul_fp}SampleSheet_original.csv"
    fi

    echo -e "A dummy SampleSheet.csv will be created for recording in BaseSpace\n"
    echo "[Header]," > "${bs_ul_fp}SampleSheet.csv"
    echo "FileFormatVersion,2" >> "${bs_ul_fp}SampleSheet.csv"
    echo "" >> "${bs_ul_fp}SampleSheet.csv"
    echo "[BCLConvert_Data]" >> "${bs_ul_fp}SampleSheet.csv"
    echo "Sample_ID" >> "${bs_ul_fp}SampleSheet.csv"
    echo "GenericSampleID" >> "${bs_ul_fp}SampleSheet.csv"
    echo "" >> "${bs_ul_fp}SampleSheet.csv"
    echo "[Cloud_Data]" >> "${bs_ul_fp}SampleSheet.csv"
    echo "Sample_ID" >> "${bs_ul_fp}SampleSheet.csv"
    echo "GenericSampleID" >> "${bs_ul_fp}SampleSheet.csv"
  fi

  bs run upload \
    --name "${record_prefix}${instrument}_${date_string}" \
    --instrument "$instrument" \
    --concurrency "low" \
    --no-progress-bars \
    "$bs_ul_fp"
fi

# Cleanup
if [[ "$clean_up_local" == "true" ]]; then
  echo -e "-c option provided. Deleting the locally downloaded files..."
  sleep 20
  rm -rf "$extend_tmp_dir"
fi
