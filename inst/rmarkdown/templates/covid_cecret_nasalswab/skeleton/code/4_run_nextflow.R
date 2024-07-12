library(here)
library(dplyr)
library(tidyr)
library(stringr)

#This Rscript submits the relevant jobs to Nextflow once the sequencing run has been uploaded

###############
# Manual inputs
###############

ec2_tmp_fp <- "~/tmp_bs_dl/"

update_pangolin_dataset <- TRUE

update_freyja_and_cecret_pipeline <- TRUE

cecret_version <- "master"

#########################
# AWS and sequencing date
#########################

#sequencing date of the run folder should match the RStudio project date
sequencing_date <- gsub("_.*", "", basename(here())) #YYYY-MM-DD

if(sequencing_date == "") {
  stop (simpleError(paste0("Please fill in the correct sequencing date or short project description in ", here("code"), "/4_run_nextflow.R")))
} else if (is.na(as.Date(sequencing_date, "%Y-%m-%d")) | nchar(sequencing_date) == 8) {
  stop (simpleError("Please enter the date into [sequencing_date] as YYYY-MM-DD"))
}

################
# Load functions
################

#this file needs to sit in a [aux_files/r_scripts/functions] directory path above this project directory
tryCatch(
  {
    source(file.path(dirname(here()), "aux_files", "r_scripts", "functions", "R_all_functions_v3.R"))
  },
  error = function(e) {
    stop (simpleError("The R_all_functions_v3.R file needs to sit in a [aux_files/r_scripts/functions] directory path above this project directory"))
  }
)

#############
# Load config
#############

#this file needs to sit in a [aux_files/r_scripts/config] directory path above this project directory
tryCatch(
  {
    source(file.path(dirname(here()), "aux_files", "r_scripts", "config", "config_variables.R"))
  },
  error = function(e) {
    stop (simpleError("The config_variables.R file needs to sit in a [aux_files/r_scripts/config] directory path above this project directory"))
  }
)

#############
# Submit jobs
#############

# If this sample sheet is missing, get it from AWS S3 bucket
sample_sheet_fn <- list.files(here("metadata", "munge"), pattern = "SampleSheet_v2.csv")

if(length(sample_sheet_fn) > 1) {
  stop(simpleError("There are more than 2 sample sheets detected!! Please delete the incorrect one"))
}

sequencer_type <- gsub("^[0-9-]*_(MiSeq|NextSeq1k2k)_.*", "\\1", sample_sheet_fn)

sequencer_regex <- case_when(sequencer_type == "MiSeq" ~ "M",
                             sequencer_type == "NextSeq1k2k" ~ "VH")

intended_sequencing_folder_regex <- paste0(gsub("^..|-", "", sequencing_date), "_", sequencer_regex, "[0-9]*_[0-9]*_[0-9A-Z-]*$")

sample_type_acronym <- gsub(paste0("^[0-9-]*_", sequencer_type, "_|_.*"), "", sample_sheet_fn)

prj_description <- gsub(paste0("^[0-9-]*_.*", sample_type_acronym, "_|_.*"), "", sample_sheet_fn)

nf_demux_samplesheet_path <- paste(s3_run_bucket, sequencing_date,
                                   tolower(paste(sequencing_date, sequencer_type, sample_type_acronym, prj_description, "nf_demux_samplesheet.csv", sep = "_")), sep = "/")

bclconvert_output_path <- paste(s3_fastq_bucket, sequencing_date, sample_type_acronym, prj_description, "processed_bclconvert", sep = "/")

# Demultiplexing
submit_screen_job(message2display = "Demultiplexing with BCLConvert",
                  ec2_login = ec2_hostname,
                  screen_session_name = "demux",
                  command2run = paste("cd ~/.tmp_screen/;",
                                      "nextflow run nf-core/demultiplex",
                                      "-c ~/.nextflow/config",
                                      "-profile", demux_profile,
                                      "-bucket-dir", paste0(s3_nextflow_work_bucket, "/demux_", sample_type_acronym, "_", sequencing_date),
                                      "-resume",
                                      "--input", nf_demux_samplesheet_path,
                                      "--outdir", bclconvert_output_path)
)

check_screen_job(message2display = "Checking BCLConvert job",
                 ec2_login = ec2_hostname,
                 screen_session_name = "demux")

# Checking the demultiplexing results
system2("aws", c("sso login"))
aws_s3_fastq_files <- system2("aws", c("s3 ls", bclconvert_output_path,
                                       "--recursive",
                                       "| grep 'R1_001.fastq.gz$'"), stdout = TRUE)

# If the aws-cli provides an SSL error on local machine, run the command through the instance
if(length(aws_s3_fastq_files) == 0) {
  aws_s3_fastq_files <- system2("ssh", c("-tt", ec2_hostname,
                                         shQuote(paste("aws s3 ls", bclconvert_output_path, "--recursive",
                                                       "| grep 'R1_001.fastq.gz$'"), type = "sh")),
                                stdout = TRUE, stderr = TRUE) %>%
    head(-1)
}

fastq_file_sizes <- aws_s3_fastq_files %>%
  str_split("\\s+") %>%
  do.call("rbind", .) %>%
  as.data.frame() %>%
  `colnames<-`(c("date", "time", "bytes", "filename")) %>%
  select(bytes, filename) %>%
  filter(!grepl("/Alignment_", filename)) %>%
  filter(!grepl("/Fastq/", filename)) %>%
  mutate(sequencing_folder = gsub(".*processed_bclconvert/", "", filename),
         sequencing_folder = gsub("/.*", "", sequencing_folder),
         filename = gsub(".*/", "", filename),
         bytes = as.numeric(bytes)) %>%
  filter(grepl(intended_sequencing_folder_regex, sequencing_folder))

if(nrow(fastq_file_sizes) == 0) {
  stop(simpleError(paste0("\nThere were no FastQ files found at path ", bclconvert_output_path,
                          "\nCheck to see if there was an issue with the demultiplexing of the run\n")))
}

instrument_run_id <- unique(fastq_file_sizes$sequencing_folder)

if(length(instrument_run_id) > 1) {
  stop(simpleWarning(paste0("\nThere are two sequencing runs that matched this date. Make sure you selected the correct sequencer!!!\n",
                            "Currently, you are pulling the sequencing run from the ", sequencer_type)))
}

undetermined_bytes <- fastq_file_sizes %>%
  filter(grepl("Undetermined", filename)) %>%
  select(bytes) %>%
  pull()

if(undetermined_bytes/sum(fastq_file_sizes$bytes) > 0.5) {
  stop(simpleError("Something might've went wrong with the demultiplexing!\nThe unassigned reads makes up more than 50% of the total reads!"))
}

fastq_path <- paste(bclconvert_output_path, instrument_run_id, sep = "/")

# If samples for the same project are split across two different runs, move the fastq files to a new bucket path before running the workflow
# system2("aws", c("s3 mv",
#                  fastq_path,
#                  paste(bclconvert_output_path, "fastq-files", sep = "/"),
#                  "--recursive",
#                  "--exclude '*'",
#                  "--include '*_R[12]_001.fastq.gz'"))
#
# fastq_path <- paste(bclconvert_output_path, "fastq-files", sep = "/")

# Download most recent Nextclade pangolin dataset
if(update_pangolin_dataset) {
  submit_screen_job(message2display = "Downloading Nextclade SARS-CoV-2 data",
                    ec2_login = ec2_hostname,
                    screen_session_name = "nextclade-dl",
                    command2run = paste("mkdir -p ~/.local/bin/;",
                                        "wget -q https://github.com/nextstrain/nextclade/releases/latest/download/nextclade-x86_64-unknown-linux-gnu -O ~/.local/bin/nextclade;",
                                        "chmod +x ~/.local/bin/nextclade;",
                                        "nextclade --version;",
                                        "nextclade dataset list --name sars-cov-2 --json > ~/nextclade-sars.json;",
                                        "nextclade dataset get --name sars-cov-2 --output-zip ~/sars.zip;",
                                        "aws s3 cp ~/nextclade-sars.json", paste0(s3_reference_bucket, "/nextclade/nextclade-sars.json;"),
                                        "aws s3 cp ~/sars.zip", paste0(s3_reference_bucket, "/nextclade/sars.zip;"),
                                        "rm ~/nextclade-sars.json ~/sars.zip")
  )

  check_screen_job(message2display = "Checking Nextclade download job",
                   ec2_login = ec2_hostname,
                   screen_session_name = "nextclade-dl")
}

# Update the Cecret pipeline; this should be done as often as possible as it also updates the freyja data used for assignment
if (update_freyja_and_cecret_pipeline) {
  submit_screen_job(message2display = "Updating Cecret pipeline",
                    ec2_login = ec2_hostname,
                    screen_session_name = "update-cecret",
                    command2run = "nextflow pull UPHL-BioNGS/Cecret -r master"
  )

  check_screen_job(message2display = "Checking Cecret update",
                   ec2_login = ec2_hostname,
                   screen_session_name = "update-cecret")
}

workflow_output_fp <- paste(s3_nextflow_output_bucket, "cecret", sample_type_acronym, paste0(sequencing_date, "_", prj_description), sep = "/")

# Cecret pipeline
submit_screen_job(message2display = "Processing data through Cecret pipeline",
                  ec2_login = ec2_hostname,
                  screen_session_name = "cecret",
                  command2run = paste("cd ~/.tmp_screen/;",
                                      "nextflow run UPHL-BioNGS/Cecret",
                                      "-profile", cecret_profile,
                                      "-bucket-dir", paste0(s3_nextflow_work_bucket, "/cecret_", sample_type_acronym, "_", sequencing_date),
                                      "-r", cecret_version,
                                      "-resume",
                                      "--reads", fastq_path,
                                      "--outdir", paste(workflow_output_fp, "processed_cecret", sep = "/"))
)

check_screen_job(message2display = "Checking Cecret job",
                 ec2_login = ec2_hostname,
                 screen_session_name = "cecret")

rstudioapi::executeCommand('activateConsole')

# Download BCLConvert files
bcl_file_download_command <- c("s3 cp", paste(s3_fastq_bucket, sequencing_date, sample_type_acronym, prj_description, sep = "/"))
bcl_file_download_param <- c("--recursive",
                             "--exclude '*'",
                             "--include '*software_versions.yml'",
                             "--include '*Demultiplex_Stats.csv'",
                             "--include '*Quality_Metrics.csv'",
                             "--include '*Top_Unknown_Barcodes.csv'",
                             "--include '*fastq.gz_fastqc_data.txt'")

aws_s3_bcl_download <- system2("aws", c(bcl_file_download_command,
                                        here("data"),
                                        bcl_file_download_param),
                               stdout = TRUE, stderr = TRUE)

# Download Cecret files
cecret_file_download_command <- c("s3 cp", workflow_output_fp)
cecret_file_download_param <- c("--recursive",
                                "--exclude '*'",
                                #                                "--include '*.consensus.fa'",
                                #                                "--include '*_filtered_R[12].fastq.gz'",
                                "--include '*_demix.tsv'",
                                "--include '*_kraken2_report.txt'",
                                "--include '*snp-dists.txt'",
                                "--include '*_amplicon_depth.csv'",
                                "--include '*nextclade.tsv'",
                                "--include '*nextclade.json'",
                                "--include '*.stats.txt'",
                                "--include '*.cov.txt'",
                                "--include '*.depth.txt'",
                                "--include '*cecret_results.csv'")

aws_s3_cecret_download <- system2("aws", c(cecret_file_download_command,
                                           here("data"),
                                           cecret_file_download_param),
                                  stdout = TRUE, stderr = TRUE)

# Download Nextclade dataset
nextclade_dataset_download_command <- c("s3 cp", paste0(s3_reference_bucket, "/nextclade/nextclade-sars.json"))

aws_s3_nextclade_dataset_download <- system2("aws", c(nextclade_dataset_download_command,
                                                      here("data", "processed_cecret", "nextclade/")),
                                             stdout = TRUE, stderr = TRUE)

# If the aws-cli provides an SSL error on local machine, run the command through the instance
if(any(grepl("fatal error", c(aws_s3_bcl_download, aws_s3_cecret_download, aws_s3_nextclade_dataset_download)))) {

  # Download BCLConvert data
  submit_screen_job(message2display = "Downloading BCLConvert data",
                    ec2_login = ec2_hostname,
                    screen_session_name = "bclconvert-dl",
                    command2run = paste("mkdir -p", paste0(ec2_tmp_fp, sequencing_date, "/data;"),
                                        "aws",
                                        paste0(bcl_file_download_command, collapse = " "),
                                        paste0(ec2_tmp_fp, sequencing_date, "/data"),
                                        paste0(bcl_file_download_param, collapse = " "))
  )

  check_screen_job(message2display = "Checking BCLConvert download job",
                   ec2_login = ec2_hostname,
                   screen_session_name = "bclconvert-dl")

  # Download Cecret data
  submit_screen_job(message2display = "Downloading Cecret data",
                    ec2_login = ec2_hostname,
                    screen_session_name = "cecret-dl",
                    command2run = paste("mkdir -p", paste0(ec2_tmp_fp, sequencing_date, "/data;"),
                                        "aws",
                                        paste0(cecret_file_download_command, collapse = " "),
                                        paste0(ec2_tmp_fp, sequencing_date, "/data"),
                                        paste0(cecret_file_download_param, collapse = " "))
  )

  check_screen_job(message2display = "Checking Cecret download job",
                   ec2_login = ec2_hostname,
                   screen_session_name = "cecret-dl")

  # Download Nextclade dataset
  submit_screen_job(message2display = "Downloading Nextclade dataset version",
                    ec2_login = ec2_hostname,
                    screen_session_name = "nextclade-version-dl",
                    command2run = paste("mkdir -p", paste0(ec2_tmp_fp, sequencing_date, "/data;"),
                                        "aws",
                                        paste0(nextclade_dataset_download_command, collapse = " "),
                                        paste0(ec2_tmp_fp, sequencing_date, "/data/","processed_cecret/", "nextclade/"))
  )

  check_screen_job(message2display = "Checking Nextclade version download job",
                   ec2_login = ec2_hostname,
                   screen_session_name = "nextclade-version-dl")

  # Download data from EC2 instance to local
  run_in_terminal(paste("scp -r", paste0(ec2_hostname, ":", ec2_tmp_fp, sequencing_date, "/data"),
                        here())
  )
}

# Download Nextflow config file for profile (use terminal because of proxy login issue)
run_in_terminal(paste("scp", paste0(ec2_hostname, ":~/.nextflow/config"),
                      here("data", "processed_cecret", "nextflow.config")),
                paste(" [On", ec2_hostname, "instance]\n",
                      "aws s3 cp ~/.nextflow/config",
                      paste0("s3://test-environment/input/", sequencing_date, "/"), "\n\n",
                      "[On local computer]\n",
                      "aws s3 cp", paste0("s3://test-environment/input/", sequencing_date, "/config"),
                      here("data", "processed_cecret", "nextflow.config"))
)

# Clean up environment
submit_screen_job(message2display = "Cleaning up EC2 run folder",
                  ec2_login = ec2_hostname,
                  screen_session_name = "delete-run",
                  command2run = paste0("rm -rf ", ec2_tmp_fp, ";",
                                       "echo Here are your files and directories at home:;",
                                       "ls ~ -GF")
)

check_screen_job(message2display = "Checking delete job",
                 ec2_login = ec2_hostname,
                 screen_session_name = "delete-run")

rstudioapi::executeCommand('activateConsole')
