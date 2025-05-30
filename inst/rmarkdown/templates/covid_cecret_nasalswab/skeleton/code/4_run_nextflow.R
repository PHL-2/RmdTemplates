library(here)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)

#This Rscript submits the relevant jobs to Nextflow once the sequencing run has been uploaded

#########################
# AWS and sequencing date
#########################

#sequencing date of the run folder should match the RStudio project date
sequencing_date <- gsub("_.*", "", basename(here())) #YYYY-MM-DD

# temporary directory to hold the sequencing run download
ec2_tmp_fp <- "~/tmp_bs_dl"

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

instrument_type <- gsub("^[0-9-]*_(MiSeq|NextSeq2000)_.*", "\\1", sample_sheet_fn)

sequencer_regex <- case_when(instrument_type == "MiSeq" ~ "M",
                             instrument_type == "NextSeq2000" ~ "VH")

intended_sequencing_folder_regex <- paste0(gsub("^..|-", "", sequencing_date), "_", sequencer_regex, "[0-9]*_[0-9]*_[0-9A-Z-]*")

nf_demux_samplesheet_path <- paste(s3_run_bucket, sequencing_date,
                                   tolower(paste(sequencing_date, instrument_type, sample_type_acronym, prj_description, "nf_demux_samplesheet.csv", sep = "_")), sep = "/")

bclconvert_output_path <- paste(s3_fastq_bucket, sequencing_date, sample_type_acronym, prj_description, "processed_bclconvert", sep = "/")

workflow_output_fp <- paste(s3_nextflow_output_bucket, "cecret", sample_type_acronym, paste0(sequencing_date, "_", prj_description), instrument_type, sep = "/")

# temporary directory to hold the screen log files
tmp_screen_fp <- paste("~", ".tmp_screen", instrument_type, paste0(sample_type_acronym, "_", pathogen_acronym), basename(here()), sep = "/")

session_suffix <- tolower(paste(instrument_type, sample_type_acronym, pathogen_acronym, basename(here()), sep = "-"))

data_output_fp <- paste0(ec2_tmp_fp, "/", session_suffix, "/data")

# Demultiplexing
demux_session <- paste0("nf-demux-", session_suffix)
demux_session_fp <- paste0(tmp_screen_fp, "/", demux_session, ";")
submit_screen_job(message2display = "Demultiplexing with BCLConvert",
                  ec2_login = ec2_hostname,
                  screen_session_name = demux_session,
                  screen_log_fp = tmp_screen_fp,
                  command2run = paste("mkdir -p", demux_session_fp,
                                      "cd", demux_session_fp,
                                      "nextflow run nf-core/demultiplex",
                                      "-c ~/.nextflow/config",
                                      "-profile", demux_profile,
                                      "-bucket-dir", paste0(s3_nextflow_work_bucket, "/demux_", sample_type_acronym, "_", sequencing_date),
                                      "-resume",
                                      "-r 1.3.2",
                                      "--input", nf_demux_samplesheet_path,
                                      "--outdir", bclconvert_output_path)
)

check_screen_job(message2display = "Checking BCLConvert job",
                 ec2_login = ec2_hostname,
                 screen_session_name = demux_session,
                 screen_log_fp = tmp_screen_fp)

# Checking the demultiplexing results
aws_s3_fastq_files <- system2("ssh", c("-tt", ec2_hostname,
                                       shQuote(paste("aws s3 ls", bclconvert_output_path, "--recursive",
                                                     "| grep 'R1_001.fastq.gz$'"), type = "sh")),
                              stdout = TRUE, stderr = TRUE) %>%
  head(-1)

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
  filter(grepl(intended_sequencing_folder_regex, sequencing_folder)) %>%
  #sometimes the Undetermined file from the sequencer gets copied over; remove these samples
  filter(!grepl("^GenericSampleID", filename)) %>%
  group_by(filename) %>%
  mutate(file_num = n()) %>%
  filter(!(file_num == 2 & grepl("^Undetermined", filename) & bytes == max(bytes))) %>%
  ungroup()

if(nrow(fastq_file_sizes) == 0) {
  stop(simpleError(paste0("\nThere were no FastQ files found at path ", bclconvert_output_path,
                          "\nCheck to see if there was an issue with the demultiplexing of the run\n")))
}

if(any(fastq_file_sizes$bytes <= 23)) {

  sample_id_no_reads <- fastq_file_sizes %>%
    mutate(filename = gsub("_.*", "", filename)) %>%
    filter(!grepl("Undetermined", filename),
           bytes <= 23)

  rm_bclconvert_session <- paste0("rm-bclconvert-", session_suffix)
  submit_screen_job(message2display = "Removing recently demultiplexed FASTQ files",
                    ec2_login = ec2_hostname,
                    screen_session_name = rm_bclconvert_session,
                    screen_log_fp = tmp_screen_fp,
                    command2run = paste("aws s3 rm", bclconvert_output_path, "--recursive")
  )

  check_screen_job(message2display = "Checking remove FASTQ files job",
                   ec2_login = ec2_hostname,
                   screen_session_name = rm_bclconvert_session,
                   screen_log_fp = tmp_screen_fp)

  stop(simpleError(paste("\nThese sample IDs had no reads:\n",
                         paste0(sample_id_no_reads$filename, collapse = ", "),
                         "\n\nPlease check that the correct IDT set for the barcodes were used",
                         "\n\nFor samples that do not have reads even with the correct barcodes, please remove them from the analysis by",
                         "\n adding their SampleID to the vector sample_w_empty_reads in *_QC_Report.Rmd",
                         "\n\nRegenerate the SampleSheet by first and then rerun this script to re-demultiplex")))
}

instrument_run_id <- unique(fastq_file_sizes$sequencing_folder)

if(length(instrument_run_id) > 1) {
  stop(simpleWarning(paste0("\nThere are two sequencing runs that matched this date. Make sure you selected the correct sequencer!!!\n",
                            "Currently, you are pulling the sequencing run from the ", instrument_type)))
}

if(remove_undetermined_file) {

  # Move the Undetermined file to another bucket path
  undetermined_mv_session <- paste0("undetermined-mv-", session_suffix)
  submit_screen_job(message2display = "Moving Undetermined files out of input filepath",
                    ec2_login = ec2_hostname,
                    screen_session_name = undetermined_mv_session,
                    screen_log_fp = tmp_screen_fp,
                    command2run = paste("aws s3 mv",
                                        paste(bclconvert_output_path, instrument_run_id, sep = "/"),
                                        paste(bclconvert_output_path, "Undetermined", instrument_run_id, sep = "/"),
                                        "--recursive",
                                        "--exclude '*'",
                                        "--include '*Undetermined_S0_*_001.fastq.gz'")
  )

  check_screen_job(message2display = "Checking Undetermined move job",
                   ec2_login = ec2_hostname,
                   screen_session_name = undetermined_mv_session,
                   screen_log_fp = tmp_screen_fp)
} else {
  undetermined_bytes <- fastq_file_sizes %>%
    filter(grepl("Undetermined", filename)) %>%
    select(bytes) %>%
    pull()

  if(undetermined_bytes/sum(fastq_file_sizes$bytes) > 0.5) {
    stop(simpleError(paste("Something might've went wrong with the demultiplexing!",
                           "The unassigned reads makes up more than 50% of the total reads!",
                           "To remove the Undetermined file from processing, set remove_undetermined_file to TRUE", sep = "\n")))
  }
}

fastq_path <- paste(bclconvert_output_path, instrument_run_id, sep = "/")

# Download most recent Nextclade dataset for lineage assignment
if (nextclade_dataset_version != "") {
  nextclade_tag <- paste("--tag", nextclade_dataset_version)
} else {
  nextclade_tag <- ""
}

update_nextclade_session <- paste0("update-nextclade-", session_suffix)
submit_screen_job(message2display = "Downloading Nextclade SARS-CoV-2 data",
                  ec2_login = ec2_hostname,
                  screen_session_name = update_nextclade_session,
                  screen_log_fp = tmp_screen_fp,
                  command2run = paste("mkdir -p ~/.local/bin/;",
                                      "wget -q https://github.com/nextstrain/nextclade/releases/latest/download/nextclade-x86_64-unknown-linux-gnu -O ~/.local/bin/nextclade;",
                                      "chmod +x ~/.local/bin/nextclade;",
                                      "nextclade --version;",
                                      "nextclade dataset get --name sars-cov-2", nextclade_tag, "--output-zip ~/sars.zip;",
                                      "aws s3 cp ~/sars.zip", paste0(s3_reference_bucket, "/nextclade/sars.zip;"),
                                      "rm ~/sars.zip")
)

check_screen_job(message2display = "Checking Nextclade download job",
                 ec2_login = ec2_hostname,
                 screen_session_name = update_nextclade_session,
                 screen_log_fp = tmp_screen_fp)

# Update the Cecret pipeline? Not always necessary if using a stable version
if (update_freyja_and_cecret_pipeline) {
  update_cecret_session <- paste0("update-cecret-", session_suffix)
  submit_screen_job(message2display = "Updating Cecret pipeline",
                    ec2_login = ec2_hostname,
                    screen_session_name = update_cecret_session,
                    screen_log_fp = tmp_screen_fp,
                    command2run = "nextflow pull UPHL-BioNGS/Cecret -r master"
  )

  check_screen_job(message2display = "Checking Cecret update",
                   ec2_login = ec2_hostname,
                   screen_session_name = update_cecret_session,
                   screen_log_fp = tmp_screen_fp)
}

# Cecret pipeline
cecret_session <- paste0("nf-cecret-", session_suffix)
cecret_session_fp <- paste0(tmp_screen_fp, "/", cecret_session, ";")
submit_screen_job(message2display = "Processing data through Cecret pipeline",
                  ec2_login = ec2_hostname,
                  screen_session_name = cecret_session,
                  screen_log_fp = tmp_screen_fp,
                  command2run = paste("mkdir -p", cecret_session_fp,
                                      "cd", cecret_session_fp,
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
                 screen_session_name = cecret_session,
                 screen_log_fp = tmp_screen_fp)

# Download BCLConvert options
bcl_file_download_command <- c("s3 cp", dirname(bclconvert_output_path))
bcl_file_download_param <- c("--recursive",
                             "--exclude '*'",
                             "--include '*/software_versions.yml'",
                             paste0("--include '*", intended_sequencing_folder_regex, "/",
                                    c("Demultiplex_Stats.csv",
                                      "Quality_Metrics.csv",
                                      "Top_Unknown_Barcodes.csv",
                                      "*fastq.gz_fastqc_data.txt"),
                                    collapse = " ",
                                    "'"))

# Download Cecret options
cecret_file_patterns <- c("software_versions.yml",
                          "_demix.tsv",
                          "_kraken2_report.txt",
                          "snp-dists.txt",
                          "_amplicon_depth.csv",
                          "nextclade.tsv",
                          ".stats.txt",
                          ".cov.txt",
                          ".depth.txt",
                          "cecret_results.csv")

cecret_file_download_command <- c("s3 cp", workflow_output_fp)
cecret_file_download_param <- c("--recursive", "--exclude '*'",
                                paste0("--include '*", cecret_file_patterns, collapse = " ", "'"))

# Download BCLConvert data
download_bclconvert_session <- paste0("down-bclconvert-", session_suffix)
submit_screen_job(message2display = "Downloading BCLConvert data",
                  ec2_login = ec2_hostname,
                  screen_session_name = download_bclconvert_session,
                  screen_log_fp = tmp_screen_fp,
                  command2run = paste("mkdir -p", paste0(data_output_fp, ";"),
                                      "aws",
                                      paste0(bcl_file_download_command, collapse = " "),
                                      data_output_fp,
                                      paste0(bcl_file_download_param, collapse = " "))
)

check_screen_job(message2display = "Checking BCLConvert download job",
                 ec2_login = ec2_hostname,
                 screen_session_name = download_bclconvert_session,
                 screen_log_fp = tmp_screen_fp)

# Download Cecret data
download_cecret_session <- paste0("down-cecret-", session_suffix)
submit_screen_job(message2display = "Downloading Cecret data",
                  ec2_login = ec2_hostname,
                  screen_session_name = download_cecret_session,
                  screen_log_fp = tmp_screen_fp,
                  command2run = paste("mkdir -p", paste0(data_output_fp, ";"),
                                      "aws",
                                      paste0(cecret_file_download_command, collapse = " "),
                                      data_output_fp,
                                      paste0(cecret_file_download_param, collapse = " "))
)

check_screen_job(message2display = "Checking Cecret download job",
                 ec2_login = ec2_hostname,
                 screen_session_name = download_cecret_session,
                 screen_log_fp = tmp_screen_fp)

# Download data from EC2 instance to local
run_in_terminal(paste("scp -r", paste0(ec2_hostname, ":", data_output_fp),
                      here())
)

# Download Nextflow config file for profile (use terminal because of proxy login issue)
# If SSH is not working, copy a nextflow.config from an older run and paste into data/processed_cecret as nextflow.config
run_in_terminal(paste("scp", paste0(ec2_hostname, ":~/.nextflow/config"),
                      here("data", "processed_cecret", "nextflow.config"))
)

# Tag intermediate processed files for deletion after 90 days
# Intermediate files are files that do not match the file patterns that are downloaded; exceptions are consensus and filtered fastq files for uploading
# Objects tagged with the tag set {Key=Fading,Value=90days} will be deleted by the bucket's lifecycle policy in 90 days
tagset <- "'TagSet=[{Key=Fading,Value=90days}]'"
tag_filename <- "s3_object_keys_2_tag.csv"

file_patterns_not_tagged <- c(cecret_file_patterns,
                              "/ivar_consensus/.*.consensus.fa",
                              "_filtered_R[12].fastq.gz")

aws_s3_cecret_intermediate_files <- system2("ssh", c("-tt", ec2_hostname,
                                                     shQuote(paste("aws s3 ls", workflow_output_fp, "--recursive",
                                                                   "| grep -ve", #reverse grep files with these patterns
                                                                   paste0(file_patterns_not_tagged, collapse = " -e ")),
                                                             type = "sh")),
                                            stdout = TRUE, stderr = TRUE) %>%
  head(-1)

if(identical(aws_s3_cecret_intermediate_files, character(0))) {
  message("\nNo intermediate files found for tagging!")
} else {

  # This is the list of files to tag for deletion
  cecret_intermediate_files <- aws_s3_cecret_intermediate_files %>%
    str_split("\\s+") %>%
    do.call("rbind", .) %>%
    as.data.frame() %>%
    `colnames<-`(c("date", "time", "bytes", "filename")) %>%
    select(bytes, filename) %>%
    mutate(bytes = as.numeric(bytes))

  cecret_intermediate_files %>%
    select(filename) %>%
    write_csv(here("data", "processed_cecret", tag_filename), col_names = FALSE)

  message("Uploading tag list to EC2")
  tag_s3_session <- paste0("nf-tag-s3-obj-", session_suffix)
  remote_tag_fp <- paste(tmp_screen_fp, tag_s3_session, sep = "/")
  mk_remote_dir(ec2_hostname, remote_tag_fp)

  run_in_terminal(paste("scp", here("data", "processed_cecret", tag_filename),
                        paste0(ec2_hostname, ":", remote_tag_fp))
  )

  # Copy the tag nextflow pipeline to EC2
  run_in_terminal(paste("scp", file.path(dirname(here()), "aux_files", "external_scripts", "nextflow", "tag_s3_objects.nf"),
                        paste0(ec2_hostname, ":", remote_tag_fp))
  )

  # Run nextflow tag objects pipeline
  submit_screen_job(message2display = "Tagging S3 objects for eventual removal",
                    ec2_login = ec2_hostname,
                    screen_session_name = tag_s3_session,
                    screen_log_fp = tmp_screen_fp,
                    command2run = paste("cd", remote_tag_fp,
                                        "nextflow run tag_s3_objects.nf",
                                        "-bucket-dir", paste0(s3_nextflow_work_bucket, "/tag_objects_", sample_type_acronym, "_", sequencing_date),
                                        #"-profile", demux_profile,
                                        #"--awscli_container", "",
                                        "--obj_key_samplesheet", tag_filename,
                                        "--s3_bucket", gsub("^s3://|/.*", "", s3_nextflow_output_bucket),
                                        "--obj_tagset", tagset)
  )

  # check_screen_job(message2display = "Checking S3 tagging job",
  #                  ec2_login = ec2_hostname,
  #                  screen_session_name = tag_s3_session,
  #                  screen_log_fp = tmp_screen_fp)
}

# Clean up environment
clean_tmp_session <- paste0("clean-tmp-", session_suffix)
submit_screen_job(message2display = "Cleaning up run from temporary folder",
                  ec2_login = ec2_hostname,
                  screen_session_name = clean_tmp_session,
                  screen_log_fp = tmp_screen_fp,
                  command2run = paste("rm -rf",
                                      paste0(ec2_tmp_fp, "/", session_suffix, "/;"),
                                      "echo 'Here are the files in the tmp directory:';",
                                      "ls", ec2_tmp_fp)
)

check_screen_job(message2display = "Checking delete job",
                 ec2_login = ec2_hostname,
                 screen_session_name = clean_tmp_session,
                 screen_log_fp = tmp_screen_fp)
