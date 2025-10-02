library(here)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)

#This Rscript submits the relevant jobs to Nextflow once the sequencing run has been uploaded

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

###############
# Load R config
###############

#this file needs to sit in a [aux_files/r_scripts/config] directory path above this project directory
tryCatch(
  {
    source(file.path(dirname(here()), "aux_files", "r_scripts", "config", "config_variables.R"))
  },
  error = function(e) {
    stop (simpleError("The config_variables.R file needs to sit in a [aux_files/r_scripts/config] directory path above this project directory"))
  }
)

############################
# Project specific variables
############################

#sequencing date of the run folder should match the RStudio project date
sequencing_date <- gsub("_.*", "", basename(here())) #YYYY-MM-DD

if(sequencing_date == "") {
  stop (simpleError(paste0("Please fill in the correct sequencing date or short project description in ", here("code"), "/4_run_nextflow.R")))
} else if (is.na(as.Date(sequencing_date, "%Y-%m-%d")) | nchar(sequencing_date) == 8) {
  stop (simpleError("Please enter the date into [sequencing_date] as YYYY-MM-DD"))
}

# If this sample sheet is missing, get it from AWS S3 bucket
sample_sheet_fn <- list.files(here("metadata", "munge"), pattern = "SampleSheet_v2.csv")

if(length(sample_sheet_fn) > 1) {
  stop(simpleError("There are more than 2 sample sheets detected!! Please delete the incorrect one"))
}

instrument_type <- gsub("^[0-9-]+_(MiSeq|NextSeq2000)_.*", "\\1", sample_sheet_fn)

# suffix for the screen session log names
session_suffix <- tolower(paste(instrument_type, sample_type_acronym, pathogen_acronym, basename(here()), sep = "-"))

sequencer_regex <- case_when(instrument_type == "MiSeq" ~ "M",
                             instrument_type == "NextSeq2000" ~ "VH")

intended_sequencing_folder_regex <- paste0(gsub("^..|-", "", sequencing_date), "_", sequencer_regex, "[0-9]*_[0-9]*_[0-9A-Z-]*")

# temporary directory to hold the screen log files and files for uploading
tmp_screen_path <- paste("~", ".tmp_screen", instrument_type, paste0(sample_type_acronym, "_", pathogen_acronym), basename(here()), sep = "/")

staging_path <- paste0(tmp_screen_path, "/staging/run_nextflow/")

data_output_path <- paste0(staging_path, "data/")

#####################################
# AWS and Nextflow specific variables
#####################################

# the Batch job definition can be created using the batch_create_nf_headnode.sh script
# the Docker image used by the job definition should be built from the Dockerfile in aux_files/docker_builds/nextflow/24.10.3
# this image needs to be uploaded to an online repository and the URI can be passed to the -i option of the script
nf_headnode_definition <- "nf-headnode"

# batch_create_nf_headnode.sh will run on the dedicated instance defined in ec2_hostname and requires aws-cli to be installed
nf_headnode_script <- paste(s3_aux_files_bucket, "external_scripts", "bash", "batch_create_nf_headnode.sh", sep = "/")

nf_config_fp <- paste(s3_aux_files_bucket, "external_scripts", "nextflow", nf_config_fn, sep = "/")

# S3 paths
nf_demux_bucket_path <- paste0(s3_nextflow_work_bucket, "/demux_", sample_type_acronym, "_", sequencing_date)
nf_demux_output_path <- paste(s3_fastq_bucket, sequencing_date, sample_type_acronym, prj_description, "processed_bclconvert", sep = "/")

nf_cecret_bucket_path <- paste0(s3_nextflow_work_bucket, "/cecret_", sample_type_acronym, "_", sequencing_date)
nf_cecret_output_path <- paste(s3_nextflow_output_bucket, "cecret", sample_type_acronym, paste0(sequencing_date, "_", prj_description), instrument_type, sep = "/")

nf_tag_s3_bucket_path <- paste0(s3_nextflow_work_bucket, "/tag_objects_", sample_type_acronym, "_", sequencing_date)

if(nchar(nf_base_container) == 0 | !is.character(nf_base_container)) {
  nf_base_container <- "staphb/ivar:1.4.1"
}

######################################
# Download and extract nextflow config
######################################

nextflow_profiles <- list(base = "phl2_main",
                          demux = "phl2_nfcore_demux",
                          cecret = tolower(paste0("phl2_cecret_", sample_type_acronym)))

nf_config <- system2("ssh", c("-tt", ec2_hostname,
                              shQuote(paste("aws s3 cp", nf_config_fp,
                                            staging_path,
                                            # add echo to buffer the cat output of the config file
                                            "&& echo && cat", paste0(staging_path, nf_config_fn), "&& echo"), type = "sh")),
                     stdout = TRUE, stderr = TRUE) %>%
  data.frame(stdout = .) %>%
  # system2 stdout appends a CR at the end of each line that needs to be removed
  mutate(stdout = gsub("\r$", "", stdout)) %>%
  filter(stdout != "",
         !grepl("Completed [0-9]|^Connection to .* closed.", stdout))

profile_write <- data.frame(to_write = character())
num_profiles_2_record <- length(nextflow_profiles)
open_brack_count <- 0
declared_profile <- FALSE
record_settings <- FALSE
record_profile <- FALSE

for(current_line in 1:nrow(nf_config)) {

  if(grepl("aws \\{$|docker \\{$", nf_config[current_line,])) {
    open_brack_count <- 0
    record_settings <- TRUE
  }

  if(grepl(paste0(nextflow_profiles, " \\{$", collapse = "|"), nf_config[current_line,])) {
    open_brack_count <- 0
    record_profile <- TRUE
    num_profiles_2_record <- num_profiles_2_record - 1

    if(!declared_profile) {
      declared_profile <- TRUE
      profile_write <- rbind(profile_write, data.frame(to_write = "profiles {"))
    }
  }

  if(grepl("\\{", nf_config[current_line,])) {
    open_brack_count <- open_brack_count + 1
  }

  if(grepl("\\}", nf_config[current_line,])) {
    open_brack_count <- open_brack_count - 1
  }

  if(record_settings | record_profile) {
    if(open_brack_count > 0) {
      profile_write <- rbind(profile_write, data.frame(to_write = nf_config[current_line,]))
    } else if(open_brack_count == 0) {
      profile_write <- rbind(profile_write, data.frame(to_write = "}"))

      if(declared_profile & num_profiles_2_record == 0) {
        profile_write <- rbind(profile_write, data.frame(to_write = "}"))
        declared_profile <- FALSE
      }

      record_settings <- FALSE
      record_profile <- FALSE
    }
  }
}

nf_config_settings <- profile_write %>%
  filter(!grepl("^\\s+}|^}$", to_write)) %>%
  mutate(stripped = str_trim(to_write),
         config_indents = get_indents(to_write),
         config_indents = factor(config_indents, levels = str_sort(unique(config_indents))),
         config_indents = paste0("config", as.numeric(config_indents)),
         row_id = row_number()) %>%
  select(-to_write) %>%
  pivot_wider(names_from = "config_indents", values_from = "stripped", values_fill = NA) %>%
  select(-row_id) %>%
  t() %>%
  as.data.frame() %>%
  fill(everything(), .direction = "down") %>%
  t() %>%
  as.data.frame() %>%
  fill(starts_with("config"), .direction = "down")

aws_region <- nf_config_settings %>%
  filter(grepl("aws", config1),
         grepl("region", config2)) %>%
  select(config2) %>%
  mutate(config2 = gsub("^.*'([/A-Za-z0-9_-]+)'$", "\\1", config2)) %>%
  pull()

ami_aws_cli <- nf_config_settings %>%
  filter(grepl("aws", config1),
         grepl("batch", config2),
         grepl("cliPath", config3)) %>%
  select(config3) %>%
  mutate(config3 = gsub("^.*'([/A-Za-z0-9_-]+)'$", "\\1", config3)) %>%
  pull()

jq_list <- nf_config_settings %>%
  filter(grepl("profiles", config1),
         grepl("process", config3),
         grepl("queue", config4)) %>%
  mutate(config4 = gsub("^.*'([/A-Za-z0-9_-]+)'$", "\\1", config4))

nf_demux_jq <- jq_list %>%
  filter(grepl(nextflow_profiles$demux, config2)) %>%
  select(config4) %>%
  pull()

nf_cecret_jq  <- jq_list %>%
  filter(grepl(nextflow_profiles$cecret, config2)) %>%
  select(config4) %>%
  pull()

nf_general_jq  <- jq_list %>%
  filter(grepl(nextflow_profiles$base, config2)) %>%
  select(config4) %>%
  pull()

nextflow_config_fp <- here("data", paste0(sequencing_date, "_nextflow.config"))

profile_write %>%
  write_csv(file = nextflow_config_fp, col_names = FALSE, quote = "none")

nf_profile_software_version <- profile_write %>%
  filter(grepl("kraken2_db|primer_set", to_write)) %>%
  separate(col = "to_write", into = c("software", "version"), sep = " = ") %>%
  mutate(pipeline = "UPHL-BioNGS-Cecret",
         version = gsub("'", "", version),
         version = gsub("/$", "", version),
         version = gsub(".*/", "", version)) %>%
  select(pipeline, software, version)

software_version_fp <- here("data", paste0(sequencing_date, "_software_version.csv"))

# Don't overwrite the software version file, if it exists
if(!file.exists(software_version_fp)) {
  nf_profile_software_version %>%
    write_csv(file = software_version_fp)
}

#######################################
# Demultiplex run and check FASTQ files
#######################################

nf_demux_samplesheet_fp <- paste(s3_run_bucket, sequencing_date,
                                   tolower(paste(sequencing_date, instrument_type, sample_type_acronym, prj_description, "nf_demux_samplesheet.csv", sep = "_")), sep = "/")

demux_session <- paste0("nf-demux-", session_suffix)
nf_headnode_screen_job(batch_job_queue = nf_demux_jq,
                       nf_bucket_dir = nf_demux_bucket_path,
                       message2display = "Demultiplexing with BCLConvert",
                       screen_session_name = demux_session,
                       command2run = paste("nextflow run nf-core/demultiplex",
                                           "-c ~/.nextflow/config",
                                           "-profile", nextflow_profiles$demux,
                                           "-bucket-dir", nf_demux_bucket_path,
                                           "-resume",
                                           "-r 1.3.2",
                                           "--input", nf_demux_samplesheet_fp,
                                           "--outdir", nf_demux_output_path)
)

check_screen_job(message2display = "Checking BCLConvert job",
                 screen_session_name = demux_session)

# Checking the demultiplexing results
message("Checking fastq file sizes...")
aws_s3_fastq_files_bucketdir <- system2("ssh", c("-tt", ec2_hostname,
                                                 shQuote(paste("aws s3 ls", nf_demux_bucket_path, "--recursive",
                                                               "| grep 'R1_001.fastq.gz$'"), type = "sh")),
                                        stdout = TRUE, stderr = TRUE) %>%
  head(-1)

fastq_file_sizes_bucketdir <- aws_s3_fastq_files_bucketdir %>%
  str_split("\\s+") %>%
  do.call("rbind", .) %>%
  as.data.frame() %>%
  `colnames<-`(c("date", "time", "bytes", "filename")) %>%
  select(bytes, filename) %>%
  filter(!grepl("/Alignment_|/Fastq/|Undetermined", filename)) %>%
  mutate(filename = gsub(".*/|_.*", "", filename),
         bytes = as.numeric(bytes)) %>%
  #sometimes the Undetermined file from the sequencer gets copied over; remove these samples
  filter(!grepl("^GenericSampleID", filename),
         bytes <= 23) %>%
  group_by(filename) %>%
  mutate(file_num = n()) %>%
  filter(!(file_num == 2 & grepl("^Undetermined", filename) & bytes == suppressWarnings(max(bytes)))) %>%
  ungroup()

if(nrow(fastq_file_sizes_bucketdir) > 0) {

  rm_work_bclconvert_session <- paste0("rm-bclconvert-workdir-", session_suffix)
  submit_screen_job(message2display = "Removing empty FASTQ files from working bucket",
                    screen_session_name = rm_work_bclconvert_session,
                    command2run = paste("aws s3 rm", nf_demux_bucket_path,
                                        "--recursive",
                                        "--exclude '*'",
                                        paste0("--include '*", fastq_file_sizes_bucketdir$filename, "*.fastq.gz*'",
                                               collapse = " "))
  )

  check_screen_job(message2display = "Checking remove workdir FASTQ job",
                   screen_session_name = rm_work_bclconvert_session)

  rm_out_bclconvert_session <- paste0("rm-bclconvert-outdir-", session_suffix)
  submit_screen_job(message2display = "Removing recently demultiplexed FASTQ files in outdir",
                    screen_session_name = rm_out_bclconvert_session,
                    command2run = paste("aws s3 rm", nf_demux_output_path,
                                        "--recursive",
                                        "--exclude '*'",
                                        paste0("--include '", intended_sequencing_folder_regex, "'"))
  )

  check_screen_job(message2display = "Checking remove outdir FASTQ job",
                   screen_session_name = rm_out_bclconvert_session)

  message("\nSampleIDs with no reads: ", paste0("\"", paste0(fastq_file_sizes_bucketdir$filename, collapse = "\", \""), "\"\n"))

  stop(simpleError(paste("Number of samples that had no reads:", nrow(fastq_file_sizes_bucketdir),
                         "\n\nThese samples are held in", nf_demux_bucket_path,
                         "\nand could not be moved to", nf_demux_output_path,
                         "\n\nIf the number of samples that have no reads corresponds to the total number of samples for the sequencing run,",
                         "\nthe wrong barcode index may have been used. Please check that the correct IDT set are used for the barcodes",
                         "\n\nOtherwise, some samples may not have reads even with the correct barcodes used (negative controls or just no DNA present in the sample)",
                         "\nFor these samples, please remove them from the analysis by adding their SampleID to the vector sample_w_empty_reads in *_QC_Report.Rmd",
                         "\n\nRegenerate the SampleSheet with the generate_metadata_.* Rscript and then rerun this script to re-demultiplex")))
}

aws_s3_fastq_files <- system2("ssh", c("-tt", ec2_hostname,
                                       shQuote(paste("aws s3 ls", nf_demux_output_path, "--recursive",
                                                     "| grep 'R1_001.fastq.gz$'"), type = "sh")),
                              stdout = TRUE, stderr = TRUE) %>%
  head(-1)

if(identical(aws_s3_fastq_files, character(0))) {
  stop(simpleError(paste("\nSomething went wrong! No FASTQ files found in", nf_demux_output_path)))
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
  filter(grepl(intended_sequencing_folder_regex, sequencing_folder)) %>%
  #sometimes the Undetermined file from the sequencer gets copied over; remove these samples
  filter(!grepl("^GenericSampleID", filename)) %>%
  group_by(filename) %>%
  mutate(file_num = n()) %>%
  filter(!(file_num == 2 & grepl("^Undetermined", filename) & bytes == suppressWarnings(max(bytes)))) %>%
  ungroup()

if(nrow(fastq_file_sizes) == 0) {
  stop(simpleError(paste0("\nThere were no FastQ files found at path ", nf_demux_output_path,
                          "\nCheck to see if there was an issue with the demultiplexing of the run\n")))
}

if(any(fastq_file_sizes$bytes <= 23)) {

  sample_id_no_reads <- fastq_file_sizes %>%
    mutate(filename = gsub("_.*", "", filename)) %>%
    filter(!grepl("Undetermined", filename),
           bytes <= 23)

  stop(simpleError(paste("\nThese SampleIDs had no reads in", nf_demux_output_path, "\n",
                         paste0("\"", paste0(sample_id_no_reads$filename, collapse = "\", \""), "\""),
                         "\n\nPlease investigate!!")))
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
                    screen_session_name = undetermined_mv_session,
                    command2run = paste("aws s3 mv",
                                        paste(nf_demux_output_path, instrument_run_id, sep = "/"),
                                        paste(nf_demux_output_path, "Undetermined", instrument_run_id, sep = "/"),
                                        "--recursive",
                                        "--exclude '*'",
                                        "--include '*Undetermined_S0_*_001.fastq.gz'")
  )

  check_screen_job(message2display = "Checking Undetermined move job",
                   screen_session_name = undetermined_mv_session)
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

fastq_path <- paste(nf_demux_output_path, instrument_run_id, sep = "/")

#########################################################
# Update the Nextclade dataset for SC2 lineage assignment
#########################################################

if (nextclade_dataset_version != "") {
  nextclade_tag <- paste("--tag", nextclade_dataset_version)
} else {
  nextclade_tag <- ""
}

update_nextclade_session <- paste0("update-nextclade-", session_suffix)

submit_screen_job(message2display = "Downloading Nextclade SARS-CoV-2 data",
                  screen_session_name = update_nextclade_session,
                  command2run = paste("cd", paste0(staging_path, ";"),
                                      "wget -q https://github.com/nextstrain/nextclade/releases/latest/download/nextclade-x86_64-unknown-linux-gnu -O ./nextclade;",
                                      "chmod +x ./nextclade;",
                                      "./nextclade --version;",
                                      "./nextclade dataset get --name sars-cov-2", nextclade_tag, "--output-zip ./sars.zip;",
                                      "aws s3 cp", paste0(s3_reference_bucket, "/nextclade/sars.zip"), "./sars_ref_file_2_check.zip;",
                                      # compare the two zip files, if they are different, upload the latest version to s3
                                      "cmp ./sars.zip ./sars_ref_file_2_check.zip",
                                      "&& echo 'Zip files are the same. No need to upload'",
                                      "|| aws s3 cp ./sars.zip", paste0(s3_reference_bucket, "/nextclade/sars.zip;"),
                                      "rm -v ./*.zip ./nextclade")
)

check_screen_job(message2display = "Checking Nextclade download job",
                 screen_session_name = update_nextclade_session)

############
# Run Cecret
############

cecret_session <- paste0("nf-cecret-", session_suffix)
nf_headnode_screen_job(batch_job_queue = nf_cecret_jq,
                       nf_bucket_dir = nf_cecret_bucket_path,
                       message2display = "Processing data through Cecret pipeline",
                       screen_session_name = cecret_session,
                       command2run = paste("nextflow run UPHL-BioNGS/Cecret",
                                           "-profile", nextflow_profiles$cecret,
                                           "-bucket-dir", nf_cecret_bucket_path,
                                           "-r", cecret_version,
                                           "-resume",
                                           "--reads", fastq_path,
                                           "--outdir", paste(nf_cecret_output_path, "processed_cecret", sep = "/"))
)

check_screen_job(message2display = "Checking Cecret job",
                 screen_session_name = cecret_session)

##################
# Download results
##################

# BCLConvert file patterns
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

# Cecret file patterns
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

cecret_file_download_param <- c("--recursive", "--exclude '*'",
                                paste0("--include '*", cecret_file_patterns, collapse = " ", "'"))

# Download BCLConvert data
download_bclconvert_session <- paste0("down-bclconvert-", session_suffix)
submit_screen_job(message2display = "Downloading BCLConvert data",
                  screen_session_name = download_bclconvert_session,
                  command2run = paste("aws s3 cp", dirname(nf_demux_output_path),
                                      data_output_path,
                                      paste0(bcl_file_download_param, collapse = " "))
)

check_screen_job(message2display = "Checking BCLConvert download job",
                 screen_session_name = download_bclconvert_session)

# Download Cecret data
download_cecret_session <- paste0("down-cecret-", session_suffix)
submit_screen_job(message2display = "Downloading Cecret data",
                  screen_session_name = download_cecret_session,
                  command2run = paste("aws s3 cp", nf_cecret_output_path,
                                      data_output_path,
                                      paste0(cecret_file_download_param, collapse = " "))
)

check_screen_job(message2display = "Checking Cecret download job",
                 screen_session_name = download_cecret_session)

# Download data from EC2 instance to local
run_in_terminal(paste("scp -r", paste0(ec2_hostname, ":", data_output_path),
                      here())
)

# Tag intermediate processed files for deletion after 90 days
# Intermediate files are files that do not match the file patterns that are downloaded; exceptions are consensus and filtered fastq files for uploading
# Objects tagged with the tag set {Key=Fading,Value=90days} will be deleted by the bucket's lifecycle policy in 90 days
tag_keys <- paste("Fading", sep = ",") #keys and values are strings of comma separated characters
tag_values <- paste("90days", sep = ",")

if(length(str_split_1(tag_keys, ",")) < length(str_split_1(tag_values, ","))) {
  stop(simpleError("Number of values cannot exceed number of keys for tagging"))
}
if(any(str_split_1(tag_keys, ",") == "")) {
  stop(simpleError("Cannot have an empty key"))
}

tag_filename <- "s3_object_keys_2_tag.csv"
nf_tag_s3_samplesheet_fp <- paste(nf_tag_s3_bucket_path, tag_filename, sep = "/")

file_patterns_not_tagged <- c(gsub("\\.", "\\\\.", cecret_file_patterns),
                              "/ivar_consensus/.*\\.consensus\\.fa",
                              "_filtered_R[12]\\.fastq\\.gz")

aws_s3_cecret_intermediate_files <- system2("ssh", c("-tt", ec2_hostname,
                                                     shQuote(paste("aws s3 ls", nf_cecret_output_path, "--recursive",
                                                                   "| grep -ve", #reverse grep files with these patterns
                                                                   paste0("'", paste0(file_patterns_not_tagged, collapse = "' -e '"), "'")),
                                                             type = "sh")),
                                            stdout = TRUE, stderr = TRUE) %>%
  head(-1)

if(!identical(aws_s3_cecret_intermediate_files, character(0))) {

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

  upload_tag_s3_session <- paste0("up-tag-s3-", session_suffix)

  run_in_terminal(paste("scp", here("data", "processed_cecret", tag_filename),
                        paste0(ec2_hostname, ":", staging_path))
  )

  # Upload s3 tag samplesheet
  submit_screen_job(message2display = "Uploading tag sheet to S3",
                    screen_session_name = upload_tag_s3_session,
                    command2run = paste("aws s3 cp",
                                        paste0(staging_path, tag_filename),
                                        nf_tag_s3_samplesheet_fp)
  )

  check_screen_job(message2display = "Checking tag sheet upload job",
                   screen_session_name = upload_tag_s3_session)

  nf_tag_s3_session <- paste0("nf-tag-s3-obj-", session_suffix)
  nf_headnode_screen_job(batch_job_queue = nf_general_jq,
                         nf_bucket_dir = nf_tag_s3_bucket_path,
                         message2display = "Tagging S3 objects for eventual removal",
                         screen_session_name = nf_tag_s3_session,
                         command2run = paste(ami_aws_cli, "s3 cp", paste(s3_aux_files_bucket, "external_scripts", "nextflow", "tag_s3_objects.nf", sep = "/"), ".;",
                                             "nextflow run tag_s3_objects.nf",
                                             "-bucket-dir", nf_tag_s3_bucket_path,
                                             "-profile", nextflow_profiles$base,
                                             "--use_container", nf_base_container,
                                             "--aws_region", aws_region,
                                             "--clipath", ami_aws_cli, #provide an ami clipath if running through AWS Batch
                                             "--obj_key_samplesheet", nf_tag_s3_samplesheet_fp,
                                             "--s3_bucket", gsub("^s3://|/.*", "", s3_nextflow_output_bucket),
                                             "--tag_keys", tag_keys,
                                             "--tag_values", tag_values)
  )

  # check_screen_job(message2display = "Checking S3 tagging job",
  #                  screen_session_name = nf_tag_s3_session)
}

# Clean up environment
clean_tmp_session <- paste0("clean-tmp-", session_suffix)
submit_screen_job(message2display = "Cleaning up run from temporary folder",
                  screen_session_name = clean_tmp_session,
                  command2run = paste("rm -rf",
                                      paste0(tmp_screen_path, "/staging/*;"),
                                      "echo 'Here are the files in the tmp directory:';",
                                      "ls", paste0(tmp_screen_path, "/staging/"))
)

check_screen_job(message2display = "Checking delete job",
                 screen_session_name = clean_tmp_session)

message("\nRscript finished successfully!")
