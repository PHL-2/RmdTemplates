library(here)
library(readr)
library(tidyverse)

lastest_seqsender_version_tested <- "v1.2.7"

#load constant variables
tryCatch(
  {
    source(file.path(dirname(here()), "aux_files", "r_scripts", "config", "config_variables.R"))
  },
  error = function(e) {
    stop (simpleError("The config_variables.R file needs to sit in a [aux_files/r_scripts/config] directory path above this project directory"))
  }
)

#load functions
tryCatch(
  {
    source(file.path(dirname(here()), "aux_files", "r_scripts", "functions", "R_all_functions_v3.R"))
  },
  error = function(e) {
    stop (simpleError("The R_all_functions_v3.R file needs to sit in a [aux_files/r_scripts/functions] directory path above this project directory"))
  }
)

ssh_seqsender_cmd <- function(pos_argument,
                              ssh_destination = ec2_hostname,
                              organism = organism_flag,
                              submission_dir = submission_path,
                              submission_name = input_submission_name,
                              config = config2use,
                              metadata = ec2_seqsender_meta_fp,
                              fasta = consensus_fasta_ul_fp,
                              is_test = test_flag) {

  docker_prepend <- paste0("docker run --rm -u $(id -u):$(id -g) -v $(pwd):/data cdcgov/seqsender:",
                           lastest_seqsender_version_tested)

  additional_flags <- paste(organism,
                            "--config_file", config,
                            "--metadata_file", metadata)

  success_msg <- unlist(lapply(c("BIOSAMPLE", "SRA", "GENBANK", "GISAID"),
                               function(x) c(paste("Creating submission files for", x),
                                             paste0("Files are stored at: /data/", paste(submission_dir, submission_name, "submission_files", x, sep = "/")))))

  if(grepl("--gisaid", pos_argument)) {
    additional_flags <- paste(additional_flags, "--fasta_file", fasta)
    success_msg <- success_msg[grepl("GISAID", success_msg)]
  }
  if(grepl("--biosample", pos_argument)) {
    success_msg <- success_msg[grepl("BIOSAMPLE", success_msg)]
    submission_msg <- paste0("Submission name: ", submission_name, "_BIOSAMPLE")
  }
  if(grepl("--sra", pos_argument)) {
    success_msg <- success_msg[grepl("SRA", success_msg)]
    submission_msg <- paste0("Submission name: ", submission_name, "_SRA")
  }
  if(grepl("--genbank", pos_argument)) {
    additional_flags <- paste(additional_flags, "--fasta_file", fasta)
    success_msg <- success_msg[grepl("GENBANK", success_msg)]
    submission_msg <- paste0("Submission name: ", submission_name, "_GENBANK")
  }
  if(grepl("^submit", pos_argument)) {
    additional_flags <- paste(additional_flags, is_test)
    success_msg <- c(success_msg,
                     "Connecting to NCBI FTP Server",
                     submission_msg,
                     paste0("Submitting '", submission_name, "'"))
  }
  if (pos_argument == "submission_status") {
    additional_flags <- ""

    success_msg <- c("Checking Submissions:",
                     paste("Submission:", submission_name),
                     "Downloading report.xml",
                     "Updating submissions complete.")
  }

  check_err_msg <- system2("ssh", c("-tt", ssh_destination,
                                    shQuote(paste(docker_prepend,
                                                  "seqsender.py", pos_argument,
                                                  "--submission_dir", gsub("/home/ec2-user/", "", submission_dir), #test
                                                  "--submission_name", submission_name,
                                                  additional_flags), type = "sh")),
                           stdout = TRUE, stderr = TRUE) %>%
    gsub("\\\r$", "", .)

  warning(simpleWarning(paste(check_err_msg, collapse = "\n")))

  if(!all(success_msg %in% check_err_msg)) {
    stop(simpleError(paste("\nSeqsender process", pos_argument, "failed!")))
  }
}

ssh_command_check <- function(stdoutput) {
  if(any(grepl("Could not resolve hostname|Operation timed out|BrokenPipeError|No such file or directory", stdoutput))) {
    stop(simpleError("SSH command to EC2 command failed"))
  }
}

#################################
# Load project specific variables
#################################

project_name <- basename(here())
aux_seqsender_fp <- file.path(dirname(here()), "aux_files", "data_submission", "seqsender")
submission_log_fp <- file.path(dirname(here()), "submission_log.csv")
fasta_seqsender_meta_fp <- here("upload", "seqsender", paste(project_name, "PHL2_seqsender_fasta_upload_orig.csv", sep = "_"))
fasta_appended_seqsender_meta_fp <- here("upload", "seqsender", paste(project_name, "PHL2_seqsender_fasta_upload.csv", sep = "_"))
fastq_seqsender_meta_fp <- here("upload", "seqsender", paste(project_name, "PHL2_seqsender_fastq_upload_orig.csv", sep = "_"))
fastq_appended_seqsender_meta_fp <- here("upload", "seqsender", paste(project_name, "PHL2_seqsender_fastq_upload.csv", sep = "_"))

#path of local seqsender config file
seqsender_config_fp <- list.files(aux_seqsender_fp, pattern = "seqsender_ww_sc2_login_config.yaml", full.names = TRUE)

#sequencing date of the run folder should match the RStudio project date
sequencing_date <- gsub("_.*", "", basename(here())) #YYYY-MM-DD

#if this sample sheet file is missing, download it manually from AWS S3 bucket
sample_sheet_fn <- list.files(here("metadata", "munge"), pattern = "SampleSheet_v2.csv")

#get sequencer of the run
instrument_type <- gsub("^[0-9-]*_(MiSeq|NextSeq2000)_.*", "\\1", sample_sheet_fn)

#get sample type
sample_type_acronym <- gsub(paste0("^[0-9-]*_", instrument_type, "_|_.*"), "", sample_sheet_fn)

authors <- gsub(", .*|,.*", "", bioinformatician_name)

#checks
if(sequencing_date == "") {
  stop (simpleError(paste0("Please fill in the correct sequencing date or short project description")))
} else if (is.na(as.Date(sequencing_date, "%Y-%m-%d")) | nchar(sequencing_date) == 8) {
  stop (simpleError("Please enter the date into [sequencing_date] as YYYY-MM-DD"))
}

if(length(sample_sheet_fn) > 1) {
  stop(simpleError("There are more than 2 sample sheets detected!! Please delete the incorrect one"))
} else if(length(sample_sheet_fn) == 0) {
  stop(simpleError("Sample sheet with suffix '_SampleSheet_v2.csv' is missing. Download this file manually from AWS S3 bucket"))
}

if(bioinformatician_name == "") {
  stop (simpleError(paste0("Please fill in bioinformatician_name in config_variables.R")))
}

####################################
# Define AWS variables and filepaths
####################################

ec2_home_fp <- "/home/ec2-user"
seqsender_upload_tmp_dir <- "tmp_data_ul"
workflow_output_fp <- paste(s3_nextflow_output_bucket, "cecret", sample_type_acronym, project_name, instrument_type, sep = "/")
submission_path <- paste(seqsender_upload_tmp_dir, sample_type_acronym, sep = "/")
input_submission_name <- paste(project_name, sample_type_acronym, sep = "_")
ec2_upload_tmp_fp <- paste(ec2_home_fp, submission_path, input_submission_name, sep = "/")
consensus_fasta_ul_fp <- paste(submission_path, input_submission_name, "processed_cecret", "merged_consensus",
                               paste0(project_name, "_PHL2_combined.fasta"), sep = "/")
consensus_fasta_fp <- paste(ec2_home_fp, consensus_fasta_ul_fp, sep = "/")
fasta_ec2_seqsender_meta_fp <- paste(submission_path, input_submission_name, "seqsender",
                                     paste(project_name, "PHL2_seqsender_fasta_upload.csv", sep = "_"), sep = "/")
fastq_ec2_seqsender_meta_fp <- paste(submission_path, input_submission_name, "seqsender",
                                     paste(project_name, "PHL2_seqsender_fastq_upload.csv", sep = "_"), sep = "/")

if(test_upload) {
  test_flag <- "--test"
  gs_client_id <- "TEST-EA76875B00C3"
  local_gs_log <- here("gisaid", paste0(project_name, "_gisaid_test.log"))
  BioProject <- bio_prj_id_test
} else {
  test_flag <- ""
  gs_client_id <- gisaid_client_id
  local_gs_log <- here("gisaid", paste0(project_name, "_gisaid.log"))
  BioProject <- bio_prj_id_ww_sc2
}

organism_flag <- "--organism COV"
config2use <- paste(submission_path, input_submission_name, "seqsender", basename(seqsender_config_fp), sep = "/")
ec2_gs_log <- paste0(ec2_upload_tmp_fp, "/submission_files/GISAID/", basename(local_gs_log))

###############################################################
# Check submission log to see if run has already been submitted
###############################################################

submission_log <- read_csv(submission_log_fp)

submission_log_test <- submission_log %>%
  filter(Submission_Name == input_submission_name,
         Submission_Type == "TEST")

submission_log_prod <- read_csv(submission_log_fp) %>%
  filter(Submission_Name == input_submission_name,
         Submission_Type == "PRODUCTION")

if(nrow(submission_log_prod) > 0 & !test_upload) {
  stop(simpleError("This run has already been submitted by seqsender. Check NCBI to see submission status"))
}

if(nrow(submission_log_test) == 0 & !test_upload) {
  stop(simpleError(paste0("This run has not been submitted as a TEST yet. Please do this first by setting 'test_upload' to TRUE\n",
                          "A completed TEST run ensures that all the components in this chunk run successfully\n",
                          "Once a TEST run has completed, rerun this chunk with 'test_upload' set to FALSE for the real submission")))
}

#######################################
# Read in fasta and fastq seqsender csv
#######################################
fasta_seqsender <- read_csv(fasta_seqsender_meta_fp) %>%
  filter(!sequence_name %in% failed_sample_fasta,
         !is.na(`bs-sample_name`)) %>%
  mutate(authors = authors,
         bioproject = BioProject)

fastq_seqsender <- read_csv(fastq_seqsender_meta_fp) %>%
  filter(!sequence_name %in% failed_sample_fasta,
         !is.na(`bs-sample_name`)) %>%
  mutate(authors = authors,
         bioproject = BioProject)

seqsender_merge_check_cols <- c("gs-covv_subm_sample_id", "bs-isolate", "bs-ww_surv_system_sample_id")

seqsender_check <- fasta_seqsender %>%
  select(fasta_name = "sequence_name",
         fasta_processing = "bs-specimen_processing_id",
         all_of(seqsender_merge_check_cols)) %>%
  merge(select(fastq_seqsender,
               fastq_name = "sequence_name",
               fastq_processing = "bs-specimen_processing_id",
               all_of(seqsender_merge_check_cols)),
        by = seqsender_merge_check_cols, all = TRUE)

concentration_cols <- c("bs-ww_surv_target_1_conc", "bs-ww_surv_target_2_conc")

conc_check <- fasta_seqsender %>%
  rbind(fastq_seqsender) %>%
  select(any_of(concentration_cols)) %>%
  pivot_longer(cols = concentration_cols) %>%
  filter(is.na(value) | value == 0)

if(nrow(fasta_seqsender) < 1) {
  stop(simpleError("No samples passed QC for data submission. Don't have to submit this run"))
}
if(!all(grepl("Merged", unlist(fasta_seqsender[, c("sequence_name", "gs-sample_name", "bs-sample_name", "bs-specimen_processing_id", "bs-gisaid_virus_name", "sra-sample_name")])))) {
  stop(simpleError("Some of the samples in the seqsender metadata sheet for GISAID submission are not the merged sample"))
}
if(any(grepl("Merged", unlist(fastq_seqsender[, c("sra-file_1", "sra-file_2")])))) {
  stop(simpleError("The fastq seqsender metadata file has a Merged sample in it. Samples to be submitted should be individual replicates"))
}
if(any(is.na(seqsender_check$fastq_name))) {
  stop(simpleError("The fasta seqsender file does not match the fastq seqsender file"))
}
if(nrow(conc_check) > 0) {
  stop(simpleError("Do not continue to submit!! Samples are either missing ddPCR concentration data or have 0 for concentration. Investigate first"))
}

fasta_seqsender %>%
  write_csv(fasta_appended_seqsender_meta_fp, na = "")

fastq_seqsender %>%
  write_csv(fastq_appended_seqsender_meta_fp, na = "")

######################
# Run seqsender on EC2
######################

cecret_file_download_command <- c("aws s3 cp", workflow_output_fp)
cecret_file_download_param <- c("--recursive",
                                "--exclude '*'",
                                #download FASTA from ivar_consensus/ instead of consensus/ because file is not folded
                                paste0("--include '*/ivar_consensus/", fasta_seqsender$sequence_name, "*.consensus.fa'"),
                                paste0("--include '*/filter/", fastq_seqsender$sequence_name, "*_filtered_R[12].fastq.gz'"))

aws_s3_cecret_download <- system2("ssh", c("-tt", ec2_hostname,
                                           shQuote(paste0("rm -rf ", ec2_upload_tmp_fp, "; ",
                                                          paste(c(cecret_file_download_command,
                                                                  ec2_upload_tmp_fp,
                                                                  cecret_file_download_param), collapse = " "), "; ",
                                                          "mkdir -p ", ec2_upload_tmp_fp, "/raw_reads; ",
                                                          "find ", ec2_upload_tmp_fp, " -name '*_filtered_R[12].fastq.gz' -exec mv -t ",
                                                          ec2_upload_tmp_fp, "/raw_reads/ {} +"),
                                                   type = "sh")),
                                  stdout = TRUE, stderr = TRUE)

ssh_command_check(aws_s3_cecret_download)

create_concat_consensus <- system2("ssh", c("-tt", ec2_hostname,
                                            shQuote(paste0("mkdir -p ", dirname(consensus_fasta_fp), "; ",
                                                           "sed -E 's/^N+|N+$|Consensus_|\\.consensus_threshold.*//g' ",
                                                           ec2_upload_tmp_fp, "/processed_cecret/ivar_consensus/* > ",
                                                           consensus_fasta_fp, ";"), type = "sh")),
                                   stdout = TRUE, stderr = TRUE)

dir.create(here("upload", "fasta"))

download_consensus <- system2("scp", c(paste0(ec2_hostname, ":", consensus_fasta_fp),
                                       here("upload", "fasta/")),
                              stdout = TRUE, stderr = TRUE)

ssh_command_check(c(create_concat_consensus, download_consensus))

upload_seqsender_files <- system2("scp", c(fasta_appended_seqsender_meta_fp,
                                           fastq_appended_seqsender_meta_fp,
                                           seqsender_config_fp,
                                           paste0(ec2_hostname, ":", ec2_upload_tmp_fp, "/seqsender")))

ssh_command_check(upload_seqsender_files)

ssh_seqsender_cmd("prep --gisaid", metadata = fasta_ec2_seqsender_meta_fp)

dir.create(here("gisaid"))

# check the seqsender file for GISAID accession numbers
# if the seqsender file already has legit GISAID accession numbers, don't resubmit because GISAID already has these samples
# have to manually submit GISAID files because the password field in seqsender does not accept all special characters
if(all(grepl("^EPI_ISL_0$", fasta_seqsender$`bs-gisaid_accession`)) | all(is.na(fasta_seqsender$`bs-gisaid_accession`))) {

  # clear these GISAID log files if exist
  file.remove(local_gs_log)

  # submit files to gisaid
  # covCLI is a program that can be downloaded from gisaid.org; add it to $PATH
  run_in_terminal(paste("ssh -tt", ec2_hostname,
                        shQuote(paste("covCLI", "upload",
                                      "--clientid", gs_client_id,
                                      "--username", gisaid_username,
                                      "--metadata", paste(ec2_upload_tmp_fp, "submission_files", "GISAID", "metadata.csv", sep = "/"),
                                      "--fasta", paste(ec2_upload_tmp_fp, "submission_files", "GISAID", "sequence.fsa", sep = "/"),
                                      "--frameshift", "catch_novel",
                                      "--log", ec2_gs_log), type = "sh")))

  download_gisaid_load <- system2("scp", c(paste0(ec2_hostname, ":", ec2_gs_log), local_gs_log))

  ssh_command_check(download_gisaid_load)
}

# read in the GISAID log file
# if upload was successful, grab the GISAID accessions for the BioSample submissions
read_gisaid_accessions <- read_csv(local_gs_log, col_names = c("code", "msg", "timestamp"))

gisaid_errors <- read_gisaid_accessions %>%
  filter(grepl("error", code, ignore.case = TRUE))

if(nrow(gisaid_errors) > 0) {
  stop(simpleError(paste0("Something went wrong! GISAID upload failed. Did you use the wrong credentials?\n",
                          "This error may also occur if these samples have been already uploaded\n",
                          "Depending on following error message, you may be able to just rerun this chunk to resubmit:\n",
                          gisaid_errors)))
}

gisaid_success <- read_gisaid_accessions %>%
  filter(grepl("epi_isl_id", code))

if(nrow(gisaid_success) != nrow(fasta_seqsender)) {
  stop(simpleError("Some GISAID samples failed to upload."))
}

gisaid_accessions <- gisaid_success %>%
  mutate(msg = gsub("^msg: ", "", msg)) %>%
  separate(msg, into = c("bs-gisaid_virus_name", "bs-gisaid_accession"), sep = "; ", extra = "merge") %>%
  #dont need GenBank note because not submitting there
  #mutate(`src-note` = paste0("GISAID virus name: ", `bs-gisaid_virus_name`, "; GISAID accession: ", `bs-gisaid_accession`)) %>%
  select(`bs-gisaid_virus_name`, `bs-gisaid_accession`)

if(!test_upload) {
  #if there is a connection error for the real submission, start appending to the file, instead of rewriting it each time. Find a way to keep the column name the same
  gisaid_accessions %>%
    write_csv(here("gisaid", paste0(project_name, "_gisaid_accessions.csv")))
}

#rewrite the fasta seqsender file with the gisaid accession numbers and remove the old prepared files
gisaid_added_seqsender <- fasta_seqsender %>%
  select(-c(`bs-gisaid_accession`)) %>%
  left_join(gisaid_accessions, by = "bs-gisaid_virus_name")

gisaid_added_seqsender %>%
  write_csv(fasta_appended_seqsender_meta_fp)

fastq_seqsender %>%
  select(-c(`bs-gisaid_accession`, `bs-gisaid_virus_name`)) %>%
  left_join(select(gisaid_added_seqsender,
                   all_of(seqsender_merge_check_cols),
                   `bs-gisaid_accession`, `bs-gisaid_virus_name`),
            by = seqsender_merge_check_cols) %>%
  write_csv(fastq_appended_seqsender_meta_fp)

reupload_seqsender_file <- system2("scp", c(fastq_appended_seqsender_meta_fp,
                                            paste0(ec2_hostname, ":", ec2_home_fp, "/", fastq_ec2_seqsender_meta_fp)))

ssh_command_check(reupload_seqsender_file)

### can uncomment the below command to create the NCBI files to check manually before submission
#ssh_seqsender_cmd("prep --biosample --sra", metadata = fastq_ec2_seqsender_meta_fp)

### submit files to biosample and sra
ssh_seqsender_cmd("submit --biosample --sra", metadata = fastq_ec2_seqsender_meta_fp)

#download submission log
download_submission_log <- system2("scp", c(paste0(ec2_hostname, ":", submission_path, "/submission_log.csv"), here()))
ssh_command_check(download_submission_log)

submitted_run <- read_csv(here("submission_log.csv")) %>%
  filter(Submission_Name == input_submission_name,
         if(test_upload) {
           Submission_Type == "TEST"
         } else {
           ubmission_Type == "PRODUCTION"
         })

all_log <- read_csv(file.path(dirname(here()), "submission_log_all.csv")) %>%
  rbind(submitted_run)

all_log %>%
  write_csv(file.path(dirname(here()), "submission_log_all.csv"))

file.remove(here("submission_log.csv"))

if(test_upload) {
  message("You have successfully completed a TEST submission. You may now submit the data for real by setting 'test_upload' to FALSE")
} else {
  message("You have successfully completed a submission to GISAID and NCBI")
}
