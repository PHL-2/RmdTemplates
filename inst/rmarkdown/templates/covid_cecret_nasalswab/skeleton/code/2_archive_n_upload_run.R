library(here)
library(dplyr)

#This Rscript uploads the sequencing run folder and related files to S3 and makes a record of the run on BaseSpace

####################
# Selected variables
####################

# temporary directory to hold the sequencing run download
ec2_tmp_fp <- "~/tmp_bs_dl"

#sequencing date of the run folder should match the RStudio project date
sequencing_date <- gsub("_.*", "", basename(here())) #YYYY-MM-DD

if(sequencing_date == "" | is.na(as.Date(sequencing_date, "%Y-%m-%d")) | nchar(sequencing_date) == 8) {
  stop(simpleError(paste0("Please use the 'YYYY-MM-DD' format for this RStudio project date. This date should correspond to the sequencing run date")))
}

###################################################
# Load functions
###################################################

#this file needs to sit in a [aux_files/r_scripts/functions] directory path above this project directory
tryCatch(
  {
    source(file.path(dirname(here()), "aux_files", "r_scripts", "functions", "R_all_functions_v3.R"))
  },
  error = function(e) {
    stop (simpleError("The R_all_functions_v3.R file needs to sit in a [aux_files/r_scripts/functions] directory path above this project directory"))
  }
)

###################################################
# Load config
###################################################

#this file needs to sit in a [aux_files/r_scripts/config] directory path above this project directory
tryCatch(
  {
    source(file.path(dirname(here()), "aux_files", "r_scripts", "config", "config_variables.R"))
  },
  error = function(e) {
    stop (simpleError("The config_variables.R file needs to sit in a [aux_files/r_scripts/config] directory path above this project directory"))
  }
)

##########################################################
# Set additional variables to check how to run this script
##########################################################

yymmdd <- gsub("^..|-", "", sequencing_date)
seq_folder_pattern <- "[0-9]*_[0-9]*_[0-9A-Z-]*"

selected_sequencer_type <- c("MiSeq", "NextSeq2000")[sequencer_select]

#these Rscripts don't account for two runs that have the same sample types, processed on the same date, on both machines, and the samples need to be processed through the same pipeline
sequencer_regex <- case_when(selected_sequencer_type == "MiSeq" ~ "M",
                             selected_sequencer_type == "NextSeq2000" ~ "VH")

sequencing_folder_regex <- paste0(yymmdd, "_", sequencer_regex, seq_folder_pattern)

s3_run_bucket_fp <- paste0(s3_run_bucket, "/", sequencing_date, "/")

# temporary directory to hold the screen log files
tmp_screen_fp <- paste("~", ".tmp_screen", selected_sequencer_type, paste0(sample_type_acronym, "_", pathogen_acronym), basename(here()), sep = "/")

session_suffix <- tolower(paste(selected_sequencer_type, sample_type_acronym, pathogen_acronym, basename(here()), sep = "-"))

tarball_script <- file.path(dirname(here()), "aux_files", "external_scripts", "bash", "create_tarball.sh")

tarball_script_options <- c("-t", ec2_tmp_fp,
                            "-i", selected_sequencer_type,
                            "-d", sequencing_date,
                            case_when(run_uploaded_2_basespace ~ "-b",
                                      TRUE ~ ""),
                            "-n", case_when(selected_sequencer_type == "MiSeq" ~ miseq_hostname,
                                            selected_sequencer_type == "NextSeq2000" ~ nextseq_hostname,
                                            TRUE ~ ""),
                            "-s", s3_run_bucket,
                            "-r", "Record__",
                            "-o",
                            "-c")

message(paste("Checking if a", selected_sequencer_type, "run on", sequencing_date, "exists on S3"))
check_run_on_s3 <- system2("ssh", c(ec2_hostname,
                                    shQuote(
                                      paste0("aws s3 ls ", s3_run_bucket_fp,
                                                   " | grep ", sequencing_folder_regex, ".tar.gz")
                                      )),
                           stdout = TRUE, stderr = TRUE) %>%
  attr("status") #if run found, returns null. if no run, returns 1

if(is.null(check_run_on_s3)) {
  stop(simpleError(paste(selected_sequencer_type, "run detected in", s3_run_bucket_fp,
                         "\nTo create a new tarball, manually delete the old file in", s3_run_bucket_fp)))
}
message("Run on S3 not found. Continuing to create tarball")

system2("ssh", c(ec2_hostname,
                 shQuote(
                   paste("mkdir -p", tmp_screen_fp)
                 ), "&&",
                 "scp", tarball_script, paste0(ec2_hostname, ":", tmp_screen_fp))
        )

submit_screen_job(message2display = "Creating tarball of the sequencing run folder",
                  ec2_login = ec2_hostname,
                  screen_session_name = paste("sequencing-tarball", session_suffix, sep = "-"),
                  screen_log_fp = tmp_screen_fp,
                  command2run = paste("bash", paste0(tmp_screen_fp, "/", basename(tarball_script)),
                                      paste(tarball_script_options, collapse = " "))
)

check_screen_job(message2display = "Checking tar job",
                 ec2_login = ec2_hostname,
                 screen_session_name = paste("sequencing-tarball", session_suffix, sep = "-"),
                 screen_log_fp = tmp_screen_fp)
