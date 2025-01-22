library(here)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
system2("aws", c("sso login"))

#This Rscript uploads the sequencing run and related files to S3

####################
# Selected variables
####################

#sequencing date of the run folder should match the RStudio project date
sequencing_date <- gsub("_.*", "", basename(here())) #YYYY-MM-DD

# temporary directory to hold the sequencing run download
ec2_tmp_fp <- "~/tmp_bs_dl"

if(sequencing_date == "") {
  stop (simpleError(paste0("Please fill in the correct sequencing date or short project description in ", here("code"), "/3_archive_upload_run.R")))
} else if (is.na(as.Date(sequencing_date, "%Y-%m-%d")) | nchar(sequencing_date) == 8) {
  stop (simpleError("Please enter the date into [sequencing_date] as YYYY-MM-DD"))
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

############################################
# Tar the sequencing folder and upload to S3
############################################

yymmdd <- gsub("^..|-", "", sequencing_date)
seq_folder_pattern <- "[0-9]*_[0-9]*_[0-9A-Z-]*"

bclconvert_sample_sheet_pattern <- "SampleSheet_v2.csv"

# If this sample sheet is missing, get it from AWS S3 bucket
sample_sheet_fn <- list.files(here("metadata", "munge"), pattern = bclconvert_sample_sheet_pattern)

if(length(sample_sheet_fn) > 1) {
  stop(simpleError("There are more than 2 SampleSheet_v2.csv files detected!! Please delete the incorrect one"))
}

sequencer_type <- gsub("^[0-9-]*_(MiSeq|NextSeq1k2k)_.*", "\\1", sample_sheet_fn)

#these Rscripts don't account for two runs that have the same sample types, processed on the same date, on both machines, and the samples need to be processed through the same pipeline
sequencer_regex <- case_when(sequencer_type == "MiSeq" ~ "M",
                             sequencer_type == "NextSeq1k2k" ~ "VH")

sequencing_folder_regex <- paste0(yymmdd, "_", sequencer_regex, seq_folder_pattern, "$")

s3_run_bucket_fp <- paste0(s3_run_bucket, "/", sequencing_date, "/")

# temporary directory to hold the screen log files
tmp_screen_fp <- paste("~", ".tmp_screen", sequencer_type, "NS_SC2", basename(here()), sep = "/")

session_suffix <- tolower(paste(sequencer_type, "ns-sc2", basename(here()), sep = "-"))

tar_command_function <- function(input_fp, output_fp = input_fp, bars = 99, use_checkpoint = TRUE) {

  # This command runs gzip separately from the tar command, in order to use the '-n' flag
  # The '-n' flag removes the time stamp from the header when compressing, which will allow generating the md5 checksum to be more consistent
  # The flags included in the tar command just tells the program to print a progress bar

  #these are the minimum files from the sequencing run that are required for demultiplexing
  minimum_bs_files_pattern <- c("'SampleSheet.csv$'",
                                "'.xml$'",
                                "'/Data/'",
                                "'/InterOp/'")
  grep_command <- paste("grep", paste("-e", minimum_bs_files_pattern, collapse = " "), " | grep -v '/Analysis/'")

  checkpoint_flags <- if(use_checkpoint) {
    paste("--checkpoint=`find", input_fp, "-type f -exec du -ak --apparent-size {} + |",
          # this 2nd grep command is only used for finding the file sizes to create the progress bar
          grep_command, "| cut -f1 | awk '{s+=$1} END {print int(s /", bars, ")}'`",
          # this flag shows a progress bar
          "--checkpoint-action='ttyout=\b=>'",
          # this flag shows the read/write speed
          #"--checkpoint-action='ttyout=%{Bytes read, Bytes written, Bytes deleted}T, Time elapsed: %ds%*\r'",
          "--record-size=1K")
  } else {
    ""
  }

  paste("cd", paste0(input_fp, ";"),

        ifelse(use_checkpoint,
               paste("echo Estimate:", paste0(c("[", rep("=", bars+1), "];"), collapse = ""),
                     "echo -n 'Progress: [ ';"),
               ""),

        # pipe find command to grep and then tar
        "find . -type f |",
        grep_command,
        "| tar",
        checkpoint_flags,
        "-cT -",
        "| gzip -n >", paste0(output_fp, ".tar.gz;"),

        ifelse(use_checkpoint,
               "echo ']\nTar completed!'",
               "")
  )
}

if(run_uploaded_2_basespace) {

  # Get the run id from BaseSpace
  bs_run <- cli_submit("bs", "list", c("runs", "-f csv")) %>%
    str_split(",") %>%
    do.call("rbind", .) %>%
    as.data.frame() %>%
    `colnames<-`(.[1, ]) %>%
    slice(-1) %>%
    filter(grepl(paste0("^", sequencing_folder_regex), Name))

  if (nrow(bs_run) != 1) {
    stop(simpleError(paste0("\nThere is no sequencing run on BaseSpace for this date: ", sequencing_date,
                            "\nCheck if the date of this Rproject matches with the uploaded sequencing run",
                            "\nThe sequencer type could also be wrong: ", sequencer_type,
                            "\nOtherwise, if you are uploading a local run, set the run_uploaded_2_basespace variable to FALSE")))
  }

  bs_run_id <- bs_run %>%
    select(Id) %>%
    pull()

  sequencing_run <- bs_run %>%
    select(Name) %>%
    pull()

  temporary_seq_run_fp <- paste0(ec2_tmp_fp, "/", sequencing_run)

}

# Run following commands in the cloud if run is on BaseSpace and EC2 access isn't an issue
if(run_uploaded_2_basespace & have_AWS_EC2_SSH_access) {

  submit_screen_job(message2display = "Creating tarball of the sequencing run folder",
                    ec2_login = ec2_hostname,
                    screen_session_name = paste("sequencing-tarball", session_suffix, sep = "-"),
                    screen_log_fp = tmp_screen_fp,
                    command2run = tar_command_function(temporary_seq_run_fp)
  )

  check_screen_job(message2display = "Checking tar job",
                   ec2_login = ec2_hostname,
                   screen_session_name = paste("sequencing-tarball", session_suffix, sep = "-"),
                   screen_log_fp = tmp_screen_fp)

  # Generate md5 checksum
  submit_screen_job(message2display = "Generating md5 checksum",
                    ec2_login = ec2_hostname,
                    screen_session_name = paste("sequencing-checksum", session_suffix, sep = "-"),
                    screen_log_fp = tmp_screen_fp,
                    command2run = paste0("cd ", ec2_tmp_fp, ";",
                                         "md5sum ", sequencing_run, ".tar.gz > ", ec2_tmp_fp, "/", sequencing_run, ".md5;",
                                         "cat ", ec2_tmp_fp, "/", sequencing_run, ".md5")
  )

  check_screen_job(message2display = "Checking md5 job",
                   ec2_login = ec2_hostname,
                   screen_session_name = paste("sequencing-checksum", session_suffix, sep = "-"),
                   screen_log_fp = tmp_screen_fp)

} else {

  dir.create(ec2_tmp_fp, recursive = TRUE)

  # Run the following commands if BaseSpace is available but not SSH (or if the run is MiSeq)
  if(sequencer_type == "MiSeq" | run_uploaded_2_basespace) {

    if(run_uploaded_2_basespace) {

      run_folder <- temporary_seq_run_fp

    } else {

      intended_miseq_folder_regex <- paste0(yymmdd, "_M", seq_folder_pattern)

      scp_command <- paste("scp -r",
                           paste0(miseq_hostname, ":", intended_miseq_folder_regex, "/"),
                           paste0(ec2_tmp_fp, ";"),
                           "sleep 2")

      run_in_terminal(scp_command)

      sequencing_run <- data.frame(run_folders = list.files(ec2_tmp_fp)) %>%
        filter(grepl(paste0(intended_miseq_folder_regex, "$"), run_folders)) %>%
        pull()

      run_folder <- paste0(ec2_tmp_fp, "/", sequencing_run)
    }

    #tar the run folder
    run_in_terminal(paste("echo 'Creating tarball of the sequencing run folder';",
                          tar_command_function(run_folder, use_checkpoint = FALSE)))

  } else if (sequencer_type == "NextSeq1k2k") {

    intended_nextseq_folder_regex <- paste0(yymmdd, "_VH", seq_folder_pattern)
    nextseq_run_fp <- "/usr/local/illumina/runs/"

    sequencing_run <- system2("ssh", c("-tt", nextseq_hostname,
                                       shQuote(paste("cd", paste0(nextseq_run_fp, ";"),
                                                     "ls | grep", intended_nextseq_folder_regex, "| tr -d '\n'"),
                                               type = "sh")),
                              stdout = TRUE)

    run_in_terminal(paste("echo 'Creating NextSeq tar file';",
                          "ssh -tt", nextseq_hostname,
                          shQuote(paste("sleep 5;",
                                        tar_command_function(paste0(nextseq_run_fp, sequencing_run),
                                                             paste0(nextseq_run_fp, sequencing_run, "/", sequencing_run)))))
    )

    scp_command <- paste("scp -r",
                         paste0(nextseq_hostname, ":", nextseq_run_fp, sequencing_run, "/", sequencing_run, ".tar.gz"),
                         paste0(ec2_tmp_fp, ";"),
                         "sleep 2")

    run_in_terminal(scp_command)

    run_folder <- paste0(ec2_tmp_fp, "/", sequencing_run)

  }
  rstudioapi::executeCommand("activateConsole")

  folder_date <- paste0("20", gsub(".*/|_.*", "", run_folder)) %>%
    as.Date(format = "%Y%m%d")

  if(folder_date != gsub("_.*", "", basename(here())) | is.na(folder_date)) {
    stop(simpleError("The run date on the sequencing folder does not match the date of this RStudio project!"))
  }

  message("\nMaking md5 file")
  md5_fp <- paste0(run_folder, ".md5")
  md5_value <- paste0(tools::md5sum(paste0(run_folder, ".tar.gz")), "  ", sequencing_run, ".tar.gz")
  message(md5_value)
  write(md5_value, file = md5_fp)

  message("Uploading local checksum and tarball files to AWS S3")
  s3_cp_md5 <- system2("aws", c("s3 cp", shQuote(path.expand(md5_fp), type = "cmd"), s3_run_bucket_fp), stdout = TRUE)
  s3_cp_run_tarball <- system2("aws", c("s3 cp", shQuote(path.expand(paste0(run_folder, ".tar.gz")), type = "cmd"), s3_run_bucket_fp), stdout = TRUE)

  # If the aws-cli provides an SSL error on local machine, run the command through the instance
  if(length(s3_cp_md5) == 0 | length(s3_cp_run_tarball) == 0) {

    mk_tmp_dir <- system2("ssh", c("-tt", ec2_hostname,
                                   shQuote(paste("mkdir -p", ec2_tmp_fp, "/", session_suffix))),
                          stdout = TRUE, stderr = TRUE)

    if(!grepl("^Connection to .* closed", mk_tmp_dir)) {
      stop(simpleError("Failed to make temporary directory in EC2 instance"))
    }

    run_in_terminal(paste("scp", md5_fp,
                          paste0(ec2_hostname, ":", ec2_tmp_fp, "/", session_suffix))
    )

    run_in_terminal(paste("scp", paste0(run_folder, ".tar.gz"),
                          paste0(ec2_hostname, ":", ec2_tmp_fp, "/", session_suffix))
    )
  }
}
rstudioapi::executeCommand("activateConsole")

sample_sheet_fp <- here("metadata", "munge", sample_sheet_fn)

nf_demux_samplesheet <- data.frame(
  id = sequencing_run,
  samplesheet = paste0(s3_run_bucket_fp, sample_sheet_fn),
  lane = "",
  flowcell = paste0(s3_run_bucket_fp, sequencing_run, ".tar.gz")
)

nfcore_demux_sample_sheet_pattern <- "nf_demux_samplesheet.csv"

nf_demux_samplesheet_fp <- here("metadata", "munge",
                                tolower(paste(sequencing_date, sequencer_type, sample_type_acronym, prj_description, nfcore_demux_sample_sheet_pattern, sep = "_")))

nf_demux_samplesheet %>%
  write_csv(file = nf_demux_samplesheet_fp)

message("Uploading samplesheets to AWS S3")
s3_cp_samplesheet <- system2("aws", c("s3 cp", shQuote(sample_sheet_fp, type = "cmd"), s3_run_bucket_fp), stdout = TRUE)
s3_cp_nf_demux_samplesheet <- system2("aws", c("s3 cp", shQuote(nf_demux_samplesheet_fp, type = "cmd"), s3_run_bucket_fp), stdout = TRUE)

if(length(s3_cp_samplesheet) == 0 | length(s3_cp_nf_demux_samplesheet) == 0) {

  run_in_terminal(paste("scp", sample_sheet_fp,
                        paste0(ec2_hostname, ":", ec2_tmp_fp, "/", session_suffix))
  )

  run_in_terminal(paste("scp", nf_demux_samplesheet_fp,
                        paste0(ec2_hostname, ":", ec2_tmp_fp, "/", session_suffix))
  )
}

if(have_AWS_EC2_SSH_access) {
  # Upload data
  submit_screen_job(message2display = "Uploading run to S3",
                    ec2_login = ec2_hostname,
                    screen_session_name = paste("upload-run", session_suffix, sep = "-"),
                    screen_log_fp = tmp_screen_fp,
                    command2run = paste("aws s3 cp",
                                        paste0(ec2_tmp_fp, "/", session_suffix),
                                        s3_run_bucket_fp,
                                        "--recursive",
                                        "--exclude '*'",
                                        paste0("--include '", sample_sheet_fn, "'"),
                                        paste0("--include '", basename(nf_demux_samplesheet_fp), "'"),
                                        paste0("--include '", sequencing_run, ".md5'"),
                                        paste0("--include '", sequencing_run, ".tar.gz'"))
  )

  check_screen_job(message2display = "Checking run upload job",
                   ec2_login = ec2_hostname,
                   screen_session_name = paste("upload-run", session_suffix, sep = "-"),
                   screen_log_fp = tmp_screen_fp)

  rstudioapi::executeCommand("activateConsole")
}
