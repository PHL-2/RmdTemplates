library(here)
library(dplyr)
library(tidyr)
library(stringr)

#This Rscript uploads the sequencing run and related files to S3

##############
# Manual input
##############

run_uploaded_2_basespace <- TRUE #is the sequencing run on basespace?

# temporary directory in ec2 to hold to sequencing run download. This directory will be deleted after running this script
ec2_tmp_fp <- "~/tmp_bs_dl/"

#sequencing date of the run folder should match the RStudio project date
sequencing_date <- gsub("_.*", "", basename(here())) #YYYY-MM-DD

if(sequencing_date == "") {
  stop (simpleError(paste0("Please fill in the correct sequencing date or short project description in ", here("code"), "/3_archive_upload_run.R")))
} else if (is.na(as.Date(sequencing_date, "%Y-%m-%d")) | nchar(sequencing_date) == 8) {
  stop (simpleError("Please enter the date into [sequencing_date] as YYYY-MM-DD"))
}

###################################################
# Load functions
###################################################

#this file needs to sit in a [aux_files/functions] directory path above this project directory
tryCatch(
  {
    source(file.path(dirname(here()), "aux_files", "functions", "R_all_functions_v3.R"))
  },
  error = function(e) {
    stop (simpleError("The R_all_functions_v3.R file needs to sit in a [aux_files/functions] directory path above this project directory"))
  }
)

###################################################
# Load config
###################################################

#this file needs to sit in a [aux_files/config] directory path above this project directory
tryCatch(
  {
    source(file.path(dirname(here()), "aux_files", "config", "config_variables.R"))
  },
  error = function(e) {
    stop (simpleError("The config_variables.R file needs to sit in a [aux_files/config] directory path above this project directory"))
  }
)

############################################
# Tar the sequencing folder and upload to S3
############################################

# If this sample sheet is missing, get it from AWS S3 bucket
sample_sheet_fn <- list.files(here("metadata", "munge"), pattern = "SampleSheet_v2.csv")

if(length(sample_sheet_fn) > 1) {
  stop(simpleError("There are more than 2 sample sheets detected!! Please delete the incorrect one"))
}

sequencer_type <- gsub("^[0-9-]*_(MiSeq|NextSeq1k2k)_.*", "\\1", sample_sheet_fn)

sample_type_acronym <- gsub(paste0("^[0-9-]*_", sequencer_type, "_|_.*"), "", sample_sheet_fn)

prj_description <- gsub(paste0("^[0-9-]*_.*", sample_type_acronym, "_|_.*"), "", sample_sheet_fn)

s3_run_bucket_fp <- paste0(s3_run_bucket, "/", sequencing_date, "/")

system2("aws", c("sso login"))

sequencing_folder_regex <- paste0(gsub("^..|-", "", sequencing_date), "_([M]{1}|[VH]{2})[0-9]*_[0-9]*_[0-9A-Z-]*$")

if(run_uploaded_2_basespace) {

  # Get the run id from BaseSpace
  bs_run <- cli_submit("bs", "list", c("runs", "-f csv")) %>%
    str_split(",") %>%
    do.call("rbind", .) %>%
    as.data.frame() %>%
    `colnames<-`(.[1, ]) %>%
    slice(-1) %>%
    filter(grepl(paste0("^", sequencing_folder_regex), Name))

  if(nrow(bs_run) > 1) {
    warning(simpleWarning(paste0("\nThere are two sequencing runs that matched this date. Make sure you selected the correct sequencer!!!\n",
                                 "Currently, you are pulling the sequencing run from the ", sequencer_type, "\n\n")))

    #these Rscripts don't account for two runs that have the same sample types, processed on the same date, on both machines, and the samples need to be processed through the same pipeline
    sequencer_regex <- case_when(sequencer_type == "MiSeq" ~ "M",
                                 sequencer_type == "NextSeq1k2k" ~ "VH")

    intended_sequencing_folder_regex <- paste0(gsub("^..|-", "", sequencing_date), "_", sequencer_regex, "[0-9]*_[0-9]*_[0-9A-Z-]*$")

    bs_run <- bs_run %>%
      filter(grepl(paste0("^", intended_sequencing_folder_regex), Name))

  }
  if (nrow(bs_run) == 0) {
    stop(simpleError(paste0("\nThere is no sequencing run on BaseSpace matching this pattern: ", intended_sequencing_folder_regex,
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

  sequencing_run_fp <- paste0(ec2_tmp_fp, sequencing_run, "/")

  # This command runs gzip separately from the tar command, in order to use the '-n' flag
  # The '-n' flag removes the time stamp from the header when compressing, which will allow generating the md5 checksum to be more consistent
  # The flags included in the tar command just tells the program to print a progress bar
  bars <- 99
  submit_screen_job(message2display = "Creating tarball of the sequencing run folder",
                    ec2_login = ec2_hostname,
                    screen_session_name = "sequencing-tarball",
                    command2run = paste("echo Estimated:", paste0(c("[", rep("=", bars+1), "];"), collapse = ""),
                                        "echo -n 'Progress:  [ ';",
                                        "tar",
                                        "--checkpoint=`du -sk --apparent-size", sequencing_run_fp, "| cut -f1 | awk '{print int($1 /", bars, ")}'`",
                                        # this flag shows a progress bar
                                        "--checkpoint-action='ttyout=\b=>'",
                                        # this flag shows the read/write speed
                                        #"--checkpoint-action='ttyout=%{Bytes read, Bytes written, Bytes deleted}T, Time elapsed: %ds%*\r'",
                                        "--record-size=1K",
                                        "-cC", sequencing_run_fp, ". |",
                                        "gzip -n >", paste0(ec2_tmp_fp, sequencing_run, ".tar.gz;"),
                                        "echo ']\nTar completed!'")
                    )

  check_screen_job(message2display = "Checking tar job",
                   ec2_login = ec2_hostname,
                   screen_session_name = "sequencing-tarball")

  # Generate md5 checksum
  submit_screen_job(message2display = "Generating md5 checksum",
                    ec2_login = ec2_hostname,
                    screen_session_name = "sequencing-checksum",
                    command2run = paste0("cd ", ec2_tmp_fp, ";",
                                         "md5sum ", sequencing_run, ".tar.gz > ", ec2_tmp_fp, sequencing_run, ".md5;",
                                         "cat ", ec2_tmp_fp, sequencing_run, ".md5")
                    )

  check_screen_job(message2display = "Checking md5 job",
                   ec2_login = ec2_hostname,
                   screen_session_name = "sequencing-checksum")

  # Upload data

  submit_screen_job(message2display = "Uploading run to S3",
                    ec2_login = ec2_hostname,
                    screen_session_name = "upload-run",
                    command2run = paste("aws s3 cp",
                                        ec2_tmp_fp,
                                        s3_run_bucket_fp,
                                        "--recursive",
                                        "--exclude '*'",
                                        paste0("--include '", sequencing_run, ".md5'"),
                                        paste0("--include '", sequencing_run, ".tar.gz'"))
                    )

  check_screen_job(message2display = "Checking run upload job",
                   ec2_login = ec2_hostname,
                   screen_session_name = "upload-run")

  rstudioapi::executeCommand('activateConsole')

} else {

  # If the run was not uploaded to BaseSpace, run this part
  # Get the path of the local run folder
  run_folder <- sequencing_folder_fp %>%
    list.files(full.names = T) %>%
    data.frame(filenames = .) %>%
    filter(grepl(format(as.Date(sequencing_date), "%y%m%d"), filenames)) %>%
    filter(grepl(sequencing_folder_regex, filenames)) %>%
    filter(!grepl("\\.tar\\.gz$|\\.md5$", filenames)) %>%
    pull()

  folder_date <- paste0("20", gsub(".*/|_.*", "", run_folder)) %>%
    as.Date(format = "%Y%m%d")

  sequencing_run <- gsub(".*/", "", run_folder)

  if(folder_date != gsub("_.*", "", basename(here())) | is.na(folder_date)) {
    stop(simpleError("The run date on the sequencing folder does not match the date of this RStudio project!"))
  }

  #tar the run folder
  message("Making sequencing run tarball")
  system2("tar", c("-czf", shQuote(paste0(run_folder, ".tar.gz"), type = "cmd"),
                   "-C", shQuote(run_folder, type = "cmd"), "."))

  message("Making md5 file")
  md5_fp <- file.path(sequencing_folder_fp, paste0(sequencing_run, ".md5"))
  paste0(sequencing_run, ".tar.gz") %>%
    paste0(tools::md5sum(file.path(sequencing_folder_fp, .)), "  ", .) %>%
    write(file = md5_fp)

  message("Uploading local checksum and tarball files to AWS S3")
  s3_cp_md5 <- system2("aws", c("s3 cp", shQuote(md5_fp, type = "cmd"), s3_run_bucket_fp), stdout = TRUE)
  s3_cp_run_tarball <- system2("aws", c("s3 cp", shQuote(paste0(run_folder, ".tar.gz"), type = "cmd"), s3_run_bucket_fp), stdout = TRUE)

  if(!all(grepl("Completed", c(s3_cp_md5, s3_cp_run_tarball), ignore.case = TRUE))) {
    stop(simpleError("Local upload of checksum and tarball files to s3 bucket failed!"))
  }

}

sample_sheet_fp <- here("metadata", "munge", sample_sheet_fn)

nf_demux_samplesheet <- data.frame(
  id = sequencing_run,
  samplesheet = paste0(s3_run_bucket_fp, sample_sheet_fn),
  lane = "",
  flowcell = paste0(s3_run_bucket_fp, sequencing_run, ".tar.gz")
)

nf_demux_samplesheet_fp <- here("metadata", "munge",
                                tolower(paste(sequencing_date, sequencer_type, sample_type_acronym, prj_description, "nf_demux_samplesheet.csv", sep = "_")))

nf_demux_samplesheet %>%
  write.csv(file = nf_demux_samplesheet_fp,
            row.names = FALSE, quote = FALSE)

message("Uploading samplesheets to AWS S3")
s3_cp_samplesheet <- system2("aws", c("s3 cp", shQuote(sample_sheet_fp, type = "cmd"), s3_run_bucket_fp), stdout = TRUE)
s3_cp_nf_demux_samplesheet <- system2("aws", c("s3 cp", shQuote(nf_demux_samplesheet_fp, type = "cmd"), s3_run_bucket_fp), stdout = TRUE)

if(!all(grepl("Completed", c(s3_cp_samplesheet, s3_cp_nf_demux_samplesheet), ignore.case = TRUE)) |
   length(s3_cp_samplesheet) == 0 |
   length(s3_cp_nf_demux_samplesheet) == 0) {
  stop(simpleError("Local upload of sample sheets to s3 bucket failed"))
}

rstudioapi::executeCommand('activateConsole')
