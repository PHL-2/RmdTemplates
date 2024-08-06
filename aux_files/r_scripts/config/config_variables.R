#############
# People info
#############
bioinformatician_name <- "First Last, Degree"
epi_name <- NA
epi_email <- NA

#################
# AWS S3 settings
#################

# s3 bucket holding the raw sequencing data (including the md5sum file and the sample sheet)
# this variable is mandatory if using the COVID_Cecret RmdTemplate
s3_run_bucket <- "s3://"
# s3 bucket holding the demultiplexed fastq files for each sample
s3_fastq_bucket <- "s3://"
# s3 bucket holding the Nextflow output files
s3_nextflow_output_bucket <- "s3://"
# s3 bucket to use as temporary storage for Nextflow intermediate files
s3_nextflow_work_bucket <- "s3://"
# s3 bucket holding reference data files
s3_reference_bucket <- "s3://"

##################
# AWS EC2 settings
##################

#host name of the ec2 instance in the .ssh/config file
ec2_hostname <- ""

# full path of where the seqsender repo was cloned on EC2 instance that holds the seqsender.py file
seqsender_fp <- "path/to/seqsender"

# full path of where to keep the fasta and fastq files for upload
seqsender_upload_tmp_fp <- "path/to/upload/folder"

####################
# Sequencer settings
####################

#host name of the NextSeq sequencer in the .ssh/config file
nextseq_hostname <- ""

###################
# Nextflow settings
###################

# Nextflow profiles to use for the demultiplexing and Cecret pipelines (should be defined in the .nextflow/config file)
demux_profile <- ""
cecret_profile <- ""

####################
# GISAID credentials
####################

gisaid_client_id <- ""
gisaid_username <- ""

################
# Local settings
################

# full path on local drive of where sequencing run folder is held
# this variable only needs to be defined if you're doing a manual upload
sequencing_folder_fp <- "path/to/local/sequencing/folder"

# full path of the shared drive to upload the pdf report
shared_drive_fp <- "path/to/shared/drive"

# full path of the directory holding the WW ddPCR run info
ddPCR_run_fp <- "path/to/local/ddpcr/folder"

# MacOS shell command
sh_exe_fp <- "sh"

# path changes on Windows system
if(Sys.info()["sysname"] == "Windows")
{

  # full path on local drive of where sequencing run folder is held
  # this variable only needs to be defined if you're doing a manual upload
  sequencing_folder_fp <- "Windows/path/to/local/sequencing/folder"

  # full path of the shared drive to upload the pdf report
  shared_drive_fp <- "Windows/path/to/shared/drive"

  # full path of where the seqsender repo was cloned
  seqsender_fp <- "Windows/path/to/seqsender"

  # full path of the git shell executable, depending on where git was installed
  # on Windows OS, it may be installed at C:/Users/user.name/AppData/Local/Programs/Git/bin/sh.exe
  sh_exe_fp <- "path/to/Git/bin/sh.exe"

}

#find s3 path variables and remove any trailing slashes
ls(pattern = "^s3_|_fp$") %>%
  mapply(assign,
         .,
         sapply(., function(x) gsub("/$", "", eval(parse(text = x)))),
         MoreArgs=list(envir=parent.frame()))
