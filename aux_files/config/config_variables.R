bioinformatician_name <- "First Last, Degree"

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
s3_nf_work_bucket <- "s3://"

#host name of the ec2 instance in the .ssh/config file
ec2_hostname <- ""

###################
# Nextflow settings
###################

# Nextflow profiles to use for the demultiplexing and Cecret pipelines (should be defined in the .nextflow/config file)
demux_profile <- ""
cecret_profile <- ""

################
# Local settings
################

# full path on local drive of where sequencing run folder is held
# this variable only needs to be defined if you're doing a manual upload
sequencing_folder_fp <- "path/to/local/sequencing/folder"

# full path of the shared drive to upload the pdf report
shared_drive_fp <- "path/to/shared/drive"

# full path of where the seqsender repo was cloned
seqsender_fp <- "path/to/seqsender"

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
