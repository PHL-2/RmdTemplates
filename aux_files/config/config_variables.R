#s3 bucket holding the sequencing runs
s3_run_bucket <- "s3://"

#location on local drive of where sequencing run folder is held
sequencing_folder_fp <- ""

#shared drive location
shared_drive_fp <- ""

#shell command
sh_exe_fp <- "sh"

bioinformatician_name <- "First Last, Degree"

#path changes on Windows system
if(Sys.info()["sysname"] == "Windows")
{
  #shared drive location
  shared_drive_fp <- ""

  #location on local drive of where sequencing run folder is held
  sequencing_folder_fp <- ""

  #full file path of the git shell executable, dependent on where git was installed
  #on a Windows OS, it may be installed at C:/Users/user.name/AppData/Local/Programs/Git/bin/sh.exe
  sh_exe_fp <- "path/to/Git/bin/sh.exe"

  #full path of where the seqsender repo was cloned
  seqsender_fp <- "path/to/seqsender"
}
