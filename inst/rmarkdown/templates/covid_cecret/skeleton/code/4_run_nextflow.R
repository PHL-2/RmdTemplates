library(here)
library(dplyr)
library(tidyr)
library(stringr)

#This Rscript submits the relevant jobs to Nextflow once the sequencing run has been uploaded

#####################
# Get sequencing date
#####################

system2("aws", c("sso login"))

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

#this file needs to sit in a [aux_files/functions] directory path above this project directory
tryCatch(
  {
    source(file.path(dirname(here()), "aux_files", "functions", "R_all_functions_v3.R"))
  },
  error = function(e) {
    stop (simpleError("The R_all_functions_v3.R file needs to sit in a [aux_files/functions] directory path above this project directory"))
  }
)

#############
# Load config
#############

#this file needs to sit in a [aux_files/config] directory path above this project directory
tryCatch(
  {
    source(file.path(dirname(here()), "aux_files", "config", "config_variables.R"))
  },
  error = function(e) {
    stop (simpleError("The config_variables.R file needs to sit in a [aux_files/config] directory path above this project directory"))
  }
)

#############
# Submit jobs
#############

# Demultiplexing
submit_screen_job(message2display = "Demultiplex with BCLConvert",
                  ec2_login = ec2_hostname,
                  screen_session_name = "demux",
                  command2run = paste("cd ~/.tmp_screen/;",
                                      "nextflow run nf-core/demultiplex",
                                      "-c ~/.nextflow/config",
                                      "-profile", demux_profile,
                                      "-bucket-dir", paste0(s3_nextflow_work_bucket, "/demux_", sequencing_date),
                                      "-resume",
                                      "--input", paste0(s3_run_bucket, "/", sequencing_date, "/", sequencing_date, "_nf_demux_samplesheet.csv"),
                                      "--outdir", paste0(s3_fastq_bucket, "/", sequencing_date, "/processed_bclconvert")))

check_screen_job(message2display = "Checking BCLConvert job",
                 ec2_login = ec2_hostname,
                 screen_session_name = "demux")

# Checking the demultiplexing results
fastq_file_sizes <- system2("aws", c("s3 ls",
                                     paste0(s3_fastq_bucket, "/", sequencing_date, "/processed_bclconvert"),
                                     "--recursive",
                                     "| grep 'R1_001.fastq.gz$'"), stdout = TRUE) %>%
  str_split("\\s+") %>%
  do.call("rbind", .) %>%
  as.data.frame() %>%
  `colnames<-`(c("date", "time", "bytes", "filename")) %>%
  select(bytes, filename) %>%
  mutate(sequencing_folder = gsub(".*processed_bclconvert/", "", filename),
         sequencing_folder = gsub("/.*", "", sequencing_folder),
         filename = gsub(".*/", "", filename),
         bytes = as.numeric(bytes))

undetermined_bytes <- fastq_file_sizes %>%
  filter(grepl("Undetermined", filename)) %>%
  select(bytes) %>%
  pull()

if(undetermined_bytes/sum(fastq_file_sizes$bytes) > 0.5) {
  stop(simpleError("Something might've went wrong with the demultiplexing!\nThe unassigned reads makes up more than 50% of the total reads!"))
}

# Download most recent Nextclade dataset
submit_screen_job(message2display = "Download Nextclade SARS-CoV-2 data",
                  ec2_login = ec2_hostname,
                  screen_session_name = "nextclade-dl",
                  command2run = paste("mkdir -p ~/.local/bin/;",
                                      "wget -q https://github.com/nextstrain/nextclade/releases/latest/download/nextclade-x86_64-unknown-linux-gnu -O ~/.local/bin/nextclade;",
                                      "chmod +x ~/.local/bin/nextclade;",
                                      "nextclade --version;",
                                      "nextclade dataset get --name sars-cov-2 --output-zip ~/sars.zip;",
                                      "aws s3 cp ~/sars.zip", paste0(s3_reference_bucket, "/nextclade/sars.zip;"),
                                      "rm ~/sars.zip"))

check_screen_job(message2display = "Checking Nextclade download job",
                 ec2_login = ec2_hostname,
                 screen_session_name = "nextclade-dl")

# Cecret pipeline
submit_screen_job(message2display = "Process data through Cecret pipeline",
                  ec2_login = ec2_hostname,
                  screen_session_name = "cecret",
                  command2run = paste("cd ~/.tmp_screen/;",
                                      "nextflow run UPHL-BioNGS/Cecret",
                                      "-profile", cecret_profile,
                                      "-bucket-dir", paste0(s3_nextflow_work_bucket, "/cecret_", sequencing_date),
                                      "-r master",
                                      "-resume",
                                      "--reads", paste0(s3_fastq_bucket, "/", sequencing_date, "/processed_bclconvert/", unique(fastq_file_sizes$sequencing_folder)),
                                      "--outdir", paste0(s3_nextflow_output_bucket, "/cecret/", sequencing_date, "_COVIDSeq/processed_cecret")))

check_screen_job(message2display = "Checking Cecret job",
                 ec2_login = ec2_hostname,
                 screen_session_name = "cecret")

# Download BCLConvert files
rstudioapi::executeCommand('activateConsole')
system2("aws", c("s3 cp",
                 paste0(s3_fastq_bucket, "/", sequencing_date),
                 here("data"),
                 "--recursive",
                 "--exclude '*'",
                 "--include '*software_versions.yml'",
                 "--include '*Demultiplex_Stats.csv'",
                 "--include '*Quality_Metrics.csv'",
                 "--include '*Top_Unknown_Barcodes.csv'",
                 "--include '*fastq.gz_fastqc_data.txt'"))

# Download Cecret files
system2("aws", c("s3 cp",
                 paste0(s3_nextflow_output_bucket, "/cecret/", sequencing_date, "_COVIDSeq"),
                 here("data"),
                 "--recursive",
                 "--exclude '*'",
                 "--include '*.consensus.fa'",
                 "--include '*_filtered_R[12].fastq.gz'",
                 "--include '*aggregated-freyja.tsv'",
                 "--include '*_kraken2_report.txt'",
                 "--include '*snp-dists.txt'",
                 "--include '*amplicon_depth.csv'",
                 "--include '*nextclade.tsv'",
                 "--include '*nextclade.json'",
                 "--include '*.stats.txt'",
                 "--include '*.cov.txt'",
                 "--include '*.depth.txt'",
                 "--include '*cecret_results.csv'"))

# Download Nextclade dataset
system2("aws", c("s3 cp",
                 paste0(s3_reference_bucket, "/nextclade/sars.zip"),
                 here("data", "processed_cecret", "nextclade")))

# Download Nextflow config file for profile
run_in_terminal(paste("scp", paste0(ec2_hostname, ":~/.nextflow/config"), here("data", "processed_cecret", "nextflow.config")))
