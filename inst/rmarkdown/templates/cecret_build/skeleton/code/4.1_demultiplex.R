#############
# Submit jobs
#############

# If this sample sheet is missing, get it from AWS S3 bucket
sample_sheet_fn <- list.files(here("metadata", "munge"), pattern = "SampleSheet_v2.csv")

if(length(sample_sheet_fn) > 1) {
  stop(simpleError("There are more than 2 sample sheets detected!! Please delete the incorrect one"))
}

sequencer_type <- gsub("^[0-9-]*_(MiSeq|NextSeq1k2k)_.*", "\\1", sample_sheet_fn)

sequencer_regex <- case_when(sequencer_type == "MiSeq" ~ "M",
                             sequencer_type == "NextSeq1k2k" ~ "VH")

intended_sequencing_folder_regex <- paste0(gsub("^..|-", "", sequencing_date), "_", sequencer_regex, "[0-9]*_[0-9]*_[0-9A-Z-]*")

sample_type_acronym <- gsub(paste0("^[0-9-]*_", sequencer_type, "_|_.*"), "", sample_sheet_fn)

prj_description <- gsub(paste0("^[0-9-]*_.*", sample_type_acronym, "_|_.*"), "", sample_sheet_fn)

nf_demux_samplesheet_path <- paste(s3_run_bucket, sequencing_date,
                                   tolower(paste(sequencing_date, sequencer_type, sample_type_acronym, prj_description, "nf_demux_samplesheet.csv", sep = "_")), sep = "/")

bclconvert_output_path <- paste(s3_fastq_bucket, sequencing_date, sample_type_acronym, prj_description, "processed_bclconvert", sep = "/")

workflow_output_fp <- paste(s3_nextflow_output_bucket, "cecret", sample_type_acronym, paste0(sequencing_date, "_", prj_description), sequencer_type, sep = "/")

data_output_fp <- paste0(ec2_tmp_fp, "/", sequencing_date, "/data")

# Demultiplexing
submit_screen_job(message2display = "Demultiplexing with BCLConvert",
                  ec2_login = ec2_hostname,
                  screen_session_name = "demux",
                  command2run = paste("cd ~/.tmp_screen/;",
                                      "nextflow run nf-core/demultiplex",
                                      "-c ~/.nextflow/config",
                                      "-profile", demux_profile,
                                      "-bucket-dir", paste0(s3_nextflow_work_bucket, "/demux_", sample_type_acronym, "_", sequencing_date),
                                      "-resume",
                                      "--input", nf_demux_samplesheet_path,
                                      "--outdir", bclconvert_output_path)
)

check_screen_job(message2display = "Checking BCLConvert job",
                 ec2_login = ec2_hostname,
                 screen_session_name = "demux")

# Checking the demultiplexing results
aws_s3_fastq_files <- system2("aws", c("s3 ls", bclconvert_output_path,
                                       "--recursive",
                                       "| grep 'R1_001.fastq.gz$'",
                                       "| grep -v 'Merged'"), stdout = TRUE)
