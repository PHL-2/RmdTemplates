# YYYY-MM-DD
date_pipeline_was_run <- ""

# N.N.N
dragen_covid_lineage_version <- ""

# N.N.N
nextclade_version <- ""

if(date_pipeline_was_run == "" | dragen_covid_lineage_version == "" | nextclade_version == "") {
  stop (simpleError("Fill out the COVIDSeq_config.R file in the sequencing run folder with the relevant information first"))
}

primer_select <-  # set variable as 1 for V3, 2 for V4, or 3 for V4.1
  
artic_primer_scheme <- c("V3", "V4", "V4.1")[primer_select]
bs_bed_file_used <- c("nCov-2019.bed", "artic_V4_primer_scheme_NC_045512.2.bed", "artic_v4.1_primers_nc045512.2.bed")[primer_select]

# minimum acceptable number of reads
min_reads <- 10000