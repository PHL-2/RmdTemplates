library(here)
library(dplyr)
library(tidyverse)
library(readxl)
library(readr)
library(stringr)

#This Rscript generates the SampleSheet for demultiplexing a run using BCLConvert and a metadata sheet for analysis of the sequencing run
#https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl_convert/bcl-convert-v3-7-5-software-guide-1000000163594-00.pdf

###################
# Default variables
###################

index_length <- "10"

sequencing_controls <- c("Water control", "Reagent control", "Mock DNA positive control")

sample_group_controls <- c("PBS", "oldWW", "ZeptoSC2")

sample_group_sites <- c("NorthEast", "SouthEast", "SouthWest")

selected_sequencer_type <- c("MiSeq", "NextSeq2000")[sequencer_select]

#sequencing date of the run folder should match the RStudio project date
sequencing_date <- gsub("_.*", "", basename(here())) #YYYY-MM-DD

# temporary directory to hold the screen log files
tmp_screen_fp <- paste("~", ".tmp_screen", selected_sequencer_type, paste0(sample_type_acronym, "_", pathogen_acronym), basename(here()), sep = "/")

session_suffix <- tolower(paste(selected_sequencer_type, sample_type_acronym, pathogen_acronym, basename(here()), sep = "-"))

# temporary directory to hold the sequencing run download
ec2_tmp_fp <- "~/tmp_bs_dl"

#file location of the nextera udi indices
barcode_fp <- file.path(dirname(here()), "aux_files", "illumina_references", "nextera-dna-udi-samplesheet-MiSeq-flex-set-a-d-2x151-384-samples.csv")

if(sequencing_date == "" | is.na(as.Date(sequencing_date, "%Y-%m-%d")) | nchar(sequencing_date) == 8) {
  stop(simpleError(paste0("Please use the 'YYYY-MM-DD' format for this RStudio project date. This date should correspond to the desired sequencing run date")))
}

create_platemap <- FALSE
create_sample_replicates <- 4 #number of biological replicates per sample for sequencing, used to create the platemap

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

########################
# Load barcode sequences
########################

#Sequences for the indices were reformmated from the sample sheet located here (just grabbed the data portion):
#https://support.illumina.com/documents/downloads/productfiles/unique-dual-indexes/nextera-dna-udi/nextera-dna-udi-lrm-samplesheet-iSeq-MiSeq-nextera-dna-flex-set-a-b-c-d-2x151-384-samples.zip
#Which is based on this product page:
#https://support.illumina.com/sequencing/sequencing_kits/idt-nextera-dna-udi/product-files.html
#Illumina Experiment Manager barcode sequences for the NextSeq2000 have different i5 sequences

barcodes <- tryCatch(
  {
    read.csv(barcode_fp, stringsAsFactors = FALSE) %>%
      mutate(idt_plate_coord = paste0(Index_Plate, "_", Index_Plate_Well)) %>%
      mutate(UDI_Index_ID = I7_Index_ID) %>%
      select(idt_plate_coord, I7_Index_ID, I5_Index_ID, UDI_Index_ID, index, index2)
  },
  error = function(e) {
    stop (simpleError("The nextera-dna-udi-samplesheet-MiSeq-flex-set-a-d-2x151-384-samples.csv file needs to sit in an [aux_files/illumina_references] directory path above this project directory"))
  }
)

#####################################
# Load the latest 5 ddPCR run results
#####################################

failed_regex <- "test|exclude"

ddPCR_files <- list.files(ddPCR_run_fp, pattern = ".*_ww_sequencing_metadata.csv", full.names = TRUE, recursive = TRUE)
ddPCR_files <- tail(ddPCR_files[!grepl(failed_regex, ddPCR_files)], 100)

ddPCR_data <- ddPCR_files %>%
  data.frame(FileName = .) %>%
  group_by(FileName) %>%
  do(read_delim(.$FileName,
                col_types = cols("sample_received_date" = col_character(),
                                 "sample_collect_date" = col_character()))) %>%
  ungroup() %>%
  mutate(ddpcr_analysis_date = as.Date(gsub(paste0(ddPCR_run_fp, "/|_.*"), "", FileName)),
         sample_received_date = as.Date(parse_date_time(sample_received_date, c("ymd", "mdy"))),
         sample_collect_date = as.Date(parse_date_time(sample_collect_date, c("ymd", "mdy"))),
         uniq_sample_name = ifelse((is.na(uniq_sample_name) | uniq_sample_name == ""),
                                   paste(sample_type_acronym, sample_received_date, sample_group, sep = "-"),
                                   uniq_sample_name)) %>%
  filter(!is.na(sample_group)) %>%
  group_by(sample_group, sample_received_date) %>%
  #get the latest run only
  filter(ddpcr_analysis_date == max(ddpcr_analysis_date)) %>%
  ungroup() %>%
  select(-FileName) %>%
  unique()

########################################
# Load the environmental samples, if any
########################################

env_fp <- max(list.files(here("metadata", "extra_metadata"), pattern = "environmental_samples.csv", full.names = TRUE))

if(!is.na(env_fp)) {

  env_data <- read_csv(env_fp) %>%
    #use the Tuesday of the sequencing week as the sample_received_date
    mutate(sample_received_date = as.Date(cut(as.POSIXct(sequencing_date), "week")) + 1) %>%
    select(uniq_sample_name = sample_name, sample_group = sample_name,
           sample_received_date, environmental_site) %>%
    #filter rows where sample_id is NA
    filter(!is.na(sample_group)) %>%
    #filter empty columns
    select(where(function(x) any(!is.na(x))),
           !matches("^\\.\\.\\."))
} else{

  env_data <- data.frame(uniq_sample_name = NA_character_,
                         sample_group = NA_character_,
                         sample_received_date = NA_character_,
                         environmental_site = NA_character_)
}

###################################
# Load index and sample info sheets
###################################

index_sheet_fp <- list.files(here("metadata", "munge"), pattern = ".xlsx", full.names = TRUE)

if(identical(index_sheet_fp, character(0))) {
  shared_index_sheet_list <- list.files(file.path(shared_drive_fp, "Sequencing_files", "3_Sample_Sheets", "wastewater", str_sub(sequencing_date, 1, 4)),
                                        pattern = "sequencing_metadata_sheet", full.names = TRUE) %>%
    data.frame(files = .) %>%
    mutate(posted_dates = gsub(".*([0-9-]{8}).*", "\\1", files))

  shared_index_fp <- shared_index_sheet_list %>%
    filter(grepl(format(as.Date(sequencing_date), format = "%m-%d-%y"), files)) %>%
    select(files) %>%
    pull()

  if(identical(shared_index_fp, character(0))) {
    if(nrow(shared_index_sheet_list) == 0) {
      stop(simpleError("\nCannot find files in the shared drive\nAre you connected to the shared drive?"))
    } else {
      posted_dates <- shared_index_sheet_list %>%
        tail(5) %>%
        select(posted_dates) %>%
        pull() %>%
        paste0(collapse = "\n")

      stop(simpleError(paste0("\nCannot find the index sheet with the expected date of ", format(as.Date(sequencing_date), format = "%m-%d-%y"),
                              "\nHere are the dates of the posted index sheets in the shared drive:\n",
                              posted_dates)))
    }
  }

  file.copy(shared_index_fp, here("metadata", "munge"))
  index_sheet_fp <- list.files(here("metadata", "munge"), pattern = ".xlsx", full.names = TRUE)
}

read_sheet <- function(fp, sheet_name) {
  tryCatch(
    {
      read_excel(fp, sheet = sheet_name) %>%
        #filter rows where sample_id is NA
        filter(!is.na(sample_name)) %>%
        #filter empty columns
        select(where(function(x) any(!is.na(x)))) %>%
        select(!matches("^\\.\\.\\.")) %>%
        mutate(across(matches("_col$|coord$"), ~ str_replace_all(., "\\d+", function(m) sprintf("%02d", as.numeric(m))))) %>%
        mutate(plate_coord = gsub("^0", "", plate_coord))
    },
    error = function(e) {
      stop (simpleError("The sample index file from the wet lab scientists may be missing or mis-formatted"))
    }
  )
}

index_sheet <- read_sheet(index_sheet_fp, "Index")
sample_info_sheet <- read_sheet(index_sheet_fp, "Sample Info")

###############
# Get run stats
###############

yymmdd <- gsub("^..|-", "", sequencing_date)
sequencer_regex <- case_when(selected_sequencer_type == "MiSeq" ~ "M",
                             selected_sequencer_type == "NextSeq2000" ~ "VH")
seq_folder_pattern <- "[0-9]*_[0-9]*_[0-9A-Z-]*"
intended_sequencing_folder_regex <- paste0(yymmdd, "_", sequencer_regex, seq_folder_pattern, "$")

record_prefix <- "Record__"

run_cd <- NA
run_q30 <- NA
run_pf <- NA
run_error <- NA

# Get the run id from BaseSpace
bs_run <- system2("ssh", c("-tt", ec2_hostname,
                           shQuote("bs list runs -f csv", type = "sh")),
                  stdout = TRUE, stderr = TRUE) %>%
  head(-1) %>%
  str_split(",") %>%
  do.call("rbind", .) %>%
  as.data.frame() %>%
  `colnames<-`(.[1, ]) %>%
  slice(-1) %>%
  filter(grepl(paste0("^", intended_sequencing_folder_regex), Name)) %>%
  filter(grepl(paste0("^", record_prefix, selected_sequencer_type, "_", sequencing_date), ExperimentName))

if(nrow(bs_run) > 1) {
  stop(simpleError("\nThere are two sequencing runs that matched this date on the same sequencer!!!\n"))
} else if (nrow(bs_run) == 0) {
  stop(simpleError(paste0("\nThere is no record on BaseSpace for this date: ", sequencing_date,
                          "\nCheck if the date of this Rproject matches with the uploaded sequencing run",
                          "\nThe sequencer type could also be wrong: ", selected_sequencer_type)))
}

bs_run_id <- bs_run %>%
  select(Id) %>%
  pull()

sequencing_run <- bs_run %>%
  select(Name) %>%
  pull()

run_stats <- system2("ssh", c("-tt", ec2_hostname,
                              shQuote(paste("bs run seqstats --id", bs_run_id), type = "sh")),
                     stdout = TRUE, stderr = TRUE) %>%
  head(-1) %>%
  list(run_stats = .) %>%
  as.data.frame()

run_cd <- run_stats %>%
  filter(grepl("SequencingStatsCompact.ClusterDensity", run_stats)) %>%
  # cluster density seems to be reported in the millions. If there is no scientific notation, divide by 1000
  mutate(run_stats = gsub(".*\\| | .*", "", run_stats),
         run_stats = ifelse(grepl("e", run_stats),
                            as.numeric(gsub("e.*", "", run_stats))*1000,
                            as.numeric(run_stats)/1000)) %>%
  pull()

run_q30 <- run_stats %>%
  filter(grepl("SequencingStatsCompact.PercentGtQ30 ", run_stats)) %>%
  mutate(run_stats = gsub(".*\\| | .*", "", run_stats),
         run_stats = as.numeric(run_stats)/100) %>%
  pull()

run_pf <- run_stats %>%
  filter(grepl("SequencingStatsCompact.PercentPf", run_stats)) %>%
  mutate(run_stats = gsub(".*\\| | .*", "", run_stats),
         run_stats = as.numeric(run_stats)) %>%
  pull()

run_error <- run_stats %>%
  filter(grepl("SequencingStatsCompact.ErrorRate ", run_stats)) %>%
  mutate(run_stats = gsub(".*\\| | .*", "", run_stats),
         run_stats = as.numeric(run_stats)/100) %>%
  pull()

#######################
# Load run sample sheet
#######################

run_samplesheet_fp <- list.files(here("metadata", "munge"), pattern = "SampleSheet", full.names = TRUE)[1]
samplesheet_exists <- file.exists(run_samplesheet_fp)

if(samplesheet_exists) {

  message("\n*****")
  message("There is already an existing 'SampleSheet' file in the metadata/munge directory")
  message("Using this sheet to generate the metadata...")
  message("*****")
  Sys.sleep(10)

} else {

  dir.create(here("metadata", "munge"), showWarnings = FALSE)
  run_samplesheet_fp <- here("metadata", "munge", "SampleSheet.csv")
  temporary_seq_run_fp <- paste0(ec2_tmp_fp, "/", session_suffix, "/", sequencing_run, "/")
  bs_dl_cmd <- paste("bs download runs --id", bs_run_id, "--output", temporary_seq_run_fp,
                     "--exclude '*' --include 'SampleSheet*'")

  # Download the run from BaseSpace onto a running EC2 instance
  download_bs_sheet_session <- paste0("down-bs-sheet-", session_suffix)
  submit_screen_job(message2display = "Downloading SampleSheet from BaseSpace",
                    ec2_login = ec2_hostname,
                    screen_session_name = download_bs_sheet_session,
                    screen_log_fp = tmp_screen_fp,
                    command2run = bs_dl_cmd
  )

  check_screen_job(message2display = "Checking BaseSpace download job",
                   ec2_login = ec2_hostname,
                   screen_session_name = download_bs_sheet_session,
                   screen_log_fp = tmp_screen_fp)

  # Get name of the final SampleSheet if there is more than 1
  list_sample_sheets <- system2("ssh", c(ec2_hostname,
                                         shQuote(
                                           paste0("ls ", temporary_seq_run_fp, "SampleSheet*")
                                         )),
                                stdout = TRUE, stderr = TRUE) %>%
    tail(1)

  # Download the SampleSheet from EC2 instance
  run_in_terminal(paste("scp",
                        paste0(ec2_hostname, ":", list_sample_sheets),
                        run_samplesheet_fp)
  )
}

run_sample_sheet <- load_sample_sheet(run_samplesheet_fp)

instrument_type <- data.frame(values = unlist(run_sample_sheet$Header)) %>%
  mutate(col_names = gsub(",.*", "", values)) %>%
  mutate(col_names = gsub(" ", "_", col_names)) %>%
  mutate(values = gsub(".*,", "", values)) %>%
  filter(grepl("instrument_type|InstrumentType", col_names, ignore.case = TRUE)) %>%
  select(values) %>%
  pull()

read_length <- data.frame(values = unlist(run_sample_sheet$Reads)) %>%
  filter(!grepl("^Index[1|2]Cycles,", values)) %>%
  mutate(values = gsub(".*,", "", values)) %>%
  pull() %>%
  unique()

if(instrument_type != selected_sequencer_type) {
  message("\n*****")
  message("The SampleSheet.csv for this ", sequencing_date, " run has the instrument set as ", instrument_type)
  message("The rest of this script will continue and processing this project as a ", instrument_type, " run")
  message("If this was not the correct sequencer used for this project, double check the sequencing date or select the appropriate selected_sequencer_type in this Rscript")
  message("*****")

  Sys.sleep(10)
}

instrument_regex <- case_when(instrument_type == "MiSeq" ~ "M",
                              instrument_type == "NextSeq2000" ~ "VH")

if(!read_length %in% c(76, 151)) {
  stop(simpleError("The read length is not 76 or 151 bp. Check the sample sheet from the sequencing run folder"))
}

####################################
# Merge index and sample info sheets
####################################

cols2merge <- c("sample_name", "plate", "plate_row", "plate_col", "plate_coord")

#merge all the sequencing sheets
metadata_sheet <- merge(index_sheet, sample_info_sheet, by = cols2merge, all = TRUE, sort = FALSE) %>%
  #add in barcodes
  merge(barcodes, by = "idt_plate_coord", all.x = TRUE, sort = FALSE) %>%
  arrange(plate, plate_col, plate_row) %>%
  #find sample group from sample_name
  mutate(sample_group = gsub("^WW-([0-9-]+)|^WW-|-Rep.*$", "", sample_name),
         sample_group = case_when(sample_group == "NE" ~ "NorthEast",
                                  sample_group == "SE" ~ "SouthEast",
                                  sample_group == "SW" ~ "SouthWest",
                                  (sample_group == "character(0)" | sample_group == "") ~ NA_character_,
                                  TRUE ~ sample_group),
         sample_id = gsub("_", "-", paste0("PHL2", "-", instrument_regex, "-", idt_plate_coord, "-", gsub("-", "", sequencing_date))),
         uniq_sample_name = gsub("-Rep[0-9]*", "", sample_name),
         sequencing_date = sequencing_date,
         prj_descrip = prj_description,
         instrument_type = instrument_type,
         read_length = read_length,
         index_length = index_length,
         run_cd = run_cd,
         run_q30 = run_q30,
         run_pf = run_pf,
         run_error = run_error,
         sample_received_date = str_extract(sample_name, pattern = "[0-9]{4}-[0,1][0-9]-[0-3][0-9]")) %>%
  select(sample_id, everything())

for(additional_sample_info in c("sample_collection_date", "PHL_sample_received_date")) {
  if(additional_sample_info %in% colnames(metadata_sheet)) {
    stop(simpleError(paste0("Investigate why the ", additional_sample_info, " column was filled out in the metadata sheet by the scientists")))
  }
}

#####################################################################
# Fill in these columns in the metadata sheet if they were left blank
#####################################################################

fill_in_columns <- c("sample_type", "sample_collected_by", "sample_received_date",
                     "organism", "host_scientific_name", "host_disease", "isolation_source",
                     "requester", "requester_email")

#add these columns as NA if missing in metadata sheet
for(x in fill_in_columns) {
  if(!grepl(paste0("^", paste0(colnames(metadata_sheet), collapse = "$|^"), "$"), x)) {
    metadata_sheet[[x]] <- NA
  }
}

#named_vector is the below character vectors containing a grepl pattern as the name and the value as characters
#col_name is the dataframe column to match against the named_vector
multi_grep <- function(named_vector, col_name) {

  ret_vector <- names(named_vector) %>%

    #loop through the names of the named vectors, use that as the grepl pattern
    sapply(., function(x)
      grepl(x, col_name)) %>%

    #if a pattern is found in the pattern, return the column name and use that to select the value
    apply(1, function(y)
      named_vector[names(y)[which(y)]]) %>%

    as.character() %>%

    gsub("character\\(0\\)", NA, .)

  if(length(ret_vector) == 0) {
    stop(simpleError("Your samples' names may not match the usual naming format. Rename the sample or adjust the named sample name vector"))
  }

  ret_vector
}

named_sample_type <- c("^Test-" = "Testing sample type",
                       "^NC-" = "Water control",
                       "^BLANK[0-9]*$|^Blank[0-9]*$" = "Reagent control",
                       "^PC[0-9]*$" = "Mock DNA positive control",
                       "^[A-Z0-9][0-9]*$" = "Nasal swab",
                       "^WW-" = "Wastewater")

extra_cols2merge <- c("uniq_sample_name", "sample_group", "sample_received_date")

metadata_sheet <- metadata_sheet %>%
  mutate(sample_type = case_when(!(is.na(sample_type) | sample_type == "") ~ sample_type,
                                 sample_name %in% env_data$sample_name ~ "Environmental control",
                                 (is.na(sample_type) | sample_type == "") ~ multi_grep(named_sample_type, sample_name),
                                 TRUE ~ NA),
         sample_collected_by = case_when(!(is.na(sample_collected_by) | sample_collected_by == "") ~ sample_collected_by,
                                         TRUE ~ "Philadelphia Water Department"),
         sample_received_date = case_when(!(is.na(sample_received_date) | as.character(sample_received_date) == "") ~ as.character(sample_received_date),
                                          #if sample collect date column is not available, grab the date from the sample_name
                                          grepl("^WW-([0-9-]{10})-", sample_name) ~ gsub("^(WW)-([0-9-]{10})-(.*)", "\\2", sample_name),
                                          #if it's a wastewater sample without a date or does not start with WW, throw an error
                                          sample_type == "Wastewater" ~ NA,
                                          #use Tuesday of the sequencing week if no date specified; older samples that are rerun should have a date manually added in on the sheet
                                          TRUE ~ as.character(as.Date(cut(as.POSIXct(sequencing_date), "week")) + 1)),
         sample_received_date = as.Date(sample_received_date),
         organism = case_when(!(is.na(organism) | organism == "") ~ organism,
                              sample_type == "Wastewater" ~ "Wastewater metagenome",
                              !is.na(sample_type) ~ sample_type,
                              TRUE ~ NA),
         host_scientific_name = case_when(!(is.na(host_scientific_name) | host_scientific_name == "") ~ host_scientific_name,
                                          !is.na(sample_type) ~ "not applicable",
                                          TRUE ~ NA),
         host_disease = case_when(!(is.na(host_disease) | host_disease == "") ~ host_disease,
                                  !is.na(sample_type) ~ "not applicable",
                                  TRUE ~ NA),
         isolation_source = case_when(!(is.na(isolation_source) | isolation_source == "") ~ isolation_source,
                                      sample_type == "Wastewater" ~ "Wastewater",
                                      grepl("control$", sample_type) ~ "Environmental",
                                      grepl("test", sample_type, ignore.case = TRUE) ~ "Test sample",
                                      TRUE ~ NA),
         requester = case_when(!(is.na(requester) | requester == "") ~ requester,
                               TRUE ~ epi_name),
         requester_email = case_when(!(is.na(requester_email) | requester_email == "") ~ requester_email,
                                     TRUE ~ epi_email),
         sample_group = case_when(grepl(paste0(sequencing_controls, collapse = "|"), sample_type) ~ sample_type,
                                  !(is.na(sample_group) | sample_group == "") ~ sample_group,
                                  TRUE ~ gsub(".*-", "", uniq_sample_name)),
         ww_group = case_when(grepl(paste0(sequencing_controls, collapse = "|"), sample_type) ~ sample_type,
                              grepl(paste0(sample_group_controls, collapse = "|"), sample_name) ~ "Wastewater control",
                              TRUE ~ "Wastewater sample")) %>%
  merge(ddPCR_data, by = extra_cols2merge, all.x = TRUE, sort = FALSE) %>%
  merge(env_data, by = extra_cols2merge, all.x = TRUE, sort = FALSE) %>%
  mutate(environmental_site = case_when(grepl(paste0(sequencing_controls, collapse = "|"), sample_type) ~ paste0(sample_name, " - ", plate_row, plate_col),
                                        grepl("Environmental control", sample_type) ~ paste0(environmental_site, " - ", plate_row, plate_col),
                                        TRUE ~ environmental_site)) %>%
  arrange(plate, plate_col, plate_row)

##########################
# Check the metadata sheet
##########################

if(all(is.na(metadata_sheet$ddpcr_analysis_date))) {
  message("\n\nWarning!!!\nSome samples selected for sequencing do not have a corresponding ddPCR date!")
  Sys.sleep(10)
}

main_sample_type <- unique(metadata_sheet$sample_type)[!grepl("control", unique(metadata_sheet$sample_type))]

if(any(is.na(main_sample_type))) {
  message("")
  stop(simpleError(paste0("This metadata sheet has NA in the sample type column!\n",
                          "Probably something went wrong with the merge of the index sheet and the epi's metadata sheet\n",
                          "Here are the samples with NA as its sample type:\n",
                          paste0(metadata_sheet[is.na(metadata_sheet$sample_type), "sample_name"], collapse = ", "))))
}

if(length(main_sample_type) > 1) {
  message("")
  stop(simpleError(paste0("This metadata sheet has more than one non-control sample type!\n",
                          "You may need to separate the metadata sheet and use the appropriate workflow for these samples types:\n",
                          paste0(main_sample_type, collapse = ", "))))
}

if(!grepl("Wastewater|Testing sample type", main_sample_type)) {
  stop(simpleError(paste0("The sample type included in the metadata sheet is not wastewater or a test sample type!\n",
                          "This may not be the appropriate workflow for this run!\n")))
}

sample_type_acronym <- case_when(main_sample_type == "Testing sample type" ~ "Test",
                                 main_sample_type == "Wastewater" ~ "WW")

if(is.na(sample_type_acronym)) {
  stop(simpleError("There's a new or misformatted sample type in the metadata sheet!"))
}

for(x in fill_in_columns) {
  if(any(is.na(metadata_sheet[[x]]))) {
    stop(simpleError(paste0("\n\nWas there a new sample included in this run?\n",
                            "There shouldn't be an NA in column ", x, "\n\n")))
  }
}

#############
# Check sheet
#############

missing_metadata_non_ctrl_samples <- metadata_sheet %>%
  filter(!grepl("control", sample_type))

missing_sample_date <- missing_metadata_non_ctrl_samples %>%
  filter(is.na(sample_received_date)) %>%
  select(sample_name) %>%
  pull() %>%
  str_sort()

#throw error
if(length(missing_sample_date) > 0) {
  message("\n*****")
  message("These non-control samples are in the sample sheet but are missing a collection date!")
  message("Something must be wrong:")
  message("*****")

  stop(simpleError(paste0(missing_sample_date, collapse = ", ")))
}

#check lowest date of sample collection
oldest_ww_date <- metadata_sheet %>%
  filter(!grepl("control", ww_group)) %>%
  select(sample_received_date) %>%
  pull() %>%
  min()

if(oldest_ww_date < seq(as.Date(sequencing_date), length=2, by="-4 month")[2]){
  warning(simpleError(paste0("\nSome samples have collection dates more than 4 months ago!!")))
}

#throw error if missing these columns
for(x in c(cols2merge, "sample_id",
           "idt_set", "idt_plate_row", "idt_plate_col", "idt_plate_coord",
           "index", "index2", "UDI_Index_ID", "I7_Index_ID", "I5_Index_ID",
           "sequencing_date", "prj_descrip", "instrument_type", "read_length", "index_length",
           "environmental_site", "sample_received_date")) {
  if(!grepl(paste0(colnames(metadata_sheet), collapse = "|"), x)) {
    stop(simpleError(paste0("\nMissing column [", x, "] in the metadata sheet template!!!")))
  }
}

#remove columns that are all empty
metadata_sheet <- metadata_sheet %>%
  select(where(function(x) any(!is.na(x))))

message("\nNumber of samples to sequence:")
print(nrow(metadata_sheet))

message("Number of barcodes?")
print(length(unique(paste0(metadata_sheet$index, metadata_sheet$index2))))

message("Dates of samples to sequence:")
print(sort(unique(metadata_sheet$sample_received_date)))

message("Sample sites to sequence:")
print(unique(metadata_sheet$sample_group))

message("What do the sample_id look like?")
print(unique(metadata_sheet$sample_id))

message("Which lanes are sequenced?")
print(unique(metadata_sheet$lane))

message("Are the barcode columns unique?")
if(length(unique(metadata_sheet$idt_plate_coord)) != dim(metadata_sheet)[1]) {
  stop(simpleError("Barcode positions are not unique!"))
}
print(length(unique(metadata_sheet$idt_plate_coord)) == dim(metadata_sheet)[1])

if(length(unique(paste0(metadata_sheet$index, metadata_sheet$index2))) != length(unique(metadata_sheet$sample_id))) {
  stop(simpleError("Differing number of samples and barcodes!"))
}

message("Are the barcodes unique?")
print(length(unique(paste0(metadata_sheet$index, metadata_sheet$index2))) == dim(metadata_sheet)[1])

message("Are the sample names unique?")
print(length(unique(metadata_sheet$sample_id)) == dim(metadata_sheet)[1])

message("Are all the forward primers found?")
print(sum(is.na(metadata_sheet$index)) == 0)

message("Are all the reverse primers found?")
print(sum(is.na(metadata_sheet$index2)) == 0)
if(sum(is.na(c(metadata_sheet$index, metadata_sheet$index2))) != 0) {
  stop(simpleError("Either the forward or reverse primers are NA!"))
}

message("Do all the sampleIDs start with a letter?")
print(all(grepl("^[A-Za-z]", metadata_sheet$sample_id)))
if(!all(grepl("^[A-Za-z]", metadata_sheet$sample_id))) {
  stop(simpleError("Some Sample IDs do not start with a letter!"))
}

message("Are there periods, underscores, or space characters in the SampleID?")
print(any(grepl(" |_|\\.", metadata_sheet$sample_id)))
if(any(grepl(" |_|\\.", metadata_sheet$sample_id))) {
  stop(simpleError("There are spaces, underscores, or periods in the Sample IDs! Please fix"))
}

missing_ddPCR_data <- missing_metadata_non_ctrl_samples %>%
  filter(if_any(ends_with("_gc/L"), is.na)) %>%
  mutate(date_group = paste(sample_received_date, "-", sample_group)) %>%
  select(date_group) %>%
  unique() %>%
  pull() %>%
  str_sort(numeric = TRUE)

#show warning
if(!identical(missing_ddPCR_data, character(0))){
  message("\n\n\n*****")
  message("These non-control samples are in the sample sheet but are missing ddPCR values!!!")
  message("Check to see if the ddPCR has been performed for these samples")
  message(paste("Add the ddPCR data to", ddPCR_run_fp, "and then rerun this script"))
  message("\nSamples in question:")
  print(missing_ddPCR_data)
  message("*****")

  Sys.sleep(15)
}

####################
# Write sample sheet
####################

samp_sheet_2_write <- metadata_sheet %>%
  filter(!sample_name %in% remove_sample_from_bcl_samplesheet) %>%
  filter(!sample_id %in% sample_w_empty_reads) %>%
  # do not include lane in the sample sheet otherwise it will only demultiplex that sample in that specified lane, not in all lanes
  rowwise() %>%
  #BCL Convert does not take Index Plate
  select(sample_id, index, index2) %>%
  rename(Sample_ID = "sample_id")

sample_sheet_fn <- paste(sequencing_date, instrument_type, sample_type_acronym, prj_description, "SampleSheet_v2.csv", sep = "_")

sample_sheet_fp <- here("metadata", "munge", sample_sheet_fn)

write_samp <- function(line2write) {
  write(paste0(line2write, collapse = ","), file = sample_sheet_fp, append = TRUE)
}

write("[Header]", file = sample_sheet_fp)
write_samp(c("FileFormatVersion", "2"))
write_samp(c("RunName", paste0(prj_description, "_", sequencing_date)))
write_samp(c("Date", sequencing_date))
write_samp(c("Workflow", "BCLConvert"))
write_samp(c("InstrumentType", instrument_type))
write_samp("")

write_samp("[Reads]")
write_samp(c("Read1Cycles", read_length))
write_samp(c("Read2Cycles", read_length))
write_samp(c("Index1Cycles", index_length))
write_samp(c("Index2Cycles", index_length))
write_samp("")

write_samp("[Sequencing_Settings]")
write_samp(c("Library Prep Kit", "COVIDSeq for Surveillance"))
write_samp(c("Index Kit", "COVIDSeq indexes_IDT for Illumina-PCR Indexes Set 1 2 3 4"))
write_samp(c("Chemistry", "Amplicon"))
write_samp("")

write_samp("[BCLConvert_Settings]")
write_samp(c("CreateFastqForIndexReads", "1"))
#don't split the lanes
write_samp(c("NoLaneSplitting", "true"))
write_samp(c("FastqCompressionFormat", "gzip"))
write_samp("")

write_samp("[BCLConvert_Data]")
write_csv(samp_sheet_2_write, file = sample_sheet_fp, col_names = TRUE, append = TRUE)

##########################################################
# Generate a separate sample sheet for nf-core/demultiplex
##########################################################

s3_run_bucket_fp <- paste0(s3_run_bucket, "/", sequencing_date, "/")

nf_demux_samplesheet <- data.frame(
  id = sequencing_run,
  samplesheet = paste0(s3_run_bucket_fp, sample_sheet_fn),
  lane = "",
  flowcell = paste0(s3_run_bucket_fp, sequencing_run, ".tar.gz")
)

nfcore_demux_sample_sheet_pattern <- "nf_demux_samplesheet.csv"

nf_demux_samplesheet_fp <- here("metadata", "munge",
                                tolower(paste(sequencing_date, instrument_type, sample_type_acronym, prj_description, nfcore_demux_sample_sheet_pattern, sep = "_")))

nf_demux_samplesheet %>%
  write_csv(file = nf_demux_samplesheet_fp)

###################################
# Upload sample sheets to S3 bucket
###################################

ec2_tmp_session_dir <- paste0(ec2_tmp_fp, "/", session_suffix, "/")

mk_remote_dir(ec2_hostname, ec2_tmp_session_dir)

run_in_terminal(paste("scp", sample_sheet_fp, nf_demux_samplesheet_fp,
                      paste0(ec2_hostname, ":", ec2_tmp_session_dir))
)

upload_samplesheet_session <- paste0("up-samplesheet-", session_suffix)
submit_screen_job(message2display = "Uploading sample sheets to S3",
                  ec2_login = ec2_hostname,
                  screen_session_name = upload_samplesheet_session,
                  screen_log_fp = tmp_screen_fp,
                  command2run = paste("aws s3 cp",
                                      ec2_tmp_session_dir,
                                      s3_run_bucket_fp,
                                      "--recursive",
                                      "--exclude '*'",
                                      paste0("--include '", sample_sheet_fn, "'"),
                                      paste0("--include '", basename(nf_demux_samplesheet_fp), "'"))
)

check_screen_job(message2display = "Checking sample sheet upload job",
                 ec2_login = ec2_hostname,
                 screen_session_name = upload_samplesheet_session,
                 screen_log_fp = tmp_screen_fp)

################################
# Write sheet to metadata folder
################################

merged_samples_metadata_sheet <- metadata_sheet %>%
  select(-c(sample_id, index, index2,
            starts_with("idt_"), starts_with("plate"), ends_with("_ID", ignore.case = FALSE))) %>%
  filter(!sample_type %in% c(sequencing_controls, "Environmental control")) %>%
  group_by(uniq_sample_name) %>%
  mutate(sample_counts = n()) %>%
  ungroup() %>%
  filter(sample_counts > 1) %>%
  select(-c(sample_name, sample_counts)) %>%
  unique() %>%
  mutate(index = "TTTTTTTTTT",
         index2 = "TTTTTTTTTT",
         UDI_Index_ID = "UDP9999",
         idt_set = "Merged",
         idt_plate_row = rep(LETTERS[1:8], 12)[row_number()],
         idt_plate_col = unlist(lapply(1:12, function (x) rep(x, 8)))[row_number()],
         idt_plate_coord = paste0(idt_set, "_", idt_plate_row, idt_plate_col),
         plate = 999,
         plate_row = idt_plate_row,
         plate_col = idt_plate_col,
         plate_coord = paste0(plate, "_", plate_row, plate_col),
         across(matches("_col$|coord$"), ~ str_replace_all(., "\\d+", function(m) sprintf("%02d", as.numeric(m)))),
         sample_id = gsub("_", "-", paste0("PHL2", "-", instrument_regex, "-", idt_plate_coord, "-", gsub("-", "", sequencing_date))),
         sample_name = paste0(uniq_sample_name, "-Merged"),
         ww_group = paste(ww_group, "- Merged"))

#does not contain PHI and accession numbers
metadata_sheet %>%
  select(-c(I7_Index_ID, I5_Index_ID)) %>%
  # NextSeq runs have index2 sequences in the reverse complement of the ones listed in the reference barcode sheet
  # however, the sample sheet used for demultiplexing needs to be the same orientation as the reference barcode sheet
  # when BCLConvert demultiplexes a NextSeq run, it will automatically reverse complement index2
  # results from the pipeline will refer to the reverse complement of index2, so this should be updated in the metadata sheet
  mutate(index2 = ifelse(instrument_type == "NextSeq2000", reverse_complement(index2), index2)) %>%
  rbind(merged_samples_metadata_sheet) %>%
  write_csv(file = here("metadata", paste0(sequencing_date, "_", prj_description, "_metadata.csv")))

bclconvert_output_final_path <- paste(s3_fastq_bucket, sequencing_date, sample_type_acronym, prj_description, "processed_bclconvert", sequencing_run, sep = "/")

#write the sample sheet for merging nextflow script
metadata_sheet %>%
  filter(!sample_name %in% remove_sample_from_bcl_samplesheet) %>%
  filter(!sample_id %in% sample_w_empty_reads) %>%
  select(fastq = sample_id, uniq_sample_name) %>%
  mutate(fastq_1 = paste0(bclconvert_output_final_path, "/", fastq, "_S", row_number(), "_R1_001.fastq.gz"),
         fastq_2 = paste0(bclconvert_output_final_path, "/", fastq, "_S", row_number(), "_R2_001.fastq.gz")) %>%
  merge(select(merged_samples_metadata_sheet, sample_id, uniq_sample_name), by = "uniq_sample_name", all = TRUE) %>%
  filter(!is.na(sample_id)) %>%
  select(sample_id, fastq_1, fastq_2) %>%
  write_csv(file = here("metadata", "munge",
                        tolower(paste(sequencing_date, instrument_type, sample_type_acronym, prj_description, "nf_concat_fastq_samplesheet.csv", sep = "_"))))

######################################################
# Create a platemap from the metadata sheet, if needed
######################################################

if(create_platemap) {
  dir.create(here("metadata", "for_scientists"))

  empty_plate <- data.frame(plate_row = unlist(lapply(LETTERS[1:8], function(x) rep(x, 12))), plate_col = sprintf("%02d", rep(1:12, 8)), plate = 1) %>%
    mutate(plate_coord = paste0(plate, "_", plate_row, plate_col),
           half_plate = !plate_row %in% LETTERS[1:4]) %>%
    arrange(half_plate, plate_col) %>%
    mutate(sample_order = row_number()) %>%
    select(plate, plate_row, plate_col, plate_coord, sample_order)

  sample_group_order <- c("ZeptoSC2", "SouthWest", "NorthEast", "SouthEast", "oldWW")

  grouped_samples <- select(metadata_sheet, sample_group, sample_received_date, uniq_sample_name) %>%
    unique() %>%
    filter(sample_group != "",
           !is.na(sample_received_date)) %>%
    mutate(sample_group = factor(sample_group, levels = sample_group_order)) %>%
    arrange(sample_received_date, sample_group) %>%
    mutate(order = 1:nrow(.)) %>%
    expand(nesting(uniq_sample_name, order), rep = paste0("-Rep", 1:create_sample_replicates)) %>%
    mutate(sample_name = paste0(uniq_sample_name, rep)) %>%
    arrange(order) %>%
    select(sample_name) %>%
    #put samples in groups of 4
    mutate(grp = (row_number() - 1) %/% 4)

  # Add in controls to the plate
  if(max(grouped_samples$grp) > 4) {
    combined_list <- grouped_samples %>%
      add_row(., sample_name = paste0("NC-pre-extract", 5:8), .before = 61) %>%
      add_row(., sample_name = c("NC-pre-extract9", "BLANK", "PC", "NC-pre-cDNA"), .before = 41) %>%
      add_row(., sample_name = paste0("NC-pre-extract", 1:4), .before = 21)
  } else {
    combined_list <- grouped_samples %>%
      add_row(., sample_name = c("NC-pre-extract5", "BLANK", "PC", "NC-pre-cDNA"), .before = 21) %>%
      add_row(., sample_name = paste0("NC-pre-extract", 1:4), .before = ceiling(median(.$grp, na.rm = TRUE))*4+1)
  }

  plate_view <- combined_list %>%
    filter(sample_name != "") %>%
    select(-grp) %>%
    rbind(data.frame(sample_name = c("NC-pre-ARTIC", "NC-pre-library"))) %>%
    mutate(sample_order = row_number()) %>%
    merge(empty_plate, by = "sample_order", all = TRUE, sort = FALSE) %>%
    mutate(sample_name = case_when(sample_order == 96 ~ "NC-corner",
                                   is.na(sample_name) ~ "",
                                   TRUE ~ sample_name)) %>%
    arrange(sample_order)

  real_plate_view <- plate_view %>%
    select(sample_name, plate_row, plate_col) %>%
    pivot_wider(names_from = "plate_col", values_from = "sample_name")

  write_csv(plate_view, file = here("metadata", "for_scientists", paste0(format(Sys.time(), "%Y%m%d"), "_combined_samples_list.csv")))

  plate_map_local_fp <- here("metadata", "for_scientists", paste0(format(Sys.time(), "%Y%m%d"), "_combined_samples_plate_map.csv"))
  write_csv(real_plate_view, file = plate_map_local_fp)

  # if the shared drive can be found, copy the platemap over
  if(file.exists(shared_drive_fp)) {
    plate_map_cp_fp <- file.path(shared_drive_fp, "Sequencing_files", "1_Plate_Maps", "wastewater", format(Sys.Date(), "%Y"),
                                 paste(format(Sys.time(), "%Y-%m-%d"), sample_type_acronym, prj_description, "Plate_Map.csv", sep = "_"))

    file.copy(plate_map_local_fp, plate_map_cp_fp, overwrite = TRUE)

  } else{
    message("\n*****")
    message("Could not access shared drive path. Plate map not copied")
    message("*****")
    Sys.sleep(5)
  }
}

message("\nRscript finished successfully!")
