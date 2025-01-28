library(here)
library(dplyr)
library(tidyr)
library(readxl)
library(readr)
library(stringr)

#This Rscript currently generates the SampleSheet for demultiplexing a run using the BCL Convert program
#https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl_convert/bcl-convert-v3-7-5-software-guide-1000000163594-00.pdf

###################
# Default variables
###################

index_length <- "10"

phi_info <- c("sample_name", "zip_char", "case_id", "breakthrough_case", "death", "hospitalized", "outbreak", "priority")

sequencer_type <- c("MiSeq", "NextSeq1k2k")[sequencer_select]

#sequencing date of the run folder should match the RStudio project date
sequencing_date <- gsub("_.*", "", basename(here())) #YYYY-MM-DD

# temporary directory to hold the screen log files
tmp_screen_fp <- paste("~", ".tmp_screen", sequencer_type, "NS_SC2", basename(here()), sep = "/")

session_suffix <- tolower(paste(sequencer_type, "ns-sc2", basename(here()), sep = "-"))

# temporary directory to hold the sequencing run download
ec2_tmp_fp <- "~/tmp_bs_dl"

#file location of the nextera udi indices
barcode_fp <- file.path(dirname(here()), "aux_files", "illumina_references", "nextera-dna-udi-samplesheet-MiSeq-flex-set-a-d-2x151-384-samples.csv")

if(sequencing_date == "" | is.na(as.Date(sequencing_date, "%Y-%m-%d")) | nchar(sequencing_date) == 8) {
  stop(simpleError(paste0("Please use the 'YYYY-MM-DD' format for this RStudio project date. This date should correspond to the desired sequencing run date")))
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

###############
# Get run stats
###############

yymmdd <- gsub("^..|-", "", sequencing_date)
seq_folder_pattern <- "[0-9]*_[0-9]*_[0-9A-Z-]*"

run_cd <- NA
run_q30 <- NA
run_pf <- NA
run_error <- NA

if(run_uploaded_2_basespace) {

  unknown_sequencer_regex <- paste0(yymmdd, "_([M]{1}|[VH]{2})", seq_folder_pattern, "$")

  # Get the run id from BaseSpace
  bs_run <- cli_submit("bs", "list", c("runs", "-f csv")) %>%
    str_split(",") %>%
    do.call("rbind", .) %>%
    as.data.frame() %>%
    `colnames<-`(.[1, ]) %>%
    slice(-1) %>%
    filter(grepl(paste0("^", unknown_sequencer_regex), Name))

  if(nrow(bs_run) > 1) {
    warning(simpleWarning(paste0("\nThere are two sequencing runs that matched this date. Make sure you selected the correct sequencer!!!\n",
                                 "Currently, you are pulling the sequencing run from the ", sequencer_type, "\n\n")))

    #these Rscripts don't account for two runs that have the same sample types, processed on the same date, on both machines, and the samples need to be processed through the same pipeline
    sequencer_regex <- case_when(sequencer_type == "MiSeq" ~ "M",
                                 sequencer_type == "NextSeq1k2k" ~ "VH")

    intended_sequencing_folder_regex <- paste0(yymmdd, "_", sequencer_regex, seq_folder_pattern, "$")

    bs_run <- bs_run %>%
      filter(grepl(paste0("^", intended_sequencing_folder_regex), Name))

  }
  if (nrow(bs_run) == 0) {
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

  run_stats <- cli_submit("bs", "run", c("seqstats", "--id", bs_run_id)) %>%
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

}

#######################
# Load run sample sheet
#######################

run_samplesheet_fp <- list.files(here("metadata", "munge"), pattern = "SampleSheet", full.names = TRUE)[1]
samplesheet_exists <- file.exists(run_samplesheet_fp)

if(samplesheet_exists) {

  message("\n*****")
  message("There is already an existing SampleSheet.csv or SampleSheet_v2.csv in the metadata/munge directory")
  message("BaseSpace run will not be re-downloaded")
  message("*****")
  Sys.sleep(10)

} else {

  run_samplesheet_fp <- here("metadata", "munge", "SampleSheet.csv")

  if(run_uploaded_2_basespace) {

    temporary_seq_run_fp <- paste0(ec2_tmp_fp, "/", sequencing_run, "/")
    bs_dl_cmd <- paste("bs download runs --id", bs_run_id, "--output", temporary_seq_run_fp)

    if(have_AWS_EC2_SSH_access) {
      # Download the run from BaseSpace onto a running EC2 instance
      submit_screen_job(message2display = "Downloading sequencing run from BaseSpace",
                        ec2_login = ec2_hostname,
                        screen_session_name = paste("basespace-run-download", session_suffix, sep = "-"),
                        screen_log_fp = tmp_screen_fp,
                        command2run = bs_dl_cmd
      )

      check_screen_job(message2display = "Checking BaseSpace download job",
                       ec2_login = ec2_hostname,
                       screen_session_name = paste("basespace-run-download", session_suffix, sep = "-"),
                       screen_log_fp = tmp_screen_fp)

      # Download the SampleSheet from EC2 instance
      run_in_terminal(paste("scp",
                            paste0(ec2_hostname, ":", temporary_seq_run_fp, "SampleSheet.csv"),
                            run_samplesheet_fp),
                      command2print = paste(" [On", ec2_hostname, "instance]\n",
                                            "aws s3 cp", paste0(temporary_seq_run_fp, "SampleSheet.csv"),
                                            paste0("s3://test-environment/input/", session_suffix, "/"), "\n\n",
                                            "[On local computer]\n",
                                            "aws s3 cp", paste0("s3://test-environment/input/", session_suffix, "/SampleSheet.csv"),
                                            run_samplesheet_fp)
      )
    } else {

      dir.create(temporary_seq_run_fp, recursive = TRUE)
      run_in_terminal(bs_dl_cmd)
      file.copy(paste0(temporary_seq_run_fp, "SampleSheet.csv"), run_samplesheet_fp)
    }

  } else if (!run_uploaded_2_basespace) {

    if(sequencer_type == "MiSeq"){

      intended_miseq_folder_regex <- paste0(yymmdd, "_M", seq_folder_pattern)

      host_samplesheet_fp <- paste0(miseq_hostname, ":", intended_miseq_folder_regex, "/SampleSheet.csv")

    } else if (sequencer_type == "NextSeq1k2k") {

      intended_nextseq_folder_regex <- paste0(yymmdd, "_VH", seq_folder_pattern)
      nextseq_run_fp <- "/usr/local/illumina/runs/"

      sequencing_run <- system2("ssh", c("-tt", nextseq_hostname,
                                         shQuote(paste("cd", paste0(nextseq_run_fp, ";"),
                                                       "ls | grep", intended_nextseq_folder_regex, "| tr -d '\n'"),
                                                 type = "sh")),
                                stdout = TRUE)

      host_samplesheet_fp <- paste0(nextseq_hostname, ":", nextseq_run_fp, sequencing_run, "/SampleSheet*.csv")
    }

    message("\n\n\n*****")
    message("Transferring SampleSheet.csv locally")
    message("*****")
    Sys.sleep(5)

    scp_command <- paste("scp",
                         host_samplesheet_fp,
                         paste(run_samplesheet_fp, ";"),
                         "sleep 2")

    run_in_terminal(scp_command)
  }
  rstudioapi::executeCommand("activateConsole")
}

run_sample_sheet <- load_sample_sheet(run_samplesheet_fp)

instrument_type <- data.frame(values = unlist(run_sample_sheet$Header)) %>%
  mutate(col_names = gsub(",.*", "", values)) %>%
  mutate(col_names = gsub(" ", "_", col_names)) %>%
  mutate(values = gsub(".*,", "", values)) %>%
  filter(grepl("instrument_type|InstrumentPlatform|InstrumentType", col_names, ignore.case = TRUE)) %>%
  select(values) %>%
  pull()

read_length <- data.frame(values = unlist(run_sample_sheet$Reads)) %>%
  filter(!grepl("^Index[1|2]Cycles,", values)) %>%
  mutate(values = gsub(".*,", "", values)) %>%
  pull() %>%
  unique()

if(instrument_type != sequencer_type) {
  message("\n*****")
  message("The SampleSheet.csv for this ", sequencing_date, " run has the instrument set as ", instrument_type)
  message("The rest of this script will continue and processing this project as a ", instrument_type, " run")
  message("If this was not the correct sequencer used for this project, double check the sequencing date or select the appropriate sequencer_type in this Rscript")
  message("*****")

  Sys.sleep(10)
}

instrument_regex <- case_when(instrument_type == "MiSeq" ~ "M",
                              instrument_type == "NextSeq1k2k" ~ "VH")

if(!read_length %in% c(76, 151)) {
  stop(simpleError("The read length is not 76 or 151 bp. Check the sample sheet from the sequencing run folder"))
}

##################
# Load index sheet
##################

index_sheet_fp <- list.files(here("metadata", "munge"), pattern = ".xlsx", full.names = TRUE)

if(length(index_sheet_fp) == 0) {
  shared_index_fp <- max(list.files(file.path(shared_drive_fp, "Sequencing_files", "3_Sample_Sheets", "nasal_swabs", str_sub(sequencing_date, 1, 4)),
                                    pattern = "sequencing_metadata_sheet", full.names = TRUE))

  if(is.na(shared_index_fp)) {
    stop(simpleError("Files in the shared drive could not be found\nAre you connected to the shared drive?"))
  }

  #get date of index sheet
  sheet_date <- as.Date(str_extract(shared_index_fp, "[0-9-]{8}"), tryFormats = c("%y%m%d", "%m%d%y", "%m-%d-%y"))

  if(abs(as.Date(sequencing_date) - sheet_date) > 5) {
    stop(simpleError(paste0("The latest index sheet found is on ", sheet_date,
                            "\nThis date is more than 5 days from the date of this RStudio project",
                            "\nSomething must be wrong")))
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

###################################################################################
# Load the metadata sheet from epidemiologists and merge with sample metadata sheet
# Make sure these sheets are not uploaded to GitHub
###################################################################################

PHL_fp <- list.files(here("metadata", "extra_metadata"), pattern = "_filtered.xlsx", full.names = TRUE, recursive = TRUE)

PHL_data <- PHL_fp %>%
  lapply(function(x) read_excel_safely(x, "PHL")) %>%
  bind_rows()

#if PHL_data exists, do the following
# if it's just a blank sheet, just save sample_name and RLU column names
if(ncol(PHL_data) == 2) {

  PHL_data <- PHL_fp %>%
    lapply(function(x) read_excel_safely(x, "PHL")) %>%
    bind_rows() %>%
    rename(sample_name = "SPECIMEN_NUMBER") %>%
    #filter rows where sample_id is NA
    filter(!is.na(sample_name)) %>%
    #make these columns character vectors
    mutate(sample_name = as.character(sample_name),
           sample_collected_by = "Philadelphia Department of Public Health")

} else if (ncol(PHL_data) > 2) {

  PHL_data <- PHL_data %>%
    rename(sample_name = "SPECIMEN_NUMBER", sample_collection_date = "SPECIMEN_DATE", gender = "GENDER", DOB = "BIRTH_DATE") %>%
    mutate(sample_collection_date = as.Date(sample_collection_date, format = "%m/%d/%Y"),
           DOB = as.Date(DOB, format = "%m/%d/%Y"),
           gender = ifelse(is.na(gender), "Unknown", gender),
           sample_collected_by = "Philadelphia Department of Public Health") %>%
    select(any_of(phi_info), sample_collection_date, sample_collected_by, DOB, age, gender, RLU) %>%
    #filter rows where sample_id is NA
    filter(!is.na(sample_name)) %>%
    #make these columns character vectors
    mutate_at(vars(any_of(c(phi_info, "gender"))), as.character)

}

#if PHL_data is not empty, do the following
if(nrow(PHL_data) > 0) {

  PHL_data <- PHL_data %>%
    #filter empty columns
    select(where(function(x) any(!is.na(x)))) %>%
    select(!matches("^\\.\\.\\.")) %>%
    #use the first day of the week (starting on Monday) as the sample_collection_date
    mutate(sample_collection_date = as.Date(cut(as.POSIXct(sample_collection_date), "week"))) %>%
    mutate(host_age_bin = cut(age, breaks = c(0, 9, as.numeric(paste0(1:6, 9)), Inf),
                              labels = c("0 - 9", paste(seq(10, 60, by = 10), "-",as.numeric(paste0(1:6, 9))), "70+"),
                              include.lowest = TRUE)) %>%
    #don't include age because it may be PHI if included with zipcode and gender
    select(-c(age, DOB)) %>%
    as.data.frame()
}

###################################################################################
# Load the metadata sheet from epidemiologists and merge with sample metadata sheet
# Make sure these sheets are not uploaded to GitHub
###################################################################################

TU_data <- PHL_fp %>%
  lapply(function(x) read_excel_safely(x, "Temple")) %>%
  lapply(function(x) mutate(x, SPECIMEN_NUMBER = as.character(SPECIMEN_NUMBER))) %>%
  bind_rows()

if(ncol(TU_data) == 2) {

  TU_data <- PHL_fp %>%
    lapply(function(x) read_excel_safely(x, "Temple")) %>%
    bind_rows() %>%
    rename(sample_name = "SPECIMEN_NUMBER", CT = "ct value") %>%
    #filter rows where sample_id is NA
    filter(!is.na(sample_name)) %>%
    mutate(sample_name = as.character(sample_name),
           sample_collected_by = "Temple University")

} else if(ncol(TU_data) > 2) {

  TU_data <- TU_data %>%
    rename(sample_name = "SPECIMEN_NUMBER",
           sample_collection_date = "Collection_date",
           CT = "ct value", gender = "GENDER") %>%
    mutate(sample_collection_date = as.Date(sample_collection_date, format = "%m/%d/%Y"),
           sample_collected_by = "Temple University") %>%
    select(any_of(phi_info), sample_collection_date, sample_collected_by, CT, age, gender) %>%
    #filter rows where sample_id is NA
    filter(!is.na(sample_name)) %>%
    #make these columns character vectors
    mutate_at(vars(any_of(c(phi_info, "gender"))), as.character)
}

if(nrow(TU_data) > 0) {

  TU_data <- TU_data %>%
    #filter empty columns
    select(where(function(x) any(!is.na(x)))) %>%
    select(!matches("^\\.\\.\\.")) %>%
    #use the first day of the week (starting on Monday) as the sample_collection_date
    mutate(sample_collection_date = as.Date(cut(as.POSIXct(sample_collection_date), "week"))) %>%
    mutate(host_age_bin = cut(age, breaks = c(0, 9, as.numeric(paste0(1:6, 9)), Inf),
                              labels = c("0 - 9", paste(seq(10, 60, by = 10), "-",as.numeric(paste0(1:6, 9))), "70+"),
                              include.lowest = TRUE)) %>%
    #don't include age because it may be PHI if included with zipcode and gender
    select(-age) %>%
    as.data.frame()
}

################################
# Load the environmental samples
################################

ENV_fp <- max(list.files(here("metadata", "extra_metadata"), pattern = "environmental_samples.csv", full.names = TRUE))

if(!is.na(ENV_fp)) {

  ENV_data <- read_csv(ENV_fp) %>%
    #use the Tuesday day of the sequencing week as the sample_collection_date
    mutate(sample_collection_date = as.Date(cut(as.POSIXct(sequencing_date), "week")) + 1,
           sample_collected_by = "Philadelphia Department of Public Health") %>%
    select(sample_name, sample_collection_date, sample_collected_by, environmental_site) %>%
    #filter rows where sample_id is NA
    filter(!is.na(sample_name)) %>%
    #filter empty columns
    select(where(function(x) any(!is.na(x)))) %>%
    select(!matches("^\\.\\.\\."))
} else{

  ENV_data <- data.frame(environmental_site = NA)
}

########################
# Load all other samples
########################

other_sheets <- unique(unlist(lapply(PHL_fp, excel_sheets)))
other_sheets <- other_sheets[!grepl("PHL|Temple", other_sheets)]

other_data <- data.frame(sample_name = NA,
                         case_id = NA,
                         sample_collection_date = NA,
                         CT = NA_real_,
                         RLU = NA_real_,
                         gender = NA,
                         age = NA_integer_,
                         zip_char = NA_character_
)

for(sheet_name in other_sheets) {

  possible_sample_names <- "SPECIMEN_NUMBER"

  possible_case_id <- "cdms_id"

  possible_zip_char <- "zip_code"

  possible_ct_value <- "ct_value"

  other_sheet <- PHL_fp %>%
    lapply(function(x) read_excel_safely(x, sheet_name)) %>%
    bind_rows() %>%
    rename_at(vars(contains(possible_sample_names)),
              ~gsub(possible_sample_names, "sample_name", ., ignore.case = TRUE)) %>%
    mutate(sample_name = gsub("\\s", "", sample_name)) %>%
    rename_at(vars(contains(possible_case_id)),
              ~gsub(possible_case_id, "case_id", ., ignore.case = TRUE)) %>%
    mutate(case_id = gsub("\\s", "", case_id)) %>%
    rename_at(vars(contains(possible_zip_char)),
              ~gsub(possible_zip_char, "zip_char", ., ignore.case = TRUE)) %>%
    rename_at(vars(contains(possible_ct_value)),
              ~gsub(possible_ct_value, "CT", ., ignore.case = TRUE)) %>%
    mutate(host_age_bin = cut(age, breaks = c(0, 9, as.numeric(paste0(1:6, 9)), Inf),
                              labels = c("0 - 9", paste(seq(10, 60, by = 10), "-",as.numeric(paste0(1:6, 9))), "70+"),
                              include.lowest = TRUE)) %>%
    as.data.frame()

  other_data <- bind_rows(other_data, other_sheet) %>%
    select(-c(sample_type, age))

}

###########################
# Merge all metadata sheets
###########################

PHL_TU_merge <- PHL_data %>%
  bind_rows(TU_data) %>%
  bind_rows(ENV_data) %>%
  bind_rows(other_data) %>%
  mutate(sample_collection_date = as.character(sample_collection_date)) %>%
  select(-any_of(c("age", "DOB")))

cols2merge <- c("sample_name", "plate", "plate_row", "plate_col", "plate_coord")

#merge all the individual sheets
metadata_sheet <- merge(index_sheet, sample_info_sheet, by = cols2merge, all = TRUE, sort = FALSE) %>%
  #add in barcodes
  merge(barcodes, by = "idt_plate_coord", all.x = TRUE, sort = FALSE) %>%
  #merge the metadata
  merge(PHL_TU_merge, by = "sample_name", all.x = TRUE, sort = FALSE) %>%
  mutate(sample_id = gsub("_", "-", paste0("PHL2", "-", instrument_regex, "-", idt_plate_coord, "-", gsub("-", "", sequencing_date))),
         sequencing_date = sequencing_date,
         prj_descrip = prj_description,
         instrument_type = instrument_type,
         read_length = read_length,
         index_length = index_length,
         run_cd = run_cd,
         run_q30 = run_q30,
         run_pf = run_pf,
         run_error = run_error) %>%
  select(sample_id, everything()) %>%
  arrange(plate, plate_col, plate_row)

#####################################################################
# Fill in these columns in the metadata sheet if they were left blank
#####################################################################

fill_in_columns <- c("sample_type", "sample_collected_by", "PHL_sample_received_date",
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
                       "^[A-Z0-9][0-9]*$" = "Nasal swab", #allow the Temple specimen IDs to be any number, once it passes 9
                       "^WW-" = "Wastewater")

metadata_sheet <- metadata_sheet %>%
  mutate(sample_type = case_when(!(is.na(sample_type) | sample_type == "") ~ sample_type,
                                 sample_name %in% ENV_data$sample_name ~ "Environmental control",
                                 (is.na(sample_type) | sample_type == "") ~ multi_grep(named_sample_type, sample_name),
                                 TRUE ~ NA)) %>%
  mutate(sample_collected_by = case_when(!(is.na(sample_collected_by) | sample_collected_by == "") ~ sample_collected_by,
                                         grepl("Water control|Reagent control|Mock DNA positive control", sample_type) ~ "Philadelphia Department of Public Health",
                                         TRUE ~ NA)) %>%
  mutate(PHL_sample_received_date = case_when(!(is.na(PHL_sample_received_date) | as.character(PHL_sample_received_date) == "") ~ as.Date(PHL_sample_received_date),
                                              #if it's a wastewater sample without a date, throw an error
                                              sample_type == "Wastewater" ~ NA,
                                              #use Tuesday of the sequencing week if no date specified; older samples that are rerun should have a date manually added in on the sheet
                                              TRUE ~ as.Date(cut(as.POSIXct(sequencing_date), "week")) + 1)) %>%
  mutate(sample_collection_date = case_when(!(is.na(sample_collection_date) | as.character(sample_collection_date) == "") ~ as.Date(sample_collection_date),
                                            #if it's a wastewater sample without a date, use the PHL sample received date
                                            sample_type == "Wastewater" ~ PHL_sample_received_date,
                                            TRUE ~ NA)) %>%
  mutate(organism = case_when(!(is.na(organism) | organism == "") ~ organism,
                              sample_type == "Nasal swab" ~ "Severe acute respiratory syndrome coronavirus 2",
                              #WW sample has to be listed as metagenome even if targeted sequencing was used
                              sample_type == "Wastewater" ~ "Wastewater metagenome",
                              !is.na(sample_type) ~ sample_type,
                              TRUE ~ NA)) %>%
  mutate(host_scientific_name = case_when(!(is.na(host_scientific_name) | host_scientific_name == "") ~ host_scientific_name,
                                          sample_type == "Nasal swab" ~ "Homo sapiens",
                                          !is.na(sample_type) ~ "not applicable",
                                          TRUE ~ NA)) %>%
  mutate(host_disease = case_when(!(is.na(host_disease) | host_disease == "") ~ host_disease,
                                  sample_type == "Nasal swab" ~ "COVID-19",
                                  !is.na(sample_type) ~ "not applicable",
                                  TRUE ~ NA)) %>%
  mutate(isolation_source = case_when(!(is.na(isolation_source) | isolation_source == "") ~ isolation_source,
                                      sample_type == "Nasal swab" ~ "Clinical",
                                      sample_type == "Wastewater" ~ "Wastewater",
                                      grepl("control$", sample_type) ~ "Environmental",
                                      grepl("test", sample_type, ignore.case = TRUE) ~ "Test sample",
                                      TRUE ~ NA)) %>%
  mutate(requester = case_when(!(is.na(requester) | requester == "") ~ requester,
                               TRUE ~ epi_name)) %>%
  mutate(requester_email = case_when(!(is.na(requester_email) | requester_email == "") ~ requester_email,
                                     TRUE ~ epi_email)) %>%
  mutate(environmental_site = case_when(grepl("Water control|Reagent control|Mock DNA positive control", sample_type) ~ paste0(sample_name, " - ", plate_row, plate_col),
                                        grepl("Environmental control", sample_type) ~ paste0(environmental_site, " - ", plate_row, plate_col),
                                        TRUE ~ environmental_site))

main_sample_type <- unique(metadata_sheet$sample_type)[!grepl("control", unique(metadata_sheet$sample_type))]

if(length(main_sample_type) > 1) {
  stop(simpleError(paste0("This metadata sheet has more than one non-control sample type!\n",
                          "You may need to separate the metadata sheet and use the appropriate workflow for these samples types:\n",
                          paste0(main_sample_type, collapse = ", "))))
  sample_type_acronym <- "Mix"
}

if(!grepl("Nasal swab|Testing sample type", main_sample_type)) {
  stop(simpleError(paste0("The sample type included in the metadata sheet is not nasal swab or a test sample type!\n",
                          "This may not be the appropriate workflow for this run!\n")))
}

sample_type_acronym <- case_when(main_sample_type == "Testing sample type" ~ "Test",
                                 main_sample_type == "Nasal swab" ~ "NS",
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

######################################################################################
# Samples in epi metadata but we don't have the samples or they could not be extracted
######################################################################################

epi_sample_not_found <- PHL_data %>%
  select(any_of("sample_name")) %>%
  rbind(select(TU_data, any_of("sample_name")))

if(ncol(epi_sample_not_found > 0)) {

  epi_sample_not_found <- epi_sample_not_found %>%
    filter(!sample_name %in% metadata_sheet$sample_name) %>%
    filter(!sample_name %in% remove_sample_from_samplesheets) %>%
    pull() %>%
    str_sort()

  #throw error
  if(length(epi_sample_not_found) > 0) {

    message("\n*****")
    message("These samples were found in the epidemiologists metadata sheet but not in our sample sheet. Something might be wrong!")
    message("Check email to see if these samples could not be located by the receiving department")
    message("Otherwise, these samples may have had an issue during extraction. Send wet lab scientists these sample names to check:")
    message("*****")

    stop(simpleError(paste0(epi_sample_not_found, collapse = ", ")))
  }

}

missing_metadata_non_ctrl_samples <- metadata_sheet %>%
  filter(!grepl("control", sample_type))

missing_sample_date <- missing_metadata_non_ctrl_samples %>%
  filter(is.na(sample_collection_date)) %>%
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

missing_sample_RLU <- missing_metadata_non_ctrl_samples %>%
  filter(sample_collected_by == "Philadelphia Department of Public Health") %>%
  filter(is.na(RLU)) %>%
  select(sample_name)

missing_sample_CT <- missing_metadata_non_ctrl_samples %>%
  filter(sample_collected_by == "Temple University") %>%
  filter(is.na(CT)) %>%
  select(sample_name)

missing_metadata_samples <- rbind(missing_sample_RLU, missing_sample_CT)

#show warning
if(nrow(missing_metadata_samples) > 0){
  message("\n*****")
  message("These non-control samples are in the sample sheet but are missing RLU or CT values")
  message("They may also be GeneXpert samples. These samples should be double checked:")
  message(paste0(pull(missing_metadata_samples), collapse = ", "))
  message("*****")

  Sys.sleep(5)
}

#############
# Check sheet
#############

#throw error if missing these columns
for(x in c(cols2merge, "sample_id",
           "idt_set", "idt_plate_row", "idt_plate_col", "idt_plate_coord",
           "index", "index2", "UDI_Index_ID", "I7_Index_ID", "I5_Index_ID",
           "sequencing_date", "prj_descrip", "instrument_type", "read_length", "index_length",
           "environmental_site", "sample_collection_date", "gender", "zip_char")) {
  if(!grepl(paste0(colnames(metadata_sheet), collapse = "|"), x)) {
    stop(simpleError(paste0("\nMissing column [", x, "] in the metadata sheet template!!!")))
  }
}

#remove columns that are all empty
metadata_sheet <- metadata_sheet %>%
  select(where(function(x) any(!is.na(x))))

#add back these columns as NA if missing (needed for the report and for seqsender)
for(x in c("qubit_conc_ng_ul", "sample_collection_date", "host_age_bin", "gender")) {
  if(!grepl(paste0(colnames(metadata_sheet), collapse = "|"), x)) {
    metadata_sheet[[x]] <- NA
  }
}

#check lowest date of sample collection
if(min(as.Date(metadata_sheet$sample_collection_date[!is.na(metadata_sheet$sample_collection_date)])) < seq(as.Date(sequencing_date), length=2, by='-6 month')[2]){
  message("Earliest sample collection date in this run:")
  message(min(as.Date(metadata_sheet$sample_collection_date[!is.na(metadata_sheet$sample_collection_date)])))
  stop(simpleError(paste0("\nSome samples have collection dates more than 6 months ago. Investigate!!")))
}

print('What do the sample_id look like?')
print(unique(metadata_sheet$sample_id))

print('Which lanes are sequenced?')
print(unique(metadata_sheet$lane))

print('Are the barcode columns unique?')
if(length(unique(metadata_sheet$idt_plate_coord)) != dim(metadata_sheet)[1]) {
  stop(simpleError("Barcode positions are not unique!"))
}
print(length(unique(metadata_sheet$idt_plate_coord)) == dim(metadata_sheet)[1])

print('Number of barcodes?')
print(length(unique(paste0(metadata_sheet$index, metadata_sheet$index2))))
print('Number of samples?')
print(length(unique(metadata_sheet$sample_id)))
if(length(unique(paste0(metadata_sheet$index, metadata_sheet$index2))) != length(unique(metadata_sheet$sample_id))) {
  stop(simpleError("Differing number of samples and barcodes!"))
}

print('Are the barcodes unique?')
print(length(unique(paste0(metadata_sheet$index, metadata_sheet$index2))) == dim(metadata_sheet)[1])

print('Are the sample names unique?')
print(length(unique(metadata_sheet$sample_id)) == dim(metadata_sheet)[1])

print('Are all the forward primers found?')
print(sum(is.na(metadata_sheet$index)) == 0)

print('Are all the reverse primers found?')
print(sum(is.na(metadata_sheet$index2)) == 0)
if(sum(is.na(c(metadata_sheet$index, metadata_sheet$index2))) != 0) {
  stop(simpleError("Either the forward or reverse primers are NA!"))
}

print('Do all the sampleIDs start with a letter?')
print(all(grepl("^[A-Za-z]", metadata_sheet$sample_id)))
if(!all(grepl("^[A-Za-z]", metadata_sheet$sample_id))) {
  stop(simpleError("Some Sample IDs do not start with a letter!"))
}

print('Are there periods, underscores, or space characters in the SampleID?')
print(any(grepl(" |_|\\.", metadata_sheet$sample_id)))
if(any(grepl(" |_|\\.", metadata_sheet$sample_id))) {
  stop(simpleError("There are spaces, underscores, or periods in the Sample IDs! Please fix"))
}

####################
# Write sample sheet
####################

samp_sheet_2_write <- metadata_sheet %>%
  filter(!sample_name %in% remove_sample_from_samplesheets) %>%
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

################################
# Write sheet to metadata folder
################################

#does not contain PHI and accession numbers
metadata_sheet %>%
  filter(!sample_name %in% remove_sample_from_samplesheets) %>%
  select(-c(I7_Index_ID, I5_Index_ID, any_of(phi_info))) %>%
  # NextSeq runs have index2 sequences in the reverse complement of the ones listed in the reference barcode sheet
  # however, the sample sheet used for demultiplexing needs to be the same orientation as the reference barcode sheet
  # when BCLConvert demultiplexes a NextSeq run, it will automatically reverse complement index2
  # results from the pipeline will refer to the reverse complement of index2, so this should be updated in the metadata sheet
  mutate(index2 = ifelse(instrument_type == "NextSeq1k2k", reverse_complement(index2), index2)) %>%
  write_csv(file = here("metadata", paste0(sequencing_date, "_", prj_description, "_metadata.csv")))

#contains PHI and accession numbers
metadata_sheet %>%
  filter(!sample_name %in% remove_sample_from_samplesheets) %>%
  select(sample_id, any_of(phi_info)) %>%
  write_csv(file = here("metadata", paste0(sequencing_date, "_", prj_description, "_PHI.csv")))
