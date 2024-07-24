library(here)
library(dplyr)
library(tidyr)
library(readxl)
library(readr)
library(stringr)

#This Rscript currently generates the SampleSheet for demultiplexing a run using the BCL Convert program
#https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl_convert/bcl-convert-v3-7-5-software-guide-1000000163594-00.pdf

##############
# Manual input
##############

run_uploaded_2_basespace <- TRUE # set this to true if the run was uploaded to BaseSpace and the data was not manually transferred to a local folder

sequencer_select <- 1 # set variable as 1 for MiSeq or 2 for NextSeq

have_AWS_EC2_SSH_access <- FALSE

remove_sample_from_bcl_samplesheet <- c("") #add in sample names to remove from demultiplexing

# temporary directory to hold the sequencing run download
ec2_tmp_fp <- "~/tmp_bs_dl/"

prj_description <- "COVIDSeq" #no spaces, should be the same as the R project

index_length <- "10"

####################
# Selected variables
####################

sequencer_type <- c("MiSeq", "NextSeq1k2k")[sequencer_select]

#sequencing date of the run folder should match the RStudio project date
sequencing_date <- gsub("_.*", "", basename(here())) #YYYY-MM-DD

#file location of the wastewater metadata
ww_meta_fp <- file.path(dirname(here()), "aux_files", "data_submission", "dcipher", "wastewater_specific_metadata.csv")

#file location of the nextera udi indices
barcode_fp <- file.path(dirname(here()), "aux_files", "illumina_references", "nextera-dna-udi-samplesheet-MiSeq-flex-set-a-d-2x151-384-samples.csv")

if(sequencing_date == "" | is.na(as.Date(sequencing_date, "%Y-%m-%d")) | nchar(sequencing_date) == 8) {
  stop(simpleError(paste0("Please use the 'YYYY-MM-DD' format for this RStudio project date. This date should correspond to the desired sequencing run date")))
}
if (prj_description == "") {
  stop (simpleError(paste0("The project description variable in ", here("code"), "/2_generate_barcodes_IDT.R is empty. Make sure it is set to the correct project workflow")))
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

#######################
# Load run sample sheet
#######################

samplesheet_exists <- file.exists(here("metadata", "munge", "SampleSheet.csv"))

sequencing_folder_regex <- paste0(gsub("^..|-", "", sequencing_date), "_([M]{1}|[VH]{2})[0-9]*_[0-9]*_[0-9A-Z-]*$")

run_cd <- NA
run_q30 <- NA
run_pf <- NA
run_error <- NA

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

  run_stats <- cli_submit("bs", "run", c("seqstats", "--id", bs_run_id))

  run_cd <- run_stats %>%
    list(run_stats = .) %>%
    as.data.frame() %>%
    filter(grepl("SequencingStatsCompact.ClusterDensity", run_stats)) %>%
    # cluster density seems to be reported in the millions. If there is no scientific notation, divide by 1000
    mutate(run_stats = gsub(".*\\| | .*", "", run_stats),
           run_stats = ifelse(grepl("e", run_stats),
                              as.numeric(gsub("e.*", "", run_stats))*1000,
                              as.numeric(run_stats)/1000)) %>%
    pull()

  run_q30 <- run_stats %>%
    list(run_stats = .) %>%
    as.data.frame() %>%
    filter(grepl("SequencingStatsCompact.PercentGtQ30 ", run_stats)) %>%
    mutate(run_stats = gsub(".*\\| | .*", "", run_stats),
           run_stats = as.numeric(run_stats)/100) %>%
    pull()

  run_pf <- run_stats %>%
    list(run_stats = .) %>%
    as.data.frame() %>%
    filter(grepl("SequencingStatsCompact.PercentPf", run_stats)) %>%
    mutate(run_stats = gsub(".*\\| | .*", "", run_stats),
           run_stats = as.numeric(run_stats)) %>%
    pull()

  run_error <- run_stats %>%
    list(run_stats = .) %>%
    as.data.frame() %>%
    filter(grepl("SequencingStatsCompact.ErrorRate ", run_stats)) %>%
    mutate(run_stats = gsub(".*\\| | .*", "", run_stats),
           run_stats = as.numeric(run_stats)/100) %>%
    pull()

  if(!samplesheet_exists) {

    temporary_seq_run_fp <- paste0(ec2_tmp_fp, sequencing_run, "/")
    bs_dl_cmd <- paste("bs download runs --id", bs_run_id, "--output", temporary_seq_run_fp)

    if(have_AWS_EC2_SSH_access) {
      # Download the run from BaseSpace onto a running EC2 instance
      submit_screen_job(message2display = "Downloading sequencing run from BaseSpace",
                        ec2_login = ec2_hostname,
                        screen_session_name = "basespace-run-download",
                        command2run = bs_dl_cmd
      )

      check_screen_job(message2display = "Checking BaseSpace download job",
                       ec2_login = ec2_hostname,
                       screen_session_name = "basespace-run-download")

      # Download the SampleSheet from EC2 instance
      run_in_terminal(paste("scp",
                            paste0(ec2_hostname, ":", temporary_seq_run_fp, "SampleSheet.csv"),
                            here("metadata", "munge")),
                      command2print = paste(" [On", ec2_hostname, "instance]\n",
                                      "aws s3 cp", paste0(temporary_seq_run_fp, "SampleSheet.csv"),
                                      paste0("s3://test-environment/input/", sequencing_date, "/"), "\n\n",
                                      "[On local computer]\n",
                                      "aws s3 cp", paste0("s3://test-environment/input/", sequencing_date, "/SampleSheet.csv"),
                                      here("metadata", "munge/"))
      )

      rstudioapi::executeCommand('activateConsole')
    } else {

      dir.create(temporary_seq_run_fp, recursive = TRUE)

      run_in_terminal(bs_dl_cmd)

      rstudioapi::executeCommand('activateConsole')

      file.copy(paste0(temporary_seq_run_fp, "SampleSheet.csv"), here("metadata", "munge", "SampleSheet.csv"))
    }
  }
} else if (!run_uploaded_2_basespace) {

  if(sequencer_type == "MiSeq"){

    #get the local run folder
    run_folder <- sequencing_folder_fp %>%
      list.files(full.names = T) %>%
      data.frame(filenames = .) %>%
      filter(grepl(format(as.Date(sequencing_date), "%y%m%d"), filenames)) %>%
      filter(grepl(sequencing_folder_regex, filenames)) %>%
      filter(!grepl("\\.tar\\.gz$|\\.md5$", filenames)) %>%
      pull()

    sequencing_run <- basename(run_folder)

    if(!samplesheet_exists) {

      message("\n\n\n*****")
      message("Copying MiSeq SampleSheet.csv")
      message("*****")
      Sys.sleep(5)

      sample_sheet_copy <- file.copy(file.path(run_folder, "SampleSheet.csv"), here("metadata", "munge", "SampleSheet.csv"))

      if(length(sample_sheet_copy) == 0) {
        stop(simpleError(paste("\nCould not find SampleSheet.csv from MiSeq run", sequencing_date, "in", sequencing_folder_fp)))
      }
    }

  } else if (sequencer_type == "NextSeq1k2k") {

    run_folder_pattern <- paste0(gsub("^..|-", "", sequencing_date))
    nextseq_run_fp <- "/usr/local/illumina/runs/"
    sequencing_run <- system2("ssh", c("-tt", nextseq_hostname,
                                       shQuote(paste("cd", paste0(nextseq_run_fp, ";"),
                                                     "ls | grep", paste0("^", run_folder_pattern, "_VH"), "| tr -d '\n'"),
                                               type = "sh")),
                              stdout = TRUE)

    if(length(sequencing_run) == 0) {
      stop(simpleError(paste("\nCould not find NextSeq run with date", sequencing_date, "on the sequencer")))
    }

    if(!samplesheet_exists) {

      message("\n\n\n*****")
      message("Transferring NextSeq SampleSheet.csv")
      message("*****")
      Sys.sleep(5)

      scp_command <- paste("scp",
                           paste0(nextseq_hostname, ":", nextseq_run_fp, sequencing_run, "/SampleSheet.csv"),
                           here("metadata", "munge", "SampleSheet.csv"))

      run_in_terminal(scp_command)
      rstudioapi::executeCommand('activateConsole')
    }
  }
}

run_samplesheet_fp <- here("metadata", "munge", "SampleSheet.csv")
run_sample_sheet <- load_sample_sheet(run_samplesheet_fp)

instrument_type <- data.frame(values = unlist(run_sample_sheet$Header)) %>%
  mutate(col_names = gsub(",.*", "", values)) %>%
  mutate(col_names = gsub(" ", "_", col_names)) %>%
  mutate(values = gsub(".*,", "", values)) %>%
  filter(grepl("instrument_type|InstrumentPlatform", col_names)) %>%
  select(values) %>%
  pull()

read_length <- data.frame(values = unlist(run_sample_sheet$Reads)) %>%
  filter(!grepl("^Index[1|2]Cycles,", values)) %>%
  mutate(values = gsub(".*,", "", values)) %>%
  pull() %>%
  unique()

if(instrument_type != sequencer_type) {
  message("\n*****")
  message("The downloaded SampleSheet.csv for this run on ", sequencing_date, " has the instrument set as the ", instrument_type)
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

#####################
# Load metadata sheet
#####################

metadata_input_fp <- list.files(here("metadata", "munge"), pattern = ".xlsx", full.names = TRUE)

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
      stop (simpleError("The sample metadata file from the wet lab scientists may be missing or mis-formatted"))
    }
  )
}

index_sheet <- read_sheet(metadata_input_fp, "Index")
sample_info_sheet <- read_sheet(metadata_input_fp, "Sample Info")

#######################################################
# Load the wastewater metadata sheet from the ddPCR run
#######################################################

ddPCR_fp <- list.files(here("metadata", "extra_metadata"), pattern = "_filtered.csv", full.names = TRUE, recursive = TRUE)

ddPCR_data <- read_csv(ddPCR_fp, show_col_types = FALSE)

# if ddPCR_data does not exist, create an empty dataframe
if(ncol(ddPCR_data) == 0) {

  ddPCR_data <- data.frame(
    sample_group = NA_character_,
    sample_collect_date = NA
  )

  # if ddPCR_data exists, reformat the columns
} else {

  ddPCR_data <- ddPCR_data %>%
    select(where(function(x) any(!is.na(x))),
           !matches("^\\.\\.\\.")) %>%
    rename(any_of(c(sample_group = "sample_id"))) %>%
    mutate(sample_group = as.character(sample_group),
           sample_collect_date = as.Date(sample_collect_date, tryFormats = c("%Y-%m-%d", "%m/%d/%y", "%m/%d/%Y"))) %>%
    as.data.frame()

}

################################
# Load the environmental samples
################################

ENV_fp <- max(list.files(here("metadata", "extra_metadata"), pattern = "environmental_samples.csv", full.names = TRUE))

if(!is.na(ENV_fp)) {

  ENV_data <- read_csv(ENV_fp) %>%
    #use the Tuesday of the sequencing week as the sample_collect_date
    mutate(sample_collect_date = as.Date(cut(as.POSIXct(sequencing_date), "week")) + 1) %>%
    select(sample_group = sample_name, sample_collect_date, environmental_site) %>%
    #filter rows where sample_id is NA
    filter(!is.na(sample_group)) %>%
    #filter empty columns
    select(where(function(x) any(!is.na(x))),
           !matches("^\\.\\.\\."))
} else{

  ENV_data <- data.frame(environmental_site = NA_character_)
}

extra_metadata_merge <- ddPCR_data %>%
  bind_rows(ENV_data) %>%
  unique()

#########################################
# Merge read in sequencing metadata sheet
#########################################

sequencing_controls <- c("Water control", "Reagent control", "Mock DNA positive control")
sample_group_controls <- c("PBS", "oldWW", "ZeptoSC2")
sample_group_sites <- c("NorthEast", "SouthEast", "SouthWest")

cols2merge <- c("sample_name", "plate", "plate_row", "plate_col", "plate_coord")

#merge all the sequencing sheets
metadata_sheet <- merge(index_sheet, sample_info_sheet, by = cols2merge, all = TRUE, sort = FALSE) %>%
  #add in barcodes
  merge(barcodes, by = "idt_plate_coord", all.x = TRUE, sort = FALSE) %>%
  arrange(plate, plate_col, plate_row) %>%
  rename(any_of(c(sample_collect_date = "sample_collection_date"))) %>%
  #find sample group from sample_name
  mutate(sample_group = gsub("^WW-([0-9]{4})-([0-9]{2})-([0-9]{2})-|^WW-|-Rep.*$", "", sample_name),
         # sample_group = as.character(lapply(sample_name,
         #                                    function(x) unlist(lapply(c(sample_group_controls, sample_group_sites),
         #                                                              function(y) y[grepl(y, x)])))),
         sample_group = ifelse((sample_group == "character(0)" | sample_group == ""), NA_character_, sample_group),
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
         run_error = run_error) %>%
  select(sample_id, everything())

#####################################################################
# Fill in these columns in the metadata sheet if they were left blank
#####################################################################

fill_in_columns <- c("sample_type", "sample_collected_by", "sample_collect_date",
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

extra_cols2merge <- c("sample_group", "sample_collect_date")

metadata_sheet <- metadata_sheet %>%
  mutate(sample_type = case_when(!(is.na(sample_type) | sample_type == "") ~ sample_type,
                                 sample_name %in% ENV_data$sample_name ~ "Environmental control",
                                 (is.na(sample_type) | sample_type == "") ~ multi_grep(named_sample_type, sample_name),
                                 TRUE ~ NA),
         sample_collected_by = case_when(!(is.na(sample_collected_by) | sample_collected_by == "") ~ sample_collected_by,
                                         TRUE ~ "Philadelphia Water Department"),
         sample_collect_date = case_when(!(is.na(sample_collect_date) | as.character(sample_collect_date) == "") ~ as.character(sample_collect_date),
                                         #if sample collect date column is not available, grab the date from the sample_name
                                         grepl("^WW-([0-9]{4})-([0-9]{2})-([0-9]{2})-", sample_name) ~ as.character(gsub("^(WW)-([0-9]{4}-[0-9]{2}-[0-9]{2})-(.*)", "\\2", sample_name)),
                                         #if it's a wastewater sample without a date or does not start with WW, throw an error
                                         sample_type == "Wastewater" ~ NA,
                                         #use Tuesday of the sequencing week if no date specified; older samples that are rerun should have a date manually added in on the sheet
                                         TRUE ~ as.character(as.Date(cut(as.POSIXct(sequencing_date), "week")) + 1)),
         sample_collect_date = as.Date(sample_collect_date),
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
                                  TRUE ~ gsub("^WW-|^Test", "", uniq_sample_name)),
         ww_group = case_when(grepl(paste0(sequencing_controls, collapse = "|"), sample_type) ~ sample_type,
                              grepl(paste0(sample_group_controls, collapse = "|"), uniq_sample_name) ~ "Wastewater control",
                              TRUE ~ "Wastewater sample")) %>%
  merge(extra_metadata_merge, by = extra_cols2merge, all.x = TRUE, sort = FALSE) %>%
  mutate(environmental_site = case_when(grepl(paste0(sequencing_controls, collapse = "|"), sample_type) ~ paste0(sample_name, " - ", plate_row, plate_col),
                                 grepl("Environmental control", sample_type) ~ paste0(environmental_site, " - ", plate_row, plate_col),
                                 TRUE ~ environmental_site)) %>%
  arrange(plate, plate_col, plate_row)

##########################
# Check the metadata sheet
##########################

main_sample_type <- unique(metadata_sheet$sample_type)[!grepl("control", unique(metadata_sheet$sample_type))]

if(length(main_sample_type) > 1) {
  stop(simpleError(paste0("This metadata sheet has more than one non-control sample type!\n",
                          "You may need to separate the metadata sheet and use the appropriate workflow for these samples types:\n",
                          paste0(main_sample_type, collapse = ", "))))
  sample_type_acronym <- "Mix"
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

#######################################################################################################
# Check for the samples in the ddPCR metadata sheet that didn't appear in the sequencing metadata sheet
#######################################################################################################

if(ncol(ddPCR_data) > 0) {

  ddPCR_sample_not_found <- ddPCR_data %>%
    filter(grepl(paste0(sample_group_sites, collapse = "|"), sample_group)) %>%
    merge(metadata_sheet, by = extra_cols2merge, all.x = TRUE) %>%
    filter(is.na(sample_id)) %>%
    mutate(date_n_name = paste(sample_collect_date, "-", sample_group)) %>%
    select(date_n_name) %>%
    pull() %>%
    str_sort()

  #throw error
  if(length(ddPCR_sample_not_found) > 0) {

    message("\n*****")
    message("These samples were found in the ddPCR metadata sheet but not in the sequencing sample sheet")
    message("Double check that the correct samples were sequenced:")
    message("*****")
    stop(simpleError(paste0(ddPCR_sample_not_found, collapse = ", ")))

  }
}

missing_metadata_non_ctrl_samples <- metadata_sheet %>%
  filter(!grepl("control", sample_type))

missing_sample_date <- missing_metadata_non_ctrl_samples %>%
  filter(is.na(sample_collect_date)) %>%
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

missing_ddPCR_data <- missing_metadata_non_ctrl_samples %>%
  filter(if_any(ends_with("_gc/L"), is.na)) %>%
  select(sample_name)

#show warning
if(nrow(missing_ddPCR_data) > 0){
  message("\n\n\n*****")
  message("These non-control samples are in the sample sheet but are missing ddPCR values")
  message("Check to see if the ddPCR has been performed for these samples:")
  message(paste0(pull(missing_ddPCR_data), collapse = ", "))
  message("*****")

  Sys.sleep(15)
}

#check lowest date of sample collection
oldest_ww_date <- metadata_sheet %>%
  filter(!grepl("control", ww_group)) %>%
  select(sample_collect_date) %>%
  pull() %>%
  min()

if(oldest_ww_date < seq(as.Date(sequencing_date), length=2, by='-4 month')[2]){
  warning(simpleError(paste0("\nSome samples have collection dates more than 4 months ago!!")))
}

#############
# Check sheet
#############

#throw error if missing these columns
for(x in c(cols2merge, "sample_id",
           "idt_set", "idt_plate_row", "idt_plate_col", "idt_plate_coord",
           "index", "index2", "UDI_Index_ID", "I7_Index_ID", "I5_Index_ID",
           "sequencing_date", "prj_descrip", "instrument_type", "read_length", "index_length",
           "environmental_site", "sample_collect_date")) {
  if(!grepl(paste0(colnames(metadata_sheet), collapse = "|"), x)) {
    stop(simpleError(paste0("\nMissing column [", x, "] in the metadata sheet template!!!")))
  }
}

#remove columns that are all empty
metadata_sheet <- metadata_sheet %>%
  select(where(function(x) any(!is.na(x))))

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
  filter(!sample_name %in% remove_sample_from_bcl_samplesheet) %>%
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
  mutate(index2 = ifelse(instrument_type == "NextSeq1k2k", reverse_complement(index2), index2)) %>%
  rbind(merged_samples_metadata_sheet) %>%
  write_csv(file = here("metadata", paste0(sequencing_date, "_", prj_description, "_metadata.csv")))

bclconvert_output_final_path <- paste(s3_fastq_bucket, sequencing_date, sample_type_acronym, prj_description, "processed_bclconvert", sequencing_run, sep = "/")

#write the sample sheet for merging nextflow script
metadata_sheet %>%
  select(fastq = sample_id, uniq_sample_name) %>%
  mutate(fastq_1 = paste0(bclconvert_output_final_path, "/", fastq, "_S", row_number(), "_R1_001.fastq.gz"),
         fastq_2 = paste0(bclconvert_output_final_path, "/", fastq, "_S", row_number(), "_R2_001.fastq.gz")) %>%
  merge(select(merged_samples_metadata_sheet, sample_id, uniq_sample_name), by = "uniq_sample_name", all = TRUE) %>%
  filter(!is.na(sample_id)) %>%
  select(sample_id, fastq_1, fastq_2) %>%
  write_csv(file = here("metadata", "munge",
                        tolower(paste(sequencing_date, instrument_type, sample_type_acronym, prj_description, "nf_concat_fastq_samplesheet.csv", sep = "_"))))
