library(here)
library(dplyr)
library(readxl)
library(readr)
library(stringr)
library(openxlsx)

#This Rscript filters out low RLU values from the metadata sheet received from the epidemiologists
#Send this filtered sheet to the epidemiologists and scientists

###################
# Default variables
###################

# number of unspecified environmental swabs to add to plate
enviro_number <- 10

################
# Load functions
################

#this file needs to sit in a [aux_files/r_scripts/functions] directory path above this project directory
tryCatch(
  {
    source(file.path(dirname(here()), "aux_files", "r_scripts", "functions", "R_all_functions_v3.R"))
  },
  error = function(e) {
    stop (simpleError("The R_all_functions_v3.R file needs to sit in a [aux_files/r_scripts/functions] directory path above this project directory"))
  }
)

#############
# Load config
#############

#this file needs to sit in a [aux_files/r_scripts/config] directory path above this project directory
tryCatch(
  {
    source(file.path(dirname(here()), "aux_files", "r_scripts", "config", "config_variables.R"))
  },
  error = function(e) {
    stop (simpleError("The config_variables.R file needs to sit in a [aux_files/r_scripts/config] directory path above this project directory"))
  }
)

###################################################
# Load the RLU report
# Make sure these sheets are not uploaded to GitHub
###################################################

# Look for this harvest report in extra_metadata folder
harvest_fn <- "COVID_SEQ_[0-9]*.csv"
harvest_fp <- list.files(here("metadata", "extra_metadata"), pattern = paste0("^", harvest_fn, "$"), full.names = TRUE)

# If the harvest file does not exist, grab it from the shared drive
# the harvest report is automatically generated each Monday and deposited onto the shared drive
if(length(harvest_fp) == 0) {

  #path of latest harvest report
  shared_harvest_fp <- max(list.files(file.path(shared_drive_fp, "Sequencing_harvest_reports"), pattern = paste0("^", harvest_fn, "$"), full.names = TRUE))

  #get date of harvest report
  date_harvest <- as.Date(str_extract(shared_harvest_fp, "[0-9]{6}"), tryFormats = c("%y%m%d"))

  if((Sys.Date() - date_harvest) > 5) {
    message(paste("The most recent COVIDSeq COPIA report found is "), date_harvest)
  }

  harvest_fp <- here("metadata", "extra_metadata", basename(shared_harvest_fp))

  file.copy(shared_harvest_fp, harvest_fp)

}

harvest_data <- read_csv(harvest_fp) %>%
  select_all(~gsub(" ", "_", .)) %>%
  mutate(Result = case_when(grepl("^positive$|^detected$", Result, ignore.case = TRUE) ~ "Positive",
                            grepl("^negative$|^not detected$|^presumptive negative$", Result, ignore.case = TRUE) ~ "Negative",
                            TRUE ~ NA))

if(any(is.na(harvest_data$Result))) {

  na_result <- harvest_data %>%
    filter(is.na(Result)) %>%
    select(Result) %>%
    unique() %>%
    pull()

  stop(simpleError(paste("The COPIA report has a new value in the result column that has not been previously captured.",
                         "Please add this new value to the template:",
                         paste(na_result, collapse = ", "), sep = "\n")))
}

RLU_data <- harvest_data %>%
  filter(Test_Name == "SARSCoV2-1") %>%
  rename(SPECIMEN_NUMBER = "Sample_ID", SPECIMEN_DATE = "Collection_Date", RLU = "Result", BIRTH_DATE = "Patient_DOB_(MM/dd/YYYY)") %>%
  select(SPECIMEN_NUMBER, SPECIMEN_DATE, BIRTH_DATE, RLU) %>%
  mutate(BIRTH_DATE = as.Date(BIRTH_DATE, format = "%m/%d/%Y"), SPECIMEN_DATE = as.Date(SPECIMEN_DATE, format = "%m/%d/%Y")) %>%
  #filter rows where sample_id is NA
  filter(!is.na(SPECIMEN_NUMBER))

#################
# Get POC samples
#################

point_of_care_samples <- harvest_data %>%
  filter(Test_Name == "POCT 4 Plex Sars-CoV-2") %>%
  select(Sample_ID) %>%
  pull()

###################################################
# Get GeneXpert samples
###################################################

GX_samples <- harvest_data %>%
  filter(Test_Name == "GeneXpert 4Plex Sars-CoV-2") %>%
  select(Sample_ID) %>%
  pull()

###################################################
# Get Respiratory Panel samples
###################################################

Resp_samples <- harvest_data %>%
  filter(Test_Name == "Severe Acute Respiratory Syndrome Coronavirus 2 (SARS-CoV-2)") %>%
  select(Sample_ID) %>%
  pull()

###################################################
# Get samples from Medical Examiner Office
###################################################

MEO_samples <- harvest_data %>%
  filter(Test_Name == "MEO POC Test Result") %>%
  select(Sample_ID) %>%
  pull()

########################################################################
# Load the metadata sheet from epidemiologists and merge with RLU values
########################################################################

PHL_all_fp <- list.files(here("metadata", "extra_metadata"), pattern = "PHLspecimens.*.xlsx", full.names = TRUE)
PHL_fp <- PHL_all_fp[!grepl("_filtered.xlsx$", PHL_all_fp)]

PHL_data <- PHL_fp %>%
  lapply(function(x) read_excel_safely(x, sheet = "PHL", skip_row = 1) %>%
           cbind(filename = x) %>%
           as.data.frame(stringsAsFactors = FALSE)) %>%
  bind_rows()

if(ncol(PHL_data) == 1) {
  PHL_data <- data.frame(SPECIMEN_NUMBER = "", RLU = "")
} else {
  PHL_data <- PHL_data %>%
    #filter rows where sample id is NA
    filter(!is.na(SPECIMEN_NUMBER)) %>%
    merge(RLU_data, by = c("SPECIMEN_NUMBER", "BIRTH_DATE", "SPECIMEN_DATE"), all.x = TRUE) %>%
    select(SPECIMEN_DATE, FacCode, agecoll, case_id, SPECIMEN_NUMBER,
           FIRST_NAME, LAST_NAME, BIRTH_DATE, age, zip_char, GENDER,
           breakthrough_case, death, hospitalized, outbreak, priority, RLU, filename) %>%
    mutate(SPECIMEN_DATE = format(SPECIMEN_DATE, "%m/%d/%Y"), BIRTH_DATE = format(BIRTH_DATE, "%m/%d/%Y"))
}

if(any(is.na(PHL_data[PHL_data$SPECIMEN_NUMBER %in% RLU_data$SPECIMEN_NUMBER, "RLU"]))) {

  no_RLU <- PHL_data %>%
    filter(SPECIMEN_NUMBER %in% RLU_data$SPECIMEN_NUMBER) %>%
    filter(is.na(RLU)) %>%
    select(SPECIMEN_NUMBER) %>%
    pull()

  stop(simpleError(paste("Serious error! These samples have an RLU value but did not get added to PHL_data", no_RLU, sep = "\n")))

}

########################
# Where samples are from
########################

missing_sample_with_RLU <- PHL_data %>%
  filter(!is.na(RLU)) %>%
  filter(SPECIMEN_NUMBER %in% point_of_care_samples) %>%
  select(SPECIMEN_NUMBER) %>%
  pull()

if(length(missing_sample_with_RLU) > 0) {

  message("\nThese samples actually have RLU values, meaning they were tested in house: ")
  message(paste0(missing_sample_with_RLU, collapse = ", "))
  message("These samples will NOT be excluded")

  point_of_care_samples <- point_of_care_samples[!grepl(paste0(missing_sample_with_RLU, collapse = "|"), point_of_care_samples)]

}

message("\nThese are low RLU samples less than 1000: ")
message(paste0(PHL_data[PHL_data$RLU < 1000 & !is.na(PHL_data$RLU), "SPECIMEN_NUMBER"], collapse = ", "))

message("\nThese are Health Center samples with missing RLU values: ")
message(paste0(point_of_care_samples[point_of_care_samples %in% PHL_data$SPECIMEN_NUMBER], collapse = ", "))

message("\nThese are MEO samples that we don't have: ")
message(paste0(MEO_samples[MEO_samples %in% PHL_data$SPECIMEN_NUMBER], collapse = ", "))

message("\nThese are GeneXpert samples with missing RLU values: ")
message(paste0(GX_samples[GX_samples %in% PHL_data$SPECIMEN_NUMBER], collapse = ", "))

message("\nThese are Respiratory Panel samples with missing RLU values: ")
message(paste0(Resp_samples[Resp_samples %in% PHL_data$SPECIMEN_NUMBER], collapse = ", "))

epi_sample_not_found <- PHL_data %>%
  filter(is.na(RLU)) %>%
  select(SPECIMEN_NUMBER) %>%
  filter(!SPECIMEN_NUMBER %in% c(point_of_care_samples, GX_samples, Resp_samples, MEO_samples)) %>%
  pull() %>%
  str_sort()

if(length(epi_sample_not_found) > 0) {

  stop(simpleError(paste("These samples were found in the epidemiologists metadata sheet but are missing RLU values: ",
                         paste0(epi_sample_not_found, collapse = ", "),
                         "Check these samples on Harvest", sep = "\n")))
}

#####################
# Load Temple samples
#####################

potential_ct_col_names <- c("ct value", "CT value", "CT values", "CTvalue", "CTvalues",
                            "ct values", "ctvalue", "ctvalues")

TU_data <- PHL_fp %>%
  lapply(function(x) read_excel_safely(x, sheet = "Temple") %>%
           cbind(filename = x) %>%
           as.data.frame(stringsAsFactors = FALSE)) %>%
  bind_rows()

if(ncol(TU_data) == 1) {
  TU_data <- data.frame(SPECIMEN_NUMBER = "") %>%
    mutate(`ct value` = "")
} else {
  TU_data <- TU_data %>%
    #filter rows where sample_id is NA
    filter(!is.na(SPECIMEN_NUMBER)) %>%
    #mutate(Collection_date = format(Collection_date, "%m/%d/%Y")) %>%
    rename_with(~ "ct value", any_of(potential_ct_col_names)) %>%
    as.data.frame()
}

################
# Filter samples
################

PHL_samples_removed <- PHL_data %>%
  filter(SPECIMEN_NUMBER != "") %>%
  filter(RLU < 1000 | SPECIMEN_NUMBER %in% c(MEO_samples, missing_samples)) %>%
  select(SPECIMEN_NUMBER) %>%
  pull()

filtered_PHL_data <- PHL_data %>%
  #remove low RLU samples
  filter(RLU >= 1000 | is.na(RLU)) %>%
  filter(!SPECIMEN_NUMBER %in% PHL_samples_removed)

TU_samples_removed <- TU_data %>%
  filter(SPECIMEN_NUMBER != "") %>%
  filter(`ct value` > 33 | SPECIMEN_NUMBER %in% c(missing_samples)) %>%
  select(SPECIMEN_NUMBER) %>%
  pull()

filtered_TU_data <- TU_data %>%
  #remove low RLU samples
  filter(`ct value` <= 33 | is.na(`ct value`)) %>%
  filter(!SPECIMEN_NUMBER %in% TU_samples_removed)

##############################################################################
# Write samples to Excel to send back to epidemiologist and wet lab scientists
##############################################################################

excel_data <- list(PHL = filtered_PHL_data, Temple = filtered_TU_data)

other_sheets <- lapply(PHL_fp, excel_sheets) %>%
  unlist() %>%
  unique()

other_sheets <- other_sheets[!grepl("PHL|Temple", other_sheets)]

other_samples <- data.frame(sample_name = "")

for(sheet_name in other_sheets) {

  possible_sample_names <- "SPECIMEN_NUMBER"

  other_data <- PHL_fp %>%
    lapply(function(x) read_excel_safely(x, sheet_name) %>%
             cbind(filename = x) %>%
             as.data.frame(stringsAsFactors = FALSE)) %>%
    bind_rows() %>%
    mutate_at(vars(contains(possible_sample_names)), ~gsub("\\s", "", .)) %>%
    as.data.frame()

  other_samples <- other_data %>%
    rename_at(vars(contains(possible_sample_names)),
              ~gsub(possible_sample_names, "sample_name", ., ignore.case = TRUE)) %>%
    arrange(across(starts_with("ct value"))) %>%
    arrange(across(starts_with("RLU"), desc)) %>%
    select(sample_name) %>%
    rbind(other_samples)

  other_data <- list(other_data)

  names(other_data) <- sheet_name

  excel_data <- append(excel_data, other_data)

}

for(file in PHL_fp) {

  filtered_fp <- gsub(".xlsx$", "_filtered.xlsx", file)

  filtered_excel <- excel_data %>%
    lapply(function(y) filter(y, if_any(matches("filename"), ~.x == file)) %>%
             select(-any_of("filename")))

  message(paste("Writing to"), filtered_fp)
  write.xlsx(filtered_excel, file = filtered_fp)
}

#########################################################################################
# Make a preliminary platemap for scientists with environmental samples and rerun samples
#########################################################################################

dir.create(here("metadata", "for_scientists"), showWarnings = FALSE)

empty_plate <- data.frame(plate_row = unlist(lapply(LETTERS[1:8], function(x) rep(x, 12))), plate_col = sprintf("%02d", rep(1:12, 8)), plate = 1) %>%
  mutate(plate_coord = paste0(plate, "_", plate_row, plate_col)) %>%
  arrange(plate_col) %>%
  mutate(sample_order = row_number()) %>%
  select(plate, plate_row, plate_col, plate_coord, sample_order)

PHL_samples <- filtered_PHL_data %>%
  rename(sample_name = "SPECIMEN_NUMBER") %>%
  arrange(desc(RLU))

TU_samples <- filtered_TU_data %>%
  rename(sample_name = "SPECIMEN_NUMBER") %>%
  arrange(`ct value`)

#if the shared drive can be accessed, copy the environmental swabs over
if(file.exists(shared_drive_fp)) {

  shared_environ_fp <- max(list.files(file.path(shared_drive_fp, "Sequencing_files", "2_Enviromental_samples", format(Sys.Date(), "%Y")),
                                      pattern = "^[0-9]*-[0-9]*-[0-9]*", full.names = TRUE, recursive = TRUE))

  environmental_file_date <- as.Date(gsub("_.*", "", basename(shared_environ_fp)))

  if(is.na(shared_environ_fp)) {
    shared_environ_fp <- list.files(file.path(shared_drive_fp, "Sequencing_files", "2_Enviromental_samples"),
                                    pattern = "^YYYY-MM-DD", full.names = TRUE, recursive = TRUE)

    environmental_file_date <- Sys.Date()
  }

  #if the date of the latest environmental samples is within 5 days of sequencing request, use this file
  if((Sys.Date() - environmental_file_date) < 5) {
    file.copy(shared_environ_fp, here("metadata", "extra_metadata"))
  }
}

environmental_samples_fp <- list.files(here("metadata", "extra_metadata"), pattern = "*_Environmental_Swab.*.xlsx", full.names = TRUE)

if(length(environmental_samples_fp) > 0) {
  enviro_samples <- environmental_samples_fp %>%
    data.frame(FileName = .) %>%
    group_by(FileName) %>%
    do(read_excel(.$FileName, col_names = TRUE)) %>%
    ungroup() %>%
    select(-FileName) %>%
    rename(sample_name = 1, environmental_site = 2) %>%
    select(sample_name, environmental_site)
} else {
  enviro_samples <- data.frame(sample_name = paste0("ENV", 1:enviro_number), environmental_site = paste0("ENV", 1:enviro_number))
  message("\nEnvironmental swab file was not found")
  Sys.sleep(5)
}

older_samples_fp <- list.files(here("metadata", "extra_metadata", "prev_run"), pattern = "_filtered.xlsx", full.names = TRUE)

if(length(older_samples_fp) > 0) {

  older_PHL_samples <- older_samples_fp %>%
    lapply(function(x) read_excel_safely(x, "PHL"))

  older_TU_samples <- older_samples_fp %>%
    lapply(function(x) read_excel_safely(x, "Temple"))

  older_samples <- bind_rows(older_PHL_samples, older_TU_samples) %>%
    rename(sample_name = "SPECIMEN_NUMBER") %>%
    mutate(PHL_sample = grepl("^H", sample_name)) %>%
    arrange(-PHL_sample, sample_name) %>%
    select(sample_name) %>%
    filter(sample_name != "") %>%
    filter(!is.na(sample_name))

} else {
  older_samples <- data.frame(sample_name = "")
}

older_samples_removed <- older_samples %>%
  filter(sample_name != "") %>%
  filter(sample_name %in% missing_samples) %>%
  pull()

older_samples <- older_samples %>%
  filter(!sample_name %in% older_samples_removed)

########################
# Report missing samples
########################

message("\nNumber of samples removed:")
message(length(PHL_samples_removed)+length(TU_samples_removed)+length(older_samples_removed))
message("\nNew samples removed:")
message(paste(c(PHL_samples_removed, TU_samples_removed), collapse = ", "))
message("\nRerun samples removed:")
message(paste(c(older_samples_removed), collapse = ", "))

#################
# Combine samples
#################

combined_list <- select(PHL_samples, sample_name) %>%
  rbind(select(TU_samples, sample_name)) %>%
  rbind(older_samples) %>%
  rbind(other_samples) %>%
  rbind(select(enviro_samples, sample_name)) %>%
  filter(sample_name != "") %>%
  #put samples in groups of 8
  mutate(grp = (row_number() - 1) %/% 8)

combined_list_first_half <- data.frame(sample_name = "NC-pre-extract", grp = 0) %>%
  rbind(combined_list) %>%
  filter(grp < 7)

combined_list_second_half <- data.frame(sample_name = "NC-pre-extract", grp = 7) %>%
  rbind(combined_list) %>%
  filter(grp >= 7)

if(nrow(combined_list_second_half) == 1) {
  combined_list_second_half <- data.frame(sample_name = "", grp = "")
}

plate_view <- combined_list_first_half %>%
  rbind(combined_list_second_half) %>%
  filter(sample_name != "") %>%
  group_by(grp) %>%
  group_modify(~ add_row(.x, sample_name = "NC-pre-extract")) %>%
  ungroup() %>%
  select(-grp) %>%
  mutate(number = cumsum(duplicated(sample_name)) + 1) %>%
  mutate(sample_name = ifelse(sample_name == "NC-pre-extract", paste0(sample_name, number), sample_name)) %>%
  select(-number) %>%
  rbind(data.frame(sample_name = c("BLANK", "PC", "NC-pre-cDNA", "NC-pre-ARTIC", "NC-pre-library"))) %>%
  mutate(sample_order = row_number()) %>%
  merge(empty_plate, by = "sample_order", all = TRUE) %>%
  mutate(sample_name = case_when(sample_order == 96 ~ "NC-corner",
                                 is.na(sample_name) ~ "",
                                 TRUE ~ sample_name)) %>%
  arrange(sample_order)

real_plate_view <- plate_view %>%
  select(sample_name, plate_row, plate_col) %>%
  tidyr::pivot_wider(names_from = "plate_col", values_from = "sample_name")

write_csv(plate_view, file = here("metadata", "for_scientists", paste0(format(Sys.time(), "%Y%m%d"), "_combined_samples_list.csv")))

plate_map_local_fp <- here("metadata", "for_scientists", paste0(format(Sys.time(), "%Y%m%d"), "_combined_samples_plate_map.csv"))
write_csv(real_plate_view, file = plate_map_local_fp)

if(copy_platemap) {

  if(file.exists(shared_drive_fp)) {

    plate_map_cp_fp <- file.path(shared_drive_fp, "Sequencing_files", "1_Plate_Maps", "nasal_swabs", format(Sys.Date(), "%Y"),
                                 paste(format(Sys.time(), "%Y-%m-%d"), sample_type_acronym, pathogen_acronym, "Plate_Map.csv", sep = "_"))
    dir.create(dirname(plate_map_cp_fp), showWarnings = FALSE)
    file.copy(plate_map_local_fp, plate_map_cp_fp, overwrite = TRUE)

  } else{

    message("\n*****")
    message("Could not access shared drive path. Plate map not copied")
    message("*****")
    Sys.sleep(5)

  }
}

write_csv(enviro_samples, file = here("metadata", "extra_metadata", paste0(format(Sys.time(), "%Y%m%d"), "_environmental_samples.csv")))
