library(here)
library(dplyr)
library(readxl)
library(readr)
library(stringr)
library(openxlsx)

#This Rscript filters out low RLU values from the metadata sheet received from the epidemiologists
#Send this filtered sheet to the epidemiologists and scientists

###################################################
# Load functions
###################################################

#this file needs to sit in a [aux_files/functions] directory path above this project directory
tryCatch(
  {
    source(file.path(dirname(here()), "aux_files", "functions", "R_all_functions_v3.R"))
  },
  error = function(e) {
    stop (simpleError("The R_all_functions_v3.R file needs to sit in a [aux_files/functions] directory path above this project directory"))
  }
)

###################################################
# Load the RLU report
# Make sure these sheets are not uploaded to GitHub
###################################################

#the HARVEST report is automatically generated each Monday and deposited onto a shared drive
RLU_file_name <- "COVID_harvest_report.csv"
date_RLU_file_name <- paste0(format(Sys.time(), "%Y%m%d"), "_", RLU_file_name)

shared_RLU_fp <- list.files("//city.phila.local/shares/Health/PHL/Admin/Sequencing_harvest_reports", pattern = paste0("^", date_RLU_file_name, "$"), full.names = TRUE)

file.copy(shared_RLU_fp, file.path(dirname(shared_RLU_fp), date_RLU_file_name))
file.copy(shared_RLU_fp, file.path(here("metadata", "extra_metadata"), date_RLU_file_name))

RLU_fp <- here("metadata", "extra_metadata", date_RLU_file_name)

RLU_data <- read_csv(RLU_fp) %>%
  filter(Test == "SARSCoV2-1") %>%
  rename(SPECIMEN_NUMBER = "Sample ID", SPECIMEN_DATE = "Draw Date", RLU = "Num Res", BIRTH_DATE = "DOB") %>%
  select(SPECIMEN_NUMBER, SPECIMEN_DATE, BIRTH_DATE, RLU) %>%
  mutate(BIRTH_DATE = as.Date(BIRTH_DATE, format = "%m/%d/%Y"), SPECIMEN_DATE = as.Date(SPECIMEN_DATE, format = "%m/%d/%Y")) %>%
  #filter rows where sample_id is NA
  filter(!is.na(SPECIMEN_NUMBER)) %>%
  #filter empty columns
  select(where(function(x) any(!is.na(x)))) %>%
  select(!matches("^\\.\\.\\."))

###################################################
# Get HC1 samples
###################################################

HC1_samples <- read_csv(RLU_fp) %>%
  filter(Test == "POCT 4 Plex Sars-CoV-2") %>%
  select(`Sample ID`) %>%
  pull()

###################################################
# Get GeneXpert samples
###################################################

GX_samples <- read_csv(RLU_fp) %>%
  filter(Test == "GeneXpert 4Plex Sars-CoV-2") %>%
  select(`Sample ID`) %>%
  pull()

###################################################
# Get Respiratory Panel samples
###################################################

Resp_samples <- read_csv(RLU_fp) %>%
  filter(Test == "Severe Acute Respiratory Syndrome Coronavirus 2 (SARS-CoV-2)") %>%
  select(`Sample ID`) %>%
  pull()

########################################################################
# Load the metadata sheet from epidemiologists and merge with RLU values
########################################################################

PHL_all_fp <- list.files(here("metadata", "extra_metadata"), pattern = "PHLspecimens.*.xlsx", full.names = TRUE)
PHL_fp <- PHL_all_fp[!grepl("_filtered.xlsx$", PHL_all_fp)]

PHL_data <- read_excel(PHL_fp, skip = 1, sheet = "PHL") %>%
  #filter rows where sample id is NA
  filter(!is.na(SPECIMEN_NUMBER)) %>%
  merge(RLU_data, by = c("SPECIMEN_NUMBER", "BIRTH_DATE", "SPECIMEN_DATE"), all.x = TRUE) %>%
  select(SPECIMEN_DATE, FacCode, agecoll, case_id, SPECIMEN_NUMBER,
         FIRST_NAME, LAST_NAME, BIRTH_DATE, age, zip_char, GENDER,
         breakthrough_case, death, hospitalized, outbreak, priority, RLU) %>%
  mutate(SPECIMEN_DATE = format(SPECIMEN_DATE, "%m/%d/%Y"), BIRTH_DATE = format(BIRTH_DATE, "%m/%d/%Y"))

if(any(is.na(PHL_data[PHL_data$SPECIMEN_NUMBER %in% RLU_data$SPECIMEN_NUMBER, "RLU"]))) {

  no_RLU <- PHL_data %>%
    filter(SPECIMEN_NUMBER %in% RLU_data$SPECIMEN_NUMBER) %>%
    filter(is.na(RLU)) %>%
    select(SPECIMEN_NUMBER) %>%
    pull()

  stop(simpleError(paste("Serious error! These samples have an RLU value but did not get added to PHL_data", no_RLU, sep = "\n")))

}

################################
# Remove Health Center 1 samples
################################

missing_sample_with_RLU <- PHL_data %>%
  filter(!is.na(RLU)) %>%
  filter(SPECIMEN_NUMBER %in% HC1_samples) %>%
  select(SPECIMEN_NUMBER) %>%
  pull()

if(length(missing_sample_with_RLU) > 0) {

  message("\nThese samples actually have RLU values, meaning they were tested in house: ")
  message(paste0(missing_sample_with_RLU, collapse = ", "))
  message("These samples will NOT be excluded")

  HC1_samples <- HC1_samples[!grepl(paste0(missing_sample_with_RLU, collapse = "|"), HC1_samples)]

}

message("\nThese are low RLU samples less than 1000: ")
message(paste0(PHL_data[PHL_data$RLU < 1000 & !is.na(PHL_data$RLU), "SPECIMEN_NUMBER"], collapse = ", "))

message("\nThese are HC1 samples with missing RLU values: ")
message(paste0(HC1_samples[HC1_samples %in% PHL_data$SPECIMEN_NUMBER], collapse = ", "))

message("\nThese are GeneXpert samples with missing RLU values: ")
message(paste0(GX_samples[GX_samples %in% PHL_data$SPECIMEN_NUMBER], collapse = ", "))

message("\nThese are Respiratory Panel samples with missing RLU values: ")
message(paste0(Resp_samples[Resp_samples %in% PHL_data$SPECIMEN_NUMBER], collapse = ", "))

epi_sample_not_found <- PHL_data %>%
  filter(is.na(RLU)) %>%
  select(SPECIMEN_NUMBER) %>%
  filter(!SPECIMEN_NUMBER %in% c(HC1_samples, GX_samples, Resp_samples)) %>%
  pull() %>%
  str_sort()

if(length(epi_sample_not_found) > 0) {

  stop(simpleError(paste("These samples were found in the epidemiologists metadata sheet but are missing RLU values: ",
                         paste0(epi_sample_not_found, collapse = ", "),
                         "Check these samples on Harvest", sep = "\n")))
}

samples_removed <- PHL_data %>%
  filter(RLU < 1000 | SPECIMEN_NUMBER %in% HC1_samples)

message("\nNumber of samples removed: ")
message(nrow(samples_removed))

##############################################################################
# Write samples to Excel to send back to epidemiologist and wet lab scientists
##############################################################################

filtered_PHL_data <- PHL_data %>%
  #remove low RLU samples
  filter(RLU >= 1000 | is.na(RLU)) %>%
  filter(!SPECIMEN_NUMBER %in% HC1_samples)

potential_ct_col_names <- c("CT value", "CT values", "CTvalue", "CTvalues",
                            "ct values", "ctvalue", "ctvalues")

TU_data <- read_excel(PHL_fp, sheet = "Temple") %>%
  #filter rows where sample_id is NA
  filter(!is.na(SPECIMEN_NUMBER)) %>%
  mutate(Collection_date = format(Collection_date, "%m/%d/%Y")) %>%
  rename_with(~ "ct value", any_of(potential_ct_col_names)) %>%
  as.data.frame()

excel_data <- list(PHL = filtered_PHL_data, Temple = TU_data)

other_sheets <- excel_sheets(PHL_fp)[!grepl("PHL|Temple", excel_sheets(PHL_fp))]

other_samples <- data.frame(sample_name = "")

for(sheet_name in other_sheets) {

  possible_sample_names <- "SPECIMEN_NUMBER"

  other_data <- PHL_fp %>%
    lapply(function(x) read_excel_safely(x, sheet_name)) %>%
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

write.xlsx(excel_data, file = gsub(".xlsx$", "_filtered.xlsx", PHL_fp))

################################################################################
# Make a preliminary platemap for scientists (no environmental samples included)
################################################################################

dir.create(here("metadata", "for_scientists"))

empty_plate <- data.frame(plate_row = unlist(lapply(LETTERS[1:8], function(x) rep(x, 12))), plate_col = sprintf("%02d", rep(1:12, 8)), plate = 1) %>%
  mutate(plate_coord = paste0(plate, "_", plate_row, plate_col)) %>%
  arrange(plate_col) %>%
  mutate(sample_order = row_number()) %>%
  select(plate, plate_row, plate_col, plate_coord, sample_order)

PHL_samples <- filtered_PHL_data %>%
  rename(sample_name = "SPECIMEN_NUMBER") %>%
  arrange(-RLU)

TU_samples <- TU_data %>%
  rename(sample_name = "SPECIMEN_NUMBER") %>%
  arrange(`ct value`)

environmental_samples_fp <- list.files(here("metadata", "extra_metadata"), pattern = "enviro.*.xlsx", full.names = TRUE, ignore.case = TRUE)

if(length(environmental_samples_fp) > 0) {
  enviro_samples <- environmental_samples_fp %>%
    data_frame(FileName = .) %>%
    group_by(FileName) %>%
    do(read_excel(.$FileName, col_names = TRUE)) %>%
    ungroup() %>%
    select(-FileName) %>%
    rename(sample_name = 1, environmental_site = 2) %>%
    select(sample_name, environmental_site)
} else {
  enviro_samples <- data.frame(sample_name = paste0("ENV", 1:8), environmental_site = paste0("ENV", 1:8))
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
    select(sample_name)

} else {
  older_samples <- data.frame(sample_name = "")
}

combined_list <- select(PHL_samples, sample_name) %>%
  rbind(select(TU_samples, sample_name)) %>%
  rbind(older_samples) %>%
  rbind(other_samples) %>%
  rbind(select(enviro_samples, sample_name)) %>%
  filter(sample_name != "") %>%
  #put samples in groups of 8
  mutate(grp = (row_number() - 1) %/% 8)

combined_list_first_half <- data.frame(sample_name = "NC", grp = 0) %>%
  rbind(combined_list) %>%
  filter(grp < 7)

combined_list_second_half <- data.frame(sample_name = "NC", grp = 7) %>%
  rbind(combined_list) %>%
  filter(grp >= 7)

if(nrow(combined_list_second_half) == 1) {
  combined_list_second_half <- data.frame(sample_name = "", grp = "")
}

plate_view <- combined_list_first_half %>%
  rbind(combined_list_second_half) %>%
  filter(sample_name != "") %>%
  group_by(grp) %>%
  group_modify(~ add_row(.x, sample_name = "NC")) %>%
  ungroup() %>%
  select(-grp) %>%
  mutate(number = cumsum(duplicated(sample_name)) + 1) %>%
  mutate(sample_name = ifelse(sample_name == "NC", paste0(sample_name, number), sample_name)) %>%
  select(-number) %>%
  rbind(data.frame(sample_name = c("BLANK", "PC"))) %>%
  mutate(sample_order = row_number()) %>%
  merge(empty_plate, by = "sample_order", all = TRUE) %>%
  mutate(sample_name = case_when(sample_order == 96 ~ "NC_corner",
                                 is.na(sample_name) ~ "",
                                 TRUE ~ sample_name))

real_plate_view <- plate_view %>%
  select(sample_name, plate_row, plate_col) %>%
  tidyr::pivot_wider(names_from = "plate_col", values_from = "sample_name")

write_csv(plate_view, file = here("metadata", "for_scientists", paste0(format(Sys.time(), "%Y%m%d"), "_combined_samples_list.csv")))

plate_map_local_fp <- here("metadata", "for_scientists", paste0(format(Sys.time(), "%Y%m%d"), "_combined_samples_plate_map.csv"))
write_csv(real_plate_view, file = plate_map_local_fp)

plate_map_cp_fp <- file.path("//city.phila.local/shares/Health/PHL/Admin/Sequencing Action plan updated/Plate Maps",
                             paste0(format(Sys.time(), "%Y-%m-%d"), "_Plate_Map.csv"))

file.copy(plate_map_local_fp, plate_map_cp_fp, overwrite = TRUE)

write_csv(enviro_samples, file = here("metadata", "extra_metadata", paste0(format(Sys.time(), "%Y%m%d"), "_environmental_samples.csv")))
