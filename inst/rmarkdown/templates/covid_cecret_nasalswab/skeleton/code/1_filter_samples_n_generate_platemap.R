library(here)
library(dplyr)
library(readxl)
library(readr)
library(stringr)
library(openxlsx)

#This Rscript appends COVID test results to the metadata sheet received from the epidemiologists
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

###############
# Load R config
###############

#this file needs to sit in a [aux_files/r_scripts/config] directory path above this project directory
tryCatch(
  {
    source(file.path(dirname(here()), "aux_files", "r_scripts", "config", "config_variables.R"))
  },
  error = function(e) {
    stop (simpleError("The config_variables.R file needs to sit in a [aux_files/r_scripts/config] directory path above this project directory"))
  }
)

############################
# Load COVID report from OEL
############################

# Look for this report in extra_metadata folder
covid_test_fn <- "COVID_SEQ_[0-9]*.csv"
covid_test_fp <- list.files(here("metadata", "extra_metadata"), pattern = paste0("^", covid_test_fn, "$"), full.names = TRUE)

# If the file does not exist, grab it from the shared drive
# The report should be automatically generated each week and deposited into the shared drive
if(length(covid_test_fp) == 0) {

  #path of latest covid report
  shared_covid_test_fp <- max(list.files(file.path(shared_drive_fp, "Sequencing_OEL_reports"), pattern = paste0("^", covid_test_fn, "$"), full.names = TRUE))

  #get date of report
  report_date <- as.Date(str_extract(shared_covid_test_fp, "[0-9]{6}"), tryFormats = c("%y%m%d"))

  if((Sys.Date() - report_date) > 5) {
    message(paste("\nThe most recent COVID report found on the shared drive is from", report_date))
  }

  covid_test_fp <- here("metadata", "extra_metadata", basename(shared_covid_test_fp))

  file.copy(shared_covid_test_fp, covid_test_fp)

}

covid_test_data <- read_csv(covid_test_fp) %>%
  select_all(~tolower(gsub(" ", "_", .))) %>%
  select(SPECIMEN_NUMBER = "sample_id", SPECIMEN_DATE = "collection_date", BIRTH_DATE = "patient_dob_(mm/dd/yyyy)", # columns used to merge
         ordering_location, clinical_test_name = "test_name", clinical_test_result = "result") %>%
  mutate(SPECIMEN_DATE = format(as.Date(SPECIMEN_DATE, format = "%m/%d/%Y"), "%m/%d/%Y"),
         BIRTH_DATE = format(as.Date(BIRTH_DATE, format = "%m/%d/%Y"), "%m/%d/%Y"),
         clinical_test_result = case_when(grepl("^positive$|^detected$", clinical_test_result, ignore.case = TRUE) ~ "positive",
                                          grepl("^negative$|^not detected$|^presumptive negative$", clinical_test_result, ignore.case = TRUE) ~ "negative",
                                          grepl("^inconclusive|^invalid|^no result|^error", clinical_test_result, ignore.case = TRUE) ~ "inconclusive",
                                          grepl("^cancel", clinical_test_result, ignore.case = TRUE) ~ "canceled",
                                          TRUE ~ NA)) %>%
  #filter rows where sample_id is NA
  filter(!is.na(SPECIMEN_NUMBER))

if(any(is.na(covid_test_data$clinical_test_result))) {

  na_result <- covid_test_data %>%
    filter(is.na(clinical_test_result)) %>%
    select(clinical_test_result) %>%
    unique() %>%
    pull()

  stop(simpleError(paste("The OEL report has a new value in the result column that has not been previously captured",
                         "Please add this new value to the template:",
                         paste(na_result, collapse = ", "), sep = "\n")))
}

test_name_msg <- setNames(c("POCT 4 Plex Sars-CoV-2",                                       # HC POC
                            "GeneXpert 4Plex Sars-CoV-2",                                   # Cepheid
                            "SARS-CoV2 Result",                                             # Panther (no RLU data until Harvest is merged with OEL)
                            "Severe Acute Respiratory Syndrome Coronavirus 2 (SARS-CoV-2)", # Respiratory panel
                            "SARS CoV 2 RNA, RT PCR",                                       # Quest
                            "MEO POC Test Result",                                          # MEO
                            "COVID-19"),                                                    # Old test
                          c("POC tested samples from the Health Centers",
                            "GeneXpert samples",
                            "Panther samples (non-numeric data)",
                            "Respiratory panel samples",
                            "Quest samples",
                            "MEO samples (not received)",
                            "Old COVID test samples"))

if(!all(covid_test_data$clinical_test_name %in% test_name_msg)) {

  new_test <- covid_test_data %>%
    filter(!clinical_test_name %in% test_name_msg) %>%
    select(clinical_test_name) %>%
    unique() %>%
    pull()

  stop(simpleError(paste("The OEL report has a new test that has not been accounted for",
                         "Please add this new test to the template:",
                         paste(new_test, collapse = ", "), sep = "\n")))
}

##########################################################################
# Load the metadata sheet from epidemiologists and merge with COVID report
##########################################################################

PHL_all_fp <- list.files(here("metadata", "extra_metadata"), pattern = "PHLspecimens.*.xlsx", full.names = TRUE)
PHL_fp <- PHL_all_fp[!grepl("_filtered.xlsx$", PHL_all_fp)]

PHL_data <- PHL_fp %>%
  lapply(function(x) read_excel_safely(x, sheet = "PHL", skip_row = 1) %>%
           cbind(filename = x) %>%
           as.data.frame(stringsAsFactors = FALSE)) %>%
  bind_rows()

if(ncol(PHL_data) == 1) {

  PHL_data <- data.frame(SPECIMEN_NUMBER = "", RLU = "")

  show_samps <- FALSE

} else {

  cols2merge <- c("SPECIMEN_NUMBER", "SPECIMEN_DATE", "BIRTH_DATE")

  PHL_data <- PHL_data %>%
    #filter rows where sample id is NA
    filter(!is.na(SPECIMEN_NUMBER)) %>%
    mutate(RLU = NA, # when Harvest data is available in the report, merge the numerical data into RLU
           #convert to date and then to string to remove any ambiguity
           SPECIMEN_DATE = format(as.Date(SPECIMEN_DATE, format = "%m/%d/%Y"), "%m/%d/%Y"),
           BIRTH_DATE = format(as.Date(BIRTH_DATE, format = "%m/%d/%Y"), "%m/%d/%Y")) %>%
    left_join(covid_test_data, by = cols2merge)

  #PHL_data checks
  non_positive_PHL_results <- PHL_data %>%
    filter(clinical_test_result != "positive" | is.na(clinical_test_result) | is.na(clinical_test_name))

  if(nrow(non_positive_PHL_results) > 0) {

    stop(simpleError(paste("Some samples do not have a positive result or any result",
                           "It's possible that the epi's metadata are missing the following columns so the clinical results could not be merged:",
                           paste0(cols2merge, collapse = ", "),
                           "\nThe latest date found in the clinical test metadata sheet from OEL is:", max(as.Date(covid_test_data$SPECIMEN_DATE, format = "%m/%d/%Y"), na.rm = TRUE),
                           paste0("\nIf the test results are not missing, check with ", epi_name, " to see if these samples should still be sequenced"),
                           "If you want to proceed with sequencing these samples without the clinical data, comment out this stop message",
                           "\nSample(s) in question:", paste0(non_positive_PHL_results[, "SPECIMEN_NUMBER"], collapse = ", "), sep = "\n")))
  }

  if(any(c(is.na(PHL_data$GENDER) | PHL_data$GENDER == "" | is.na(PHL_data$age) | PHL_data$age == ""), na.rm = TRUE)) {

    missing_meta <- PHL_data %>%
      filter(is.na(GENDER) | GENDER == "" | is.na(age) | age == "") %>%
      select(SPECIMEN_NUMBER) %>%
      pull()

    stop(simpleError(paste("Some samples are missing patient age or gender\n",
                           "Please fill in the missing information\n",
                           "Sample(s) in question:\n",
                           paste0(missing_meta, collapse = ", "))))
  }

  show_samps <- TRUE
}

#############################
# Where are the samples from?
#############################

if(show_samps) {

  message("\nWhere are the requested samples from?")

  for(test in test_name_msg) {
    message(paste0("\n", names(test_name_msg)[test_name_msg == test], ":"))

    samples2print <- PHL_data %>%
      filter(clinical_test_name == test) %>%
      select(SPECIMEN_NUMBER, SPECIMEN_DATE, ordering_location, clinical_test_result, RLU) %>%
      mutate(low_RLU_sample_removed = case_when(is.na(RLU) ~ "FALSE", #non-Panther samples should not be removed
                                                RLU < 1000 ~ "TRUE",
                                                TRUE ~ "FALSE"))

    if(nrow(samples2print) > 0) {
      samples2print %>%
        print.AsIs()
    }
  }

  message("")

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
    mutate(Collection_date = format(Collection_date, "%m/%d/%Y")) %>%
    rename_with(~ "ct value", any_of(potential_ct_col_names)) %>%
    as.data.frame()
}

################
# Filter samples
################

meo_samples <- covid_test_data %>%
  filter(clinical_test_result == "MEO POC Test Result") %>%
  select(SPECIMEN_NUMBER) %>%
  pull()

PHL_samples_2_remove <- PHL_data %>%
  filter(SPECIMEN_NUMBER != "") %>%
  #remove low RLU samples, samples in the missing_samples vector, and MEO samples
  filter(RLU < 1000 | SPECIMEN_NUMBER %in% c(missing_samples, meo_samples)) %>%
  select(SPECIMEN_NUMBER) %>%
  pull()

filtered_PHL_data <- PHL_data %>%
  filter(!SPECIMEN_NUMBER %in% PHL_samples_2_remove)

TU_samples_2_remove <- TU_data %>%
  filter(SPECIMEN_NUMBER != "") %>%
  #remove high CT samples
  filter(`ct value` > 33 | SPECIMEN_NUMBER %in% missing_samples) %>%
  select(SPECIMEN_NUMBER) %>%
  pull()

filtered_TU_data <- TU_data %>%
  filter(!SPECIMEN_NUMBER %in% TU_samples_2_remove)

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
  arrange(desc(RLU), sample_name)

TU_samples <- filtered_TU_data %>%
  rename(sample_name = "SPECIMEN_NUMBER") %>%
  arrange(`ct value`, sample_name)

#if the shared drive can be accessed, copy the environmental swabs over
if(file.exists(shared_drive_fp)) {

  shared_environ_fp <- suppressWarnings(max(list.files(file.path(shared_drive_fp, "Sequencing_files", "2_Enviromental_samples", format(Sys.Date(), "%Y")),
                                                       pattern = "^[0-9]*-[0-9]*-[0-9]*", full.names = TRUE, recursive = TRUE)))

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

older_samples_2_remove <- older_samples %>%
  filter(sample_name != "") %>%
  filter(sample_name %in% missing_samples) %>%
  pull()

older_samples <- older_samples %>%
  filter(!sample_name %in% older_samples_2_remove)

########################
# Report missing samples
########################

message("\nNumber of samples removed:")
message(length(PHL_samples_2_remove)+length(TU_samples_2_remove)+length(older_samples_2_remove))
message("\nNew samples removed:")
message(paste(c(PHL_samples_2_remove, TU_samples_2_remove), collapse = ", "))
message("\nRerun samples removed:")
message(paste(c(older_samples_2_remove), collapse = ", "))

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
  # removing last water control after samples due to request
  head(-1) %>%
  select(-grp) %>%
  mutate(number = cumsum(duplicated(sample_name)) + 1) %>%
  mutate(sample_name = ifelse(sample_name == "NC-pre-extract", paste0(sample_name, number), sample_name)) %>%
  select(-number) %>%
  rbind(data.frame(sample_name = c("PC", "NC-pre-cDNA", "NC-pre-ARTIC", "NC-pre-library"))) %>%
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

message("\nRscript finished successfully!")
