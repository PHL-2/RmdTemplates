library(here)
library(dplyr)
library(readxl)
library(readr)
library(stringr)
library(openxlsx)

#TODO: Include GeneXpert report and Repiratory Panel report and provide messages if found

#This Rscript filters out low RLU values from the metadata sheet received from the epidemiologists
#Send this filtered sheet to the scientists

###################################################
# Load the RLU report
# Make sure these sheets are not uploaded to GitHub
###################################################

RLU_fp <- list.files(here("metadata", "extra_metadata"), pattern = ".CSV", full.names = TRUE)

RLU_data <- read_csv(RLU_fp) %>%
  filter(Test == "SARSCoV2-1") %>%
  rename(SPECIMEN_NUMBER = "Sample ID", SPECIMEN_DATE = "Draw Date", RLU = "Num Res", BIRTH_DATE = "DOB") %>%
  #rename(SPECIMEN_NUMBER = "Sample ID", SPECIMEN_DATE = "Draw Date", RLU = "Num Res", GENDER = "Sex", age = "Age", BIRTH_DATE = "DOB") %>%
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

receiving_samples_not_found <- read_csv(RLU_fp) %>%
  filter(Test == "POCT 4 Plex Sars-CoV-2") %>%
  select(`Sample ID`) %>%
  pull()

########################################################################
# Load the metadata sheet from epidemiologists and merge with RLU values
########################################################################

PHL_all_fp <- list.files(here("metadata", "extra_metadata"), pattern = ".xlsx", full.names = TRUE)
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
  filter(SPECIMEN_NUMBER %in% receiving_samples_not_found) %>%
  select(SPECIMEN_NUMBER) %>%
  pull()

if(length(missing_sample_with_RLU) > 0) {

  message("\nThese samples actually have RLU values, meaning they were tested in house: ")
  message(paste0(missing_sample_with_RLU, collapse = ", "))
  message("These samples will NOT be excluded")

  receiving_samples_not_found <- receiving_samples_not_found[!grepl(paste0(missing_sample_with_RLU, collapse = "|"), receiving_samples_not_found)]

}

message("\nThese are HC1 samples with missing RLU values: ")
message(paste0(receiving_samples_not_found[receiving_samples_not_found %in% PHL_data$SPECIMEN_NUMBER], collapse = ", "))

epi_sample_not_found <- PHL_data %>%
  filter(is.na(RLU)) %>%
  select(SPECIMEN_NUMBER) %>%
  filter(!SPECIMEN_NUMBER %in% receiving_samples_not_found) %>%
  pull() %>%
  str_sort()

if(length(epi_sample_not_found) > 0) {

  stop(simpleError(paste("These samples were found in the epidemiologists metadata sheet but are missing RLU values: ",
                         paste0(epi_sample_not_found, collapse = ", "),
                         "Check the email to see if these samples could not be located by the receiving department and add them to receiving_samples_not_found",
                         "Otherwise, check these samples on Harvest. They may have been tested on the GeneXpert, which we will keep", sep = "\n")))

}

##############################################################################
# Write samples to Excel to send back to epidemiologist and wet lab scientists
##############################################################################

filtered_PHL_data <- PHL_data %>%
  #remove low RLU samples
  filter(RLU >= 1000 | is.na(RLU)) %>%
  filter(!SPECIMEN_NUMBER %in% receiving_samples_not_found)

TU_data <- read_excel(PHL_fp, sheet = "Temple") %>%
  #filter rows where sample_id is NA
  filter(!is.na(SPECIMEN_NUMBER)) %>%
  mutate(Collection_date = format(Collection_date, "%m/%d/%Y")) %>%
  as.data.frame()

excel_data <- list(PHL = filtered_PHL_data, Temple = TU_data)

write.xlsx(excel_data, file = gsub(".xlsx$", "_filtered.xlsx", PHL_fp))
