library(here)
library(dplyr)
library(readxl)
library(readr)
library(stringr)
library(openxlsx)

#This Rscript filters out low RLU values from the metadata sheet received from the epidemiologists
#Send this filtered sheet to the scientists

##############
# Manual input
##############

# add the sample names that could not be retrieved from our receiving department as a string; this list from an email by Morris replying to Jasmine
receiving_samples_not_found <- '' %>%
  str_split(pattern = ", ") %>%
  unlist()

###################################################
# Load the RLU report
# Make sure these sheets are not uploaded to GitHub
###################################################

RLU_fp <- list.files(here("metadata", "extra_metadata"), pattern = ".csv", full.names = TRUE)

RLU_data <- read_csv(RLU_fp) %>%
  filter(`SARS-CoV2 Result` == "POSITIVE") %>%
  rename(SPECIMEN_NUMBER = "Sample ID", SPECIMEN_DATE = "Draw Date", RLU = "SARSCoV2-1", GENDER = "Sex", age = "Age", BIRTH_DATE = "DOB") %>%
  select(SPECIMEN_NUMBER, SPECIMEN_DATE, BIRTH_DATE, age, GENDER, RLU) %>%
  mutate(GENDER = case_when(GENDER == "M" ~ "Male",
                            GENDER == "F" ~ "Female",
                            TRUE ~ NA_character_)) %>%
  mutate(BIRTH_DATE = as.Date(BIRTH_DATE, format = "%m/%d/%Y"), SPECIMEN_DATE = as.Date(SPECIMEN_DATE, format = "%m/%d/%Y")) %>%
  mutate(age = ifelse(grepl("mo", age), 0, age)) %>%
  #filter rows where sample_id is NA
  filter(!is.na(SPECIMEN_NUMBER)) %>%
  #filter empty columns
  select(where(function(x) any(!is.na(x)))) %>%
  select(!matches("^\\.\\.\\."))

########################################################################
# Load the metadata sheet from epidemiologists and merge with RLU values
########################################################################

PHL_all_fp <- list.files(here("metadata", "extra_metadata"), pattern = ".xlsx", full.names = TRUE)
PHL_fp <- PHL_all_fp[!grepl("_filtered.xlsx$", PHL_all_fp)]

PHL_data <- read_excel(PHL_fp, skip = 1, sheet = "PHL") %>%
  #filter rows where sample id is NA
  filter(!is.na(SPECIMEN_NUMBER)) %>%
  merge(RLU_data, by = c("SPECIMEN_NUMBER", "BIRTH_DATE", "age", "GENDER", "SPECIMEN_DATE"), all.x = TRUE) %>%
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
message(paste0(receiving_samples_not_found, collapse = ", "))

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
