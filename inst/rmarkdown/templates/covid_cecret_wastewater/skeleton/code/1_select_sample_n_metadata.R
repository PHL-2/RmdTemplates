library(here)
library(dplyr)
library(tidyverse)
library(readxl)
library(readr)
library(stringr)

#This Rscript adds the relevant metadata fields to wastewater samples for sequencing

##############
# Manual input
##############

create_sample_replicates <- 4 #enter a positive integer for the number of biological replicates created during concentration and extraction

# set this to TRUE to copy the platemap to the shared drive
copy_platemap <- TRUE

sample_type_acronym <- "WW" #use WW for wastewater samples

prj_description <- "COVIDSeq" #no spaces, should be the same as the R project

# number of unspecified environmental swabs to add to plate
#enviro_number <- 10

################
# Load functions
################

#this file needs to sit in a [aux_files/functions] directory path above this project directory
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

#this file needs to sit in a [aux_files/config] directory path above this project directory
tryCatch(
  {
    source(file.path(dirname(here()), "aux_files", "r_scripts", "config", "config_variables.R"))
  },
  error = function(e) {
    stop (simpleError("The config_variables.R file needs to sit in a [aux_files/r_scripts/config] directory path above this project directory"))
  }
)

#####################################
# Load the latest 5 ddPCR run results
#####################################

failed_regex <- "test|exclude"

ddPCR_files <- list.files(ddPCR_run_fp, pattern = ".*_ww_sequencing_metadata.csv", full.names = TRUE, recursive = TRUE)
ddPCR_files <- tail(ddPCR_files[!grepl(failed_regex, ddPCR_files)], 10)

ddPCR_data <- ddPCR_files %>%
  data_frame(FileName = .) %>%
  group_by(FileName) %>%
  do(read_delim(.$FileName,
                show_col_types = FALSE,
                col_types = cols("sample_received_date" = col_date()))) %>%
  ungroup() %>%
  mutate(ddpcr_analysis_date = as.Date(gsub(paste0(ddPCR_run_fp, "/|_.*"), "", FileName))) %>%
  group_by(sample_group, sample_received_date) %>%
  #get the latest run only
  filter(ddpcr_analysis_date == max(ddpcr_analysis_date)) %>%
  ungroup() %>%
  select(-FileName)

##########################
# Load the selection sheet
##########################

project_name <- basename(here())
sequencing_date <- as.Date(gsub("_.*", "", project_name))

select_fp <- list.files(here("metadata", "extra_metadata"), pattern = "*.csv", full.names = TRUE)
select_fp <- select_fp[!grepl("_filtered.csv$", select_fp)]

if(length(select_fp) == 0) {
  stop (simpleError("\n\nThere is no sample selection sheet in extra_metadata!"))
}

selection_data <- lapply(select_fp, read_csv) %>%
  do.call(rbind, .) %>%
  mutate(sample_received_date = as.Date(sample_received_date, tryFormats = c("%Y-%m-%d", "%m/%d/%y", "%m/%d/%Y"))) %>%
  select(any_of("sample_group"), sample_received_date) %>%
  #filter empty columns
  select(where(function(x) any(!is.na(x))),
         !matches("^\\.\\.\\.")) %>%
  merge(ddPCR_data, by = c("sample_group", "sample_received_date"), all.x = TRUE, sort = FALSE) %>%
  mutate(uniq_sample_name = ifelse((is.na(uniq_sample_name) | uniq_sample_name == ""),
                              paste("WW", sample_received_date, sample_group, sep = "-"),
                              uniq_sample_name))

if(all(is.na(selection_data$ddpcr_analysis_date))) {
  message("\n\nWarning!!!\nSamples selected for sequencing do not have a corresponding ddPCR date!")
  Sys.sleep(10)
}

###################################################
# Write selected samples and their metadata to file
###################################################

selection_data %>%
  write_csv(here("metadata", "extra_metadata", paste0(basename(here()), "_filtered.csv")))

#########################################################################################
# Make a preliminary platemap for scientists with environmental samples and rerun samples
#########################################################################################

dir.create(here("metadata", "for_scientists"))

empty_plate <- data.frame(plate_row = unlist(lapply(LETTERS[1:8], function(x) rep(x, 12))), plate_col = sprintf("%02d", rep(1:12, 8)), plate = 1) %>%
  mutate(plate_coord = paste0(plate, "_", plate_row, plate_col)) %>%
  arrange(plate_col) %>%
  mutate(sample_order = row_number()) %>%
  select(plate, plate_row, plate_col, plate_coord, sample_order)

#if the shared drive can be accessed, copy the environmental swabs over
# if(file.exists(shared_drive_fp)) {
#
#   shared_environ_fp <- max(list.files(file.path(shared_drive_fp, "Sequencing Action plan updated", "Enviromental_samples"),
#                                       pattern = "^[0-9]*-[0-9]*-[0-9]*", full.names = TRUE, recursive = TRUE))
#
#   environmental_file_date <- as.Date(gsub("_.*", "", basename(shared_environ_fp)))
#
#   #if the date of the latest environmental samples is within 5 days of sequencing request, use this file
#   if((Sys.Date() - environmental_file_date) < 5) {
#     file.copy(shared_environ_fp, here("metadata", "extra_metadata"))
#   }
# }

# environmental_samples_fp <- list.files(here("metadata", "extra_metadata"), pattern = "*_Environmental_Swab.*.xlsx", full.names = TRUE)
#
#
# if(length(environmental_samples_fp) > 0) {
#   enviro_samples <- environmental_samples_fp %>%
#     data.frame(FileName = .) %>%
#     group_by(FileName) %>%
#     do(read_excel(.$FileName, col_names = TRUE)) %>%
#     ungroup() %>%
#     select(-FileName) %>%
#     rename(any_of(c(uniq_sample_name = "sample_name",
#                     uniq_sample_name = 1,
#                     environmental_site = 2))) %>%
#     select(uniq_sample_name, environmental_site)
# } else {
#   enviro_samples <- data.frame(uniq_sample_name = paste0("ENV", 1:enviro_number), environmental_site = paste0("ENV", 1:enviro_number))
#   message("\nEnvironmental swab file was not found")
#   Sys.sleep(5)
# }

older_samples_fp <- list.files(here("metadata", "extra_metadata", "prev_run"), pattern = "_filtered.csv", full.names = TRUE)

if(length(older_samples_fp) > 0) {

  older_ddPCR_samples <- older_samples_fp %>%
    lapply(read_csv) %>%
    select(sample_group, sample_received_date, uniq_sample_name) %>%
    filter(sample_group != "",
           !is.na(sample_group))

} else {
  older_samples <- data.frame(sample_group = "", sample_received_date = "", uniq_sample_name = "")
}

#################
# Combine samples
#################

sample_group_order <- c("ZeptoSC2", "SouthWest", "NorthEast", "SouthEast", "oldWW")

combined_list <- select(selection_data, sample_group, sample_received_date, uniq_sample_name) %>%
  rbind(older_samples) %>%
  filter(sample_group != "",
         !is.na(sample_received_date)) %>%
  mutate(sample_group = factor(sample_group, levels = sample_group_order)) %>%
  arrange(sample_received_date, sample_group) %>%
  mutate(order = 1:nrow(.)) %>%
  expand(nesting(uniq_sample_name, order), rep = paste0("-Rep", 1:create_sample_replicates)) %>%
  mutate(sample_name = paste0(uniq_sample_name, rep)) %>%
  arrange(order) %>%
  select(sample_name) %>%
  #put samples in groups of 8
  mutate(grp = (row_number() - 1) %/% 8) %>%
  add_row(., sample_name = paste0("NC-pre-extract", 1:8), .before = ceiling(median(.$grp))*8+1) %>%
  filter(sample_name != "") %>%
  select(-grp)

plate_view <- combined_list %>%
  rbind(data.frame(sample_name = c("NC-pre-extract9", "BLANK", "PC", "NC-pre-cDNA", "NC-pre-ARTIC", "NC-pre-library"))) %>%
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

if(copy_platemap) {

  if(file.exists(shared_drive_fp)) {

    plate_map_cp_fp <- file.path(shared_drive_fp, "Sequencing Action plan updated", "1_Plate_Maps", "wastewater", format(Sys.Date(), "%Y"),
                                 paste(format(Sys.time(), "%Y-%m-%d"), sample_type_acronym, prj_description, "Plate_Map.csv", sep = "_"))
    file.copy(plate_map_local_fp, plate_map_cp_fp, overwrite = TRUE)

  } else{

    message("\n*****")
    message("Could not access shared drive path. Plate map not copied")
    message("*****")
    Sys.sleep(5)

  }
}

#write_csv(enviro_samples, file = here("metadata", "extra_metadata", paste0(format(Sys.time(), "%Y%m%d"), "_environmental_samples.csv")))

############################
# Report samples to sequence
############################

message("\nNumber of samples to sequence:")
message(nrow(selection_data))
message("\nDate of samples to sequence:")
message(paste0(unique(selection_data$sample_received_date), collapse = "\n"))
message("\nSample sites to sequence:")
message(paste0(unique(selection_data$sample_group), collapse = "\n"))
message("\nNumber of rerun samples to sequence:")
if(length(older_samples_fp) > 0) {
  message(nrow(older_samples))
} else {
  message(0)
}
