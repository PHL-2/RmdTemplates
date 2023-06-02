library(here)
library(dplyr)
library(readxl)
library(readr)
library(stringr)

#This Rscript is currently written to generating the SampleSheet for the Local Run Manager Module on the MiSeq
#https://support.illumina.com/downloads/local-run-manager-generate-fastq-module-v3.html
#BCL Convert may require a different SampleSheet format
#https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl_convert/bcl-convert-v3-7-5-software-guide-1000000163594-00.pdf

#TODO: detect any samples in the previous_run folder and change its run number to 2
#can also have an option here to add the run number to each sample
#TODO: try removing manual inputs into the sample metadata sheet

munging_fp <- here("metadata", "munge")

##############
# Manual input
##############

#s3 filepath to upload compressed sequencing run folder
s3_bucket_incoming_fp <- "s3://incoming-run"

prj_description <- "COVIDSeq" #no spaces, should be the same as the R project

instrument_select <- 1 #select 1 for MiSeq or 2 for NextSeq
instrument_type <- c("MiSeq", "NextSeq")[instrument_select]

read_length <- "76"
index_length <- "10"

phi_info <- c("sample_name", "zip_char", "case_id", "breakthrough_case", "death", "hospitalized", "outbreak", "priority")

#file location of the nextera udi indices
#don't have to change this if the file sits in a the metadata_references directory in the parent directory of the project
barcode_fp <- file.path(dirname(here()), "aux_files", "metadata_references", "nextera-dna-udi-samplesheet-MiSeq-flex-set-a-d-2x151-384-samples.csv")

#sequencing date will get grabbed from the R project name
sequencing_date <- gsub("_.*", "", basename(here())) #YYYY-MM-DD

if(sequencing_date == "" | prj_description == "" | any(is.na(instrument_type))) {
  stop (simpleError(paste0("Please fill in the sequencing date, short project description, or correct instrument in ", munging_fp, "/generate_barcodes_IDT.R")))
} else if (is.na(as.Date(sequencing_date, "%Y-%m-%d")) | nchar(sequencing_date) == 8) {
  stop (simpleError("Please enter the date into [sequencing_date] as YYYY-MM-DD"))
}

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
# Load config
###################################################

#this file needs to sit in a [aux_files/config] directory path above this project directory
tryCatch(
  {
    source(file.path(dirname(here()), "aux_files", "config", "config_variables.R"))
  },
  error = function(e) {
    stop (simpleError("The config_variables.R file needs to sit in a [aux_files/config] directory path above this project directory"))
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
    stop (simpleError("The nextera-dna-udi-samplesheet-MiSeq-flex-set-a-d-2x151-384-samples.csv file needs to sit in an [aux_files/metadata_references] directory path above this project directory"))
  }
)

if(sequencing_date == "2022-08-01") {
  barcodes <- data.frame(idt_plate_coord = paste0("Z_", LETTERS[1:8], "01"),
                         I7_Index_ID = "UDP9999", I5_Index_ID = "UDP9999", UDI_Index_ID = "UDP9999",
                         index = c("CTTCCTAGGA", "GAGGCCTATT", "GTGACACGCA", "CTGACTCTAC",
                                   "AGATCCATTA", "GTGGACAAGT", "ATTATCCACT", "AGTGTTGCAC"),
                         index2 = c("CCTAGAGTAT", "CTAGTCCGGA", "GCTTACGGAC", "ACGGCCGTCA",
                                    "ATCTCTACCA", "CCGTGGCCTT", "TACGCACGTA", "CTGGTACACG"))
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


###################################################################################
# Load the metadata sheet from epidemiologists and merge with sample metadata sheet
# Make sure these sheets are not uploaded to GitHub
###################################################################################

PHL_fp <- list.files(here("metadata", "extra_metadata"), pattern = "_filtered.xlsx", full.names = TRUE, recursive = TRUE)

PHL_data <- PHL_fp %>%
  lapply(function(x) read_excel_safely(x, "PHL")) %>%
  bind_rows()

#if PHL_data exists, do the following
if(ncol(PHL_data) > 0) {

  PHL_data <- PHL_data %>%
    mutate(SPECIMEN_DATE = as.Date(SPECIMEN_DATE, format = "%m/%d/%Y"), BIRTH_DATE = as.Date(BIRTH_DATE, format = "%m/%d/%Y")) %>%
    rename(sample_name = "SPECIMEN_NUMBER", sample_collection_date = "SPECIMEN_DATE", gender = "GENDER", DOB = "BIRTH_DATE") %>%
    select(any_of(phi_info), sample_collection_date, DOB, age, gender, RLU) %>%
    mutate(gender = ifelse(is.na(gender), "Unknown", gender)) %>%
    #filter rows where sample_id is NA
    filter(!is.na(sample_name)) %>%
    #make these columns character vectors
    mutate(across(c(sample_name, case_id, breakthrough_case, priority, gender), as.character))
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
  bind_rows()

if(ncol(TU_data) > 0) {

  TU_data <- TU_data %>%
    mutate(Collection_date = as.Date(Collection_date, format = "%m/%d/%Y")) %>%
    rename(sample_name = "SPECIMEN_NUMBER", sample_collection_date = "Collection_date", CT = "ct value", gender = "GENDER") %>%
    select(any_of(phi_info), sample_collection_date, CT, age, gender) %>%
    #filter rows where sample_id is NA
    filter(!is.na(sample_name)) %>%
    #make these columns character vectors
    mutate(across(c(sample_name, case_id, breakthrough_case, priority, gender), as.character))
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

ENV_fp <- list.files(here("metadata", "extra_metadata"), pattern = "environmental_samples.csv", full.names = TRUE)

if(length(ENV_fp) > 0) {

  ENV_data <- read_csv(ENV_fp) %>%
    #use the first day of the week (starting on Monday) as the sample_collection_date
    mutate(sample_collection_date = as.Date(cut(as.POSIXct(Sys.time()), "week"))) %>%
    select(sample_name, sample_collection_date, environmental_site) %>%
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
                         zip_char = NA_integer_
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
  mutate(sample_id = gsub("_", "-", paste0("PHL2", "-", idt_plate_coord, "-", gsub("-", "", sequencing_date)))) %>%
  select(sample_id, everything()) %>%
  arrange(plate, plate_col, plate_row) %>%
  mutate(sequencing_date = sequencing_date) %>%
  mutate(prj_descrip = prj_description) %>%
  mutate(instrument_type = instrument_type) %>%
  mutate(read_length = read_length) %>%
  mutate(index_length = index_length) %>%
  mutate(environmental_material = ifelse(grepl("swab|control", sample_type, ignore.case = TRUE), NA, environmental_material)) %>%
  mutate(collection_device = ifelse(grepl("waste water|control", sample_type, ignore.case = TRUE), NA, collection_device)) %>%
  mutate(environmental_site = case_when(grepl("Water control|Reagent control", sample_type) ~ paste0(plate_row, plate_col),
                                        grepl("Environmental control", sample_type) ~ paste0(environmental_site, " - ", plate_row, plate_col),
                                        TRUE ~ environmental_site))

if(min(as.Date(metadata_sheet$sample_collection_date[!is.na(metadata_sheet$sample_collection_date)])) < seq(as.Date(sequencing_date), length=2, by='-2 month')[2]){
  stop(simpleError(paste0("Some samples have collection dates more than 2 months ago. Investigate!!")))
}

######################################################################################
# Samples in epi metadata but we don't have the samples or they could not be extracted
######################################################################################

message("\nThese samples were found in the epidemiologists metadata sheet but not in our sample sheet. Check the email to see if these samples could not be located by the receiving department")
message("Otherwise, these samples may have had an issue during extraction. Send wet lab scientists these sample names to check")

epi_sample_not_found <- PHL_data %>%
  select(any_of("sample_name")) %>%
  rbind(select(TU_data, any_of("sample_name")))

if(ncol(epi_sample_not_found > 0)) {
  epi_sample_not_found <- epi_sample_not_found %>%
    filter(!sample_name %in% metadata_sheet$sample_name) %>%
    pull() %>%
    str_sort()
}

message(paste0(epi_sample_not_found, collapse = ", "))

missing_metadata_non_ctrl_samples <- metadata_sheet %>%
  filter(!grepl("control", sample_type)) %>%
  filter(sample_type != "Wastewater")

missing_sample_date <- missing_metadata_non_ctrl_samples %>%
  filter(is.na(sample_collection_date)) %>%
  select(sample_name)

missing_sample_RLU <- missing_metadata_non_ctrl_samples %>%
  filter(sample_collected_by == "Philadelphia Department of Public Health") %>%
  filter(is.na(RLU)) %>%
  select(sample_name)

missing_sample_CT <- missing_metadata_non_ctrl_samples %>%
  filter(sample_collected_by == "Temple University") %>%
  filter(is.na(CT)) %>%
  select(sample_name)

missing_metadata_samples <- rbind(missing_sample_date, missing_sample_RLU, missing_sample_CT)

if(nrow(missing_metadata_samples) > 0){
  stop(simpleError(paste0("These non-control samples are in the sample sheet but are missing RLU values, CT values, or collection date from the epidemiologists!\n",
                          "They may also be GeneXpert samples\n",
                          paste0(pull(missing_metadata_samples), collapse = ", "))))
}

#############
# Check sheet
#############

#throw error if missing these columns
for(x in c(cols2merge, "sample_id",
           "idt_set", "idt_plate_row", "idt_plate_col", "idt_plate_coord",
           "sample_type", "sample_collected_by", "PHL_sample_received_date",
           "organism", "host_scientific_name", "host_disease", "isolation_source",
           "index", "index2", "UDI_Index_ID", "I7_Index_ID", "I5_Index_ID",
           "sequencing_date", "prj_descrip", "instrument_type", "read_length", "index_length",
           "environmental_site", "sample_collection_date", "gender", "zip_char")) {
  if(!grepl(paste0(colnames(metadata_sheet), collapse = "|"), x)) {
    stop(simpleError(paste0("Missing column [", x, "] in the metadata sheet template!!!")))
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
if(min(as.Date(metadata_sheet$sample_collection_date[!is.na(metadata_sheet$sample_collection_date)])) < seq(as.Date(sequencing_date), length=2, by='-2 month')[2]){
  stop(simpleError(paste0("Some samples have collection dates more than 2 months ago. Investigate!!")))
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
  # do not include lane in the sample sheet otherwise it will only demultiplex that sample in that specified lane, not in all lanes
  rowwise() %>%
  #BCL Convert does not take Index Plate
  select(sample_id, index, index2) %>%
  rename(Sample_ID = "sample_id")

sample_sheet_fn <- paste0(sequencing_date, "_SampleSheet_v2.csv")

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
write_samp("")

write_samp("[BCLConvert_Data]")
write_csv(samp_sheet_2_write, file = sample_sheet_fp, col_names = TRUE, append = TRUE)

################################
# Write sheet to metadata folder
################################

#does not contain PHI and accession numbers
metadata_sheet %>%
  select(-c(I7_Index_ID, I5_Index_ID, any_of(phi_info))) %>%
  write.csv(file = here("metadata", paste0(sequencing_date, "_", prj_description, "_metadata.csv")), row.names = FALSE)

#contains PHI and accession numbers
metadata_sheet %>%
  select(sample_id, any_of(phi_info)) %>%
  write.csv(file = here("metadata", paste0(sequencing_date, "_", prj_description, "_PHI.csv")), row.names = FALSE)

###########################
#Get sequencing folder path
###########################

#get the newly added run folder
run_folder <- sequencing_folder_fp %>%
  list.files(full.names = T) %>%
  file.info() %>%
  filter(grepl(format(as.Date(sequencing_date), "%y%m%d"), rownames(.))) %>%
  rownames()

folder_date <- paste0("20", gsub(".*/|_.*", "", run_folder)) %>%
  as.Date(format = "%Y%m%d")

if(folder_date != gsub("_.*", "", basename(here()))) {
  stop(simpleError("The run date on the sequencing folder does not match the date of this RStudio project!"))
}

sequencing_run <- gsub(".*/", "", run_folder)

s3_samplesheet_fp <- file.path(s3_bucket_incoming_fp, sequencing_date, sample_sheet_fn)

nf_demux_samplesheet <- data.frame(
  id = sequencing_run,
  samplesheet = s3_samplesheet_fp,
  lane = "",
  flowcell = file.path(s3_bucket_incoming_fp, sequencing_date, paste0(sequencing_run, ".tar.gz"))
)

nf_demux_samplesheet_fp <- here("metadata", "munge", paste0(sequencing_date, "_nf_demux_samplesheet.csv"))
s3_nf_demux_samplesheet_fp <- file.path(s3_bucket_incoming_fp, sequencing_date, paste0(sequencing_date, "_nf_demux_samplesheet.csv"))

nf_demux_samplesheet %>%
  write.csv(file = nf_demux_samplesheet_fp,
            row.names = FALSE, quote = FALSE)

# aws_cli_cp_fp <- file.path(dirname(here()), "aux_files", "aws_cli_commands", "aws_cli_1_cp.sh")
#
# aws_cp_samplesheet <- cli_submit(sh_exe_fp, aws_cli_cp_fp,
#                                  c(shQuote(sample_sheet_fp, type = "sh"),
#                                    shQuote(s3_samplesheet_fp, type = "sh")))
#
# aws_cp_nf_demux_samplesheet <- cli_submit(sh_exe_fp, aws_cli_cp_fp,
#                                           c(shQuote(nf_demux_samplesheet_fp, type = "sh"),
#                                             shQuote(s3_nf_demux_samplesheet_fp, type = "sh")))
#
# if(!all(grepl("Completed", c(aws_cp_samplesheet, aws_cp_nf_demux_samplesheet), ignore.case = TRUE))) {
#   stop(simpleError("Samplesheet upload to s3 bucket failed"))
# }
