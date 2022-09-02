library(here)
library(dplyr)
library(readxl)
library(readr)
library(stringr)

#This Rscript is currently written to generating the SampleSheet for the Local Run Manager Module on the MiSeq
#https://support.illumina.com/downloads/local-run-manager-generate-fastq-module-v3.html
#BCL Convert may require a different SampleSheet format
#https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl_convert/bcl-convert-v3-7-5-software-guide-1000000163594-00.pdf

munging_fp <- here("metadata", "munge")

##############
# Manual input
##############

prj_description <- "COVIDSeq" #no spaces, should be the same as the R project

instrument_select <- 1 #select 1 for MiSeq or 2 for NextSeq
instrument_type <- c("MiSeq", "NextSeq")[instrument_select]

read_length <- "76"

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

#####################
# Load metadata sheet
#####################

metadata_input_fp <- list.files(here("metadata", "munge"), pattern = ".xlsx", full.names = TRUE)

read_sheet <- function(fp, sheet_name) {
  read_excel(fp, sheet = sheet_name) %>%
    #filter rows where sample_id is NA
    filter(!is.na(sample_name)) %>%
    #filter empty columns
    select(where(function(x) any(!is.na(x)))) %>%
    select(!matches("^\\.\\.\\.")) %>%
    mutate(across(matches("_col$|coord$"), ~ str_replace_all(., "\\d+", function(m) sprintf("%02d", as.numeric(m))))) %>%
    mutate(plate_coord = gsub("^0", "", plate_coord))
}

index_sheet <- read_sheet(metadata_input_fp, "Index")
sample_info_sheet <- read_sheet(metadata_input_fp, "Sample Info")

###################################################################################
# Load the metadata sheet from epidemiologists and merge with sample metadata sheet
# Make sure these sheets are not uploaded to GitHub
###################################################################################

PHL_fp <- list.files(here("metadata", "extra_metadata"), pattern = ".xlsx", full.names = TRUE)

PHL_data <- read_excel(PHL_fp, skip = 1) %>%
  rename(sample_name = "SPECIMEN_NUMBER", sample_collection_date = "SPECIMEN_DATE", gender = "GENDER") %>%
  select(sample_name, sample_collection_date, age, gender, zip_char, priority) %>%
  #filter rows where sample_id is NA
  filter(!is.na(sample_name)) %>%
  #filter empty columns
  select(where(function(x) any(!is.na(x)))) %>%
  select(!matches("^\\.\\.\\.")) %>%
  #use the first day of the week (starting on Monday) as the sample_collection_date
  mutate(sample_collection_date = as.Date(cut(as.POSIXct(sample_collection_date), "week"))) %>%
  mutate(host_age_bin = cut(age, breaks = c(0, 9, as.numeric(paste0(1:6, 9)), Inf),
                            labels = c("0 - 9", paste(seq(10, 60, by = 10), "-",as.numeric(paste0(1:6, 9))), "70+"),
                            include.lowest = TRUE)) %>%
  #don't include age because it may be PHI if included with zipcode and gender
  select(-age)

###################################################
# Load the monthly RLU report
# Make sure these sheets are not uploaded to GitHub
###################################################

RLU_fp <- list.files(here("metadata", "extra_metadata"), pattern = ".csv", full.names = TRUE)

RLU_data <- read_csv(RLU_fp) %>%
  rename(sample_name = "Sample ID") %>%
  select(sample_name, RLU) %>%
  #filter rows where sample_id is NA
  filter(!is.na(sample_name)) %>%
  #filter empty columns
  select(where(function(x) any(!is.na(x)))) %>%
  select(!matches("^\\.\\.\\."))

###########################
# Merge all metadata sheets
###########################

cols2merge <- c("sample_name", "plate", "plate_row", "plate_col", "plate_coord")

#merge all the individual sheets
metadata_sheet <- merge(index_sheet, sample_info_sheet, by = cols2merge, all = TRUE, sort = FALSE) %>%
  #add in barcodes
  merge(barcodes, by = "idt_plate_coord", all.x = TRUE, sort = FALSE) %>%
  #merge the metadata from epi's
  merge(PHL_data, by = "sample_name", all.x = TRUE, sort = FALSE) %>%
  #merge the RLU data
  merge(RLU_data, by = "sample_name", all.x = TRUE, sort = FALSE) %>%
  mutate(sample_id = gsub("_", "-", paste0("PHL2", "-", idt_plate_coord, "-", gsub("-", "", sequencing_date)))) %>%
  select(sample_id, everything()) %>%
  arrange(plate, plate_col, plate_row) %>%
  mutate(sequencing_date = sequencing_date) %>%
  mutate(prj_descrip = prj_description) %>%
  mutate(instrument_type = instrument_type) %>%
  mutate(read_length = read_length) %>%
  mutate(environmental_material = ifelse(grepl("swab|control", sample_type), NA, environmental_material)) %>%
  mutate(collection_device = ifelse(grepl("waste water|control", sample_type), NA, collection_device)) %>%
  #remove empty columns again
  select(where(function(x) any(!is.na(x))))

#############
# Check sheet
#############

#throw error if missing these columns
for(x in c(cols2merge, "sample_id",
           "idt_set", "idt_plate_row", "idt_plate_col", "idt_plate_coord",
           "sample_type", "sample_collected_by", "PHL_sample_received_date",
           "organism", "host_scientific_name", "host_disease", "isolation_source",
           "index", "index2", "UDI_Index_ID", "I7_Index_ID", "I5_Index_ID",
           "sequencing_date", "prj_descrip", "instrument_type", "read_length",
           "sample_collection_date", "host_age_bin", "gender", "zip_char", "priority")) {
  if(!grepl(paste0(colnames(metadata_sheet), collapse = "|"), x)) {
    stop(simpleError(paste0("Missing column [", x, "] in the metadata sheet template!!!")))
  }
}

#make these columns NA if missing
for(x in c("qubit_conc_ng_ul",
           "library_conc_ng_ul",
           "lane", "run_number")) {
  if(!grepl(paste0(colnames(metadata_sheet), collapse = "|"), x)) {
    metadata_sheet[[x]] <- NA
  }
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
  mutate(Index_Plate = which(LETTERS == idt_set)) %>%
  mutate(Index_Plate_Well = paste0(idt_plate_row, idt_plate_col)) %>%
  select(sample_id, Index_Plate, Index_Plate_Well, I7_Index_ID, index, I5_Index_ID, index2, UDI_Index_ID) %>%
  rename(Sample_ID = "sample_id")

sample_sheet_fp <- here("metadata", "munge", "SampleSheet.csv")

write_samp <- function(line2write) {
  write(paste0(line2write, collapse = ","), file = sample_sheet_fp, append = TRUE)
}

write("[Header]", file = sample_sheet_fp)
write_samp(c("Experiment Name", paste0(prj_description, "_", sequencing_date)))
write_samp(c("Date", sequencing_date))
write_samp(c("Workflow", "GenerateFASTQ"))
write_samp(c("Library Prep Kit", "COVIDSeq for Surveillance"))
write_samp(c("Index Kit", "COVIDSeq indexes_IDT for Illumina-PCR Indexes Set 1 2 3 4"))
write_samp(c("Chemistry", "Amplicon"))
write_samp(c("Instrument type", instrument_type))
write_samp("")

write_samp("[Reads]")
#writing read length twice for paired reads
write_samp(read_length)
write_samp(read_length)
write_samp("")

write_samp("[Settings]")
write_samp(c("adapter", "CTGTCTCTTATACACATCT"))
write_samp("")

write_samp("[Data]")
write_csv(samp_sheet_2_write, file = sample_sheet_fp, col_names = TRUE, append = TRUE)

################################
# Write sheet to metadata folder
################################

metadata_sheet %>%
  select(-c(I7_Index_ID, I5_Index_ID)) %>%
  write.csv(file = here("metadata", paste0(sequencing_date, "_", prj_description, ".csv")), row.names = FALSE)
