library(here)
library(dplyr)
library(readxl)
library(readr)

#This Rscript is currently written to generating the SampleSheet for the Local Run Manager Module on the MiSeq
#https://support.illumina.com/downloads/local-run-manager-generate-fastq-module-v3.html
#BCL Convert may require a different SampleSheet format
#https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl_convert/bcl-convert-v3-7-5-software-guide-1000000163594-00.pdf

munging_fp <- here("metadata", "munge")

##############
# Manual input
##############

sequencing_date <- "" #YYYY-MM-DD
prj_description <- "" #no spaces

instrument_select <- NA #select 1 for MiSeq or 2 for NextSeq
instrument_type <- c("MiSeq", "NextSeq")[instrument_select]

read_length <- "76"

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

#local index sequences
barcode_fp <- here(munging_fp, "nextera-dna-udi-samplesheet-MiSeq-flex-set-a-d-2x151-384-samples.csv")

barcodes <- read.csv(barcode_fp, stringsAsFactors = FALSE) %>%
  #change all ending 0 to another character
  mutate(Index_Plate_Well = gsub("0$", "zzz", Index_Plate_Well)) %>%
  #remove the middle 0 in the plate positions
  mutate(Index_Plate_Well = gsub("0", "", Index_Plate_Well)) %>%
  mutate(Index_Plate_Well = gsub("zzz", "0", Index_Plate_Well)) %>%
  mutate(idt_plate_coord = paste0(Index_Plate, "_", Index_Plate_Well)) %>%
  rename(UDI_Index_ID = "I7_Index_ID") %>%
  select(idt_plate_coord, UDI_Index_ID, index, index2)

#################################################
# Load metadata sheet and merge with barcode file
#################################################

metadata_input_fp <- here(munging_fp, list.files(munging_fp, pattern = ".xlsx"))

read_sheet <- function(fp, sheet_name) {
  read_excel(fp, sheet = sheet_name) %>%
    #filter rows where sample_id is NA
    filter(!is.na(sample_id)) %>%
    #filter empty columns
    select(where(function(x) any(!is.na(x)))) %>%
    select(!matches("^\\.\\.\\."))
}

qubit_sheet <- read_sheet(metadata_input_fp, "Qubit") %>%
  #mutate this column as character if exists
  mutate_at(vars(one_of('qubit_date')), as.character)
index_sheet <- read_sheet(metadata_input_fp, "Index")
library_sheet <- read_sheet(metadata_input_fp, "Library")
extra_sheet <- read_sheet(metadata_input_fp, "Extra")

cols2merge <- c("sample_id", "plate", "plate_row", "plate_col", "plate_coord")

metadata_sheet <- merge(qubit_sheet, index_sheet, by = cols2merge, all = TRUE, sort = FALSE) %>%
  merge(library_sheet, by = cols2merge, all = TRUE, sort = FALSE) %>%
  merge(extra_sheet, by = cols2merge, all = TRUE, sort = FALSE) %>%
  #add in barcodes
  merge(barcodes, by = "idt_plate_coord", all.x = TRUE) %>%
  mutate(sequencing_date = sequencing_date) %>%
  mutate(prj_descrip = prj_description) %>%
  mutate(instrument_type = instrument_type) %>%
  mutate(read_length = read_length)

#############
# Check sheet
#############

#throw error if missing these columns
for(x in c(cols2merge,
           "idt_set", "idt_plate_row", "idt_plate_col", "idt_plate_coord",
           "sample_type",
           "index", "index2",
           "sequencing_date", "prj_descrip")) {
  if(!grepl(paste0(colnames(metadata_sheet), collapse = "|"), x)) {
    stop(simpleError(paste0("Missing column [", x, "]!!! in the metadata sheet!")))
  }
}

#make these columns NA if missing
for(x in c("qubit_conc_ng_ul", "qubit_date",
           "library_conc_ng_ul",
           "lane", "run_number")) {
  if(!grepl(paste0(colnames(metadata_sheet), collapse = "|"), x)) {
    metadata_sheet[[x]] <- NA
  }
}

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
  stop(simpleError("There are spaces or periods in the Sample IDs! Please fix"))
}

####################
# Write sample sheet
####################

samp_sheet_2_write <- metadata_sheet %>%
  # do not include lane in the sample sheet otherwise it will only demultiplex that sample in that specified lane, not in all lanes
  select(sample_id, index, index2, UDI_Index_ID) %>%
  rename(Sample_ID = "sample_id")

sample_sheet_fp <- here("metadata", "munge", "SampleSheet.csv")

write_samp <- function(line2write) {
  write(paste0(line2write, collapse = ","), file = sample_sheet_fp, append = TRUE)
}

write("[Header]", file = sample_sheet_fp)
write_samp(c("Experiment Name", prj_description))
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

write.csv(metadata_sheet, file = here("metadata", paste0(sequencing_date, "_", prj_description, ".csv")), row.names = FALSE)
