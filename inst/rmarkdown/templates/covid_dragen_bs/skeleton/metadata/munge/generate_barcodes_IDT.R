library(here)
library(dplyr)
library(readxl)

#This Rscript is currently written to generating the SampleSheet for the Local Run Manager Module on the MiSeq
#https://support.illumina.com/downloads/local-run-manager-generate-fastq-module-v3.html
#BCL Convert may require a different SampleSheet format
#https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl_convert/bcl-convert-v3-7-5-software-guide-1000000163594-00.pdf

munging_fp <- here("metadata", "munge")

##############
# Manual input
##############

sequencing_date <- "2022-06-30" #YYYY-MM-DD
prj_description <- "CovidSeq" #no spaces

if(sequencing_date == "" | prj_description == "") {
  stop (simpleError(paste0("Please fill in the sequencing date or the short project description in ", munging_fp, "/generate_barcodes_IDT.R")))
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
  select(idt_plate_coord, I7_Index_ID, index, index2)

#####################
# Load metadata sheet
#####################

metadata_input_fp <- here(munging_fp, list.files(munging_fp, pattern = ".xlsx"))

read_sheet <- function(fp, sheet_name) {
  read_excel(fp, sheet = sheet_name) %>%
    #filter rows where sample_id is NA
    filter(!is.na(sample_id)) %>%
    #filter empty columns
    select(where(function(x) any(!is.na(x))))
}

qubit_sheet <- read_sheet(metadata_input_fp, "Qubit")
index_sheet <- read_sheet(metadata_input_fp, "Index")
library_sheet <- read_sheet(metadata_input_fp, "Library")
extra_sheet <- read_sheet(metadata_input_fp, "Extra")

cols2merge <- c("sample_id", "plate", "plate_row", "plate_col", "plate_coord")

metadata_sheet <- merge(qubit_sheet, index_sheet, by = cols2merge, all = TRUE, sort = FALSE) %>%
  merge(library_sheet, by = cols2merge, all = TRUE, sort = FALSE) %>%
  merge(extra_sheet, by = cols2merge, all = TRUE, sort = FALSE) %>%
  mutate(sequencing_date = sequencing_date)

#throw error if missing these columns
for(x in c(cols2merge,
         "qubit_conc_ng_ul", "qubit_date",
         "idt_set", "idt_plate_row", "idt_plate_col", "idt_plate_coord",
         "library_conc_ng_ul",
         "sample_type", "lane", "run_number")) {
  if(!grepl(paste0(colnames(metadata_sheet), collapse = "|"), x)) {
    stop(simpleError(paste0("Missing column [", x, "]!!! in the metadata sheet!")))
  }
}

####################
# Write sample sheet
####################


#print(temp)

#temp <- filter(temp, pooled)

print("Which section?")
print(all_lanes[i])
print('Are the barcode columns unique?')
print(length(unique(paste0(temp$barcode_index_set, temp$barcode_coord))) == dim(temp)[1])
print('Number of barcodes?')
print(length(unique(paste0(temp$barcode_index_set, temp$barcode_coord))))
print('Are the barcodes unique?')
print(length(unique(paste0(temp$I7_index_seq, temp$I5_index_seq))) == dim(temp)[1])
print('Are the sample names unique?')
print(length(unique(temp$SampleID)) == dim(temp)[1])
print('Are all the forward primers found?')
print(sum(is.na(sample_sheet$temp$I7_index_seq)) == 0)
print('Are all the reverse primers found?')
print(sum(is.na(sample_sheet$temp$I5_index_seq)) == 0)
print('Do all the sampleIDs start with a letter?')
print(all(grepl("^[A-Za-z]", temp$SampleID)))
print('Are there illegal characters in the SampleID?')
print(any(grepl(" |_|-", temp$SampleID)))

merger sample sheet first with barcode

have a project destiption, date, index, index2, sampleid

################################
# Write sheet to metadata folder
################################

write.csv(metadata_sheet, file = here("metadata", paste0(sequencing_date, "_", prj_description, ".csv")), row.names = FALSE)



