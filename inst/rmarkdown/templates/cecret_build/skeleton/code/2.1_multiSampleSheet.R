###############
# Combining sample sheets
###############


#### Trigger
metadata_input_fp <- list.files(path = file.path(here(), "inputs"),
                                pattern = "sequencing_metadata_sheet",
                                full.names = TRUE)

if(length(checkNumMetadataSheets) > 1) {
  source(here(), "code/2.1_multiSampleSheet.R")
}


#### Process

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



samp_sheet_2_write <- metadata_sheet %>%
  filter(!sample_name %in% remove_sample_from_bcl_samplesheet) %>%
  # do not include lane in the sample sheet otherwise it will only demultiplex that sample in that specified lane, not in all lanes
  rowwise() %>%
  #BCL Convert does not take Index Plate
  select(sample_id, index, index2) %>%
  rename(Sample_ID = "sample_id")

sample_sheet_fn <- paste(sequencing_date, instrument_type, sample_type_acronym, prj_description, "SampleSheet_v2.csv", sep = "_")

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
#don't split the lanes
write_samp(c("NoLaneSplitting", "true"))
write_samp(c("FastqCompressionFormat", "gzip"))
write_samp("")

write_samp("[BCLConvert_Data]")
write_csv(samp_sheet_2_write, file = sample_sheet_fp, col_names = TRUE, append = TRUE)

################################
# Write sheet to metadata folder
################################

merged_samples_metadata_sheet <- metadata_sheet %>%
  select(-c(sample_id, index, index2,
            starts_with("idt_"), starts_with("plate"), ends_with("_ID", ignore.case = FALSE))) %>%
  filter(!sample_type %in% c(sequencing_controls, "Environmental control")) %>%
  group_by(uniq_sample_name) %>%
  mutate(sample_counts = n()) %>%
  ungroup() %>%
  filter(sample_counts > 1) %>%
  select(-c(sample_name, sample_counts)) %>%
  unique() %>%
  mutate(index = "TTTTTTTTTT",
         index2 = "TTTTTTTTTT",
         UDI_Index_ID = "UDP9999",
         idt_set = "Merged",
         idt_plate_row = rep(LETTERS[1:8], 12)[row_number()],
         idt_plate_col = unlist(lapply(1:12, function (x) rep(x, 8)))[row_number()],
         idt_plate_coord = paste0(idt_set, "_", idt_plate_row, idt_plate_col),
         plate = 999,
         plate_row = idt_plate_row,
         plate_col = idt_plate_col,
         plate_coord = paste0(plate, "_", plate_row, plate_col),
         across(matches("_col$|coord$"), ~ str_replace_all(., "\\d+", function(m) sprintf("%02d", as.numeric(m)))),
         sample_id = gsub("_", "-", paste0("PHL2", "-", instrument_regex, "-", idt_plate_coord, "-", gsub("-", "", sequencing_date))),
         sample_name = paste0(uniq_sample_name, "-Merged"),
         ww_group = paste(ww_group, "- Merged"))

#does not contain PHI and accession numbers
metadata_sheet %>%
  select(-c(I7_Index_ID, I5_Index_ID)) %>%
  # NextSeq runs have index2 sequences in the reverse complement of the ones listed in the reference barcode sheet
  # however, the sample sheet used for demultiplexing needs to be the same orientation as the reference barcode sheet
  # when BCLConvert demultiplexes a NextSeq run, it will automatically reverse complement index2
  # results from the pipeline will refer to the reverse complement of index2, so this should be updated in the metadata sheet
  mutate(index2 = ifelse(instrument_type == "NextSeq1k2k", reverse_complement(index2), index2)) %>%
  rbind(merged_samples_metadata_sheet) %>%
  write_csv(file = here("metadata", paste0(sequencing_date, "_", prj_description, "_metadata.csv")))

bclconvert_output_final_path <- paste(s3_fastq_bucket, sequencing_date, sample_type_acronym, prj_description, "processed_bclconvert", sequencing_run, sep = "/")

#write the sample sheet for merging nextflow script
metadata_sheet %>%
  select(fastq = sample_id, uniq_sample_name) %>%
  mutate(fastq_1 = paste0(bclconvert_output_final_path, "/", fastq, "_S", row_number(), "_R1_001.fastq.gz"),
         fastq_2 = paste0(bclconvert_output_final_path, "/", fastq, "_S", row_number(), "_R2_001.fastq.gz")) %>%
  merge(select(merged_samples_metadata_sheet, sample_id, uniq_sample_name), by = "uniq_sample_name", all = TRUE) %>%
  filter(!is.na(sample_id)) %>%
  select(sample_id, fastq_1, fastq_2) %>%
  write_csv(file = here("metadata", "munge",
                        tolower(paste(sequencing_date, instrument_type, sample_type_acronym, prj_description, "nf_concat_fastq_samplesheet.csv", sep = "_"))))

