#### Created by Christopher Gu, PHD, MPH 1/31/25 ####
####### Script for reading in standard ddPCR data exported from the QX Manager ##########

d_pre <- read_csv(copy_path, id = "File")

if(sum(grepl(pattern = "Sample description", x = colnames(d_pre)))>0) {
  d_pre <- d_pre %>%
    select(
      Well, Sample = `Sample description 1`, Spiked = `Sample description 2`, Sample_date = `Sample description 3`, Extra = `Sample description 4`, SampleType, Target, `Conc(copies/uL)` = `Conc(copies/µL)`, Droplet_count = `Accepted Droplets`, Positives, Negatives) %>%
    #mutate(Sample = gsub(pattern = "dd_", replacement = "dd", x = Sample)) %>%
    separate(Sample, into = c("Type", "Process"), sep = "_", remove = FALSE) %>%
    filter(!is.na(Sample),
           Sample != "Buffer") %>%
    mutate(posDrops10k = Positives/Droplet_count * 10000 * 100,
           Process = ifelse(is.na(Process), "Non-processed", Process),
           Process = ifelse(grepl(pattern = "dil", x = Extra, ignore.case = T), "Dilution", Process),
           Process = ifelse(Spiked == "Unspiked", "Unspiked", Process),
           Process = droplevels(factor(Process, level = c("Concentrated", "Unconcentrated", "Non-processed", "Unspiked","Dilution", "Frozen", "QuickExtract"))),
           `Conc(copies/uL)` = as.numeric(`Conc(copies/uL)`),
           Lysed = as.factor(ifelse(grepl(pattern = "lyse", x = Extra, ignore.case = T), "Lysed", "Non-Lysed")),
           Sample_date = as.Date(Sample_date, tryFormats = c("%m/%d/%Y", "%Y-%m-%d")),
           SampleType = ifelse(SampleType == "NTC", "NegCtrl", SampleType),
           Target_present = Positives >= 3
    )
} else {
  d_pre <- d_pre %>%
    separate(Sample, into = c("Sample", "Spiked", "Sample_date", "Extra"), sep = "(?<=[A-Za-z])-(?!H2O|SC2)|(?<!ddPCR|PBS)-(?=[A-za-z])") %>%
    pivot_longer(
      cols = c("N2 Cp/uL","E Cp/uL","MHV Cp/uL","N2 Negatives","E Negatives","MHV Negatives","N2 Positives","E Positives","MHV Positives",),
      names_to = c("Target", "secondary"),
      names_sep = " ",
      values_to = "value"
    ) %>%
    mutate(SampleType = ifelse(grepl(pattern = "SC2", x = Sample), "PosCtrl",
                               ifelse(Sample == "PBS" | Sample == "oldWW" | Sample == "ddPCR-H2O", "NegCtrl", "Unknown"))) %>%
    select(Well, Sample, Spiked, SampleType, Sample_date, Extra, Target, secondary, value, Droplets) %>%
    pivot_wider(
      id_cols = c("Well", "Sample", "Spiked", "SampleType","Sample_date", "Extra", "Target", "Droplets"),
      names_from = "secondary",
      values_from = "value"
    ) %>%
    rename(Droplet_count = "Droplets", `Conc(copies/uL)` = "Cp/uL") %>%
    #mutate(Sample = gsub(pattern = "dd_", replacement = "dd", x = Sample)) %>%
    separate(Sample, into = c("Type", "Process"), sep = "_", remove = FALSE) %>%
    filter(!is.na(Sample),
           Sample != "Buffer") %>%
    mutate(percentage_positive = Positives/Droplet_count * 100,
           Process = ifelse(is.na(Process), "Non-processed", Process),
           Process = ifelse(grepl(pattern = "dil", x = Extra, ignore.case = T), "Dilution", Process),
           Process = ifelse(Spiked == "Unspiked", "Unspiked", Process),
           Spiked = ifelse(grepl(pattern = "Spiked", x = Type), "Unspiked", Spiked),
           Process = droplevels(factor(Process, level = c("Concentrated", "Unconcentrated", "Non-processed", "Unspiked","Dilution", "Frozen", "QuickExtract"))),
           `Conc(copies/uL)` = as.numeric(`Conc(copies/uL)`),
           Lysed = as.factor(ifelse(grepl(pattern = "lyse", x = Extra, ignore.case = T), "Lysed", "Non-Lysed")),
           Sample_date = as.Date(Sample_date, tryFormats = c("%m/%d/%Y", "%Y-%m-%d")),
           SampleType = ifelse(SampleType == "NTC", "NegCtrl", SampleType),
           Target_present = Positives >= 3
    )
}
