#### Created by Christopher Gu, PHD, MPH 1/31/25 ####
####### Script for cluster ddPCR amplitude data ##########

library(here)
library(tidyverse)
library(Ckmeans.1d.dp)

#Droplet volume used for calculations of concentration
droplet_vol <- 0.795 #nanoliter
#another R analysis script used 0.91 (https://link.springer.com/article/10.1007/s00216-015-8773-4)
# this R package uses 0.85 https://github.com/Gromgorgel/ddPCR
#QX Manager User Guide has the above volume. Can use that as the reference and put a link into the validation



#File path for sample info/meta
copyPathSampleInfo <-
  list.files(
    path = file.path(here(), "input"),
    pattern = "import",
    full.names = TRUE,
    recursive = FALSE
  )

#Loading amplitude data
droplets <- read_csv(copy_path,
           id = "FileName",
           col_names = T,
           skip = 3) %>%
  mutate(plate_coord = gsub(".*_", "", gsub("_Amplitude.csv", "", FileName))) %>%
  select(plate_coord,
         where(function(x)
           any(!is.na(x))),-FileName) %>%
  set_names(c("plate_coord", colnames(.)[(ceiling(ncol(.) / 2) + 1):ncol(.)], 2:(ceiling(ncol(.) / 2)))) %>%
  select(1:ceiling(ncol(.) / 2)) %>%
  pivot_longer(cols = -plate_coord,
               names_to = "channel",
               values_to = "amplitude") %>%
  mutate(
    amplitude = as.numeric(amplitude),
    plate_row = gsub("([A-H])[0-9]*", "\\1", plate_coord),
    plate_col = gsub("[A-H]", "", plate_coord)
  )

#Loading Sample info/metadata
sampleInfo <- read_csv(copyPathSampleInfo, col_names = T, skip = 4) %>%
  select(
    Well,
    Sample_name = `Sample description 1`,
    Spiked = `Sample description 2`,
    Sample_received = `Sample description 3`,
    Extra = `Sample description 4`,
    SampleType,
    channel = TargetName
  ) %>%
  filter(!is.na(Sample_name)) %>%
  separate(Sample_name,
           into = c("Type", "Process"),
           sep = "_")

#Perform clustering on positive control samples
thresholdDrops <- sampleInfo %>%
  filter(SampleType == "PosCtrl") %>%
  left_join(droplets, by = c("Well" = "plate_coord", "channel")) %>%
  group_by(Well,
           Type,
           Process,
           Spiked,
           Extra,
           Sample_received,
           SampleType,
           channel) %>%
  mutate(cluster = unlist(Ckmeans.1d.dp(amplitude, k = c(1, 2))["cluster"])) %>%
  group_by(Well,
           Type,
           Process,
           Spiked,
           Extra,
           Sample_received,
           SampleType,
           channel,
           cluster) %>%
  summarize(
    nAmp = n(),
    meanAmp = mean(amplitude),
    sdAmp = sd(amplitude)
  ) %>%
  mutate(
    clusterAssignment = case_when(
      meanAmp == min(meanAmp) ~ "Negatives",
      meanAmp == max(meanAmp) ~ "Positives",
      .default = NA
    )
  )

#Use positive control sample cluster amplitudes to assign clusters to rest of the data
clusterAssignment <- droplets %>%
  group_by_all() %>%
  left_join(thresholdDrops, by = c("channel")) %>%
  mutate(ampDiffPos = abs(amplitude - meanAmp)) %>%
  filter(ampDiffPos == min(ampDiffPos)) %>%
  select(-ampDiffPos,-meanAmp,-sdAmp)



#Generating work data for later analysis
d_pre <- clusterAssignment %>%
  pivot_wider(id_cols = c("Well", "type", "process", "Spiked", "Extra", "Sample_received", "channel", "SampleType"), names_from = clusterAssignment, values_from = nAmp) %>%
  rename(Target = "channel") %>%
  group_by(Well, Target) %>%
  mutate(Positives = ifelse(is.na(Positives), 0, Positives),
         Droplet_count = Negatives + Positives,
         lambda = -log(Negatives/Droplet_count), #average number of targets based on Poisson distribution
         `Conc(copies/uL)` = lambda*(1000/droplet_vol), #1000 to convert nL to uL and to account for droplet volume
         #percentage_positive = Positives/Droplet_count,
         Target_present = Positives >= 3) %>%
  select(-lambda)

# clusterAssignment %>%
#   #left_join(clusterAssignment, by = c("Well", "type", "process", "Extra", "Sample_received", "channel", "SampleType", "cluster")) %>%
#   #group_by(Well) %>%
#   #mutate(droplet_num = row_number()) %>%
#   #filter(channel == "NVO") %>%
#   ggplot(aes(y = amplitude, x = channel, fill = clusterAssignment, color = clusterAssignment)) +
#   geom_boxplot() +
#   #geom_point(position = position_jitterdodge(jitter.width = 0.2)) +
#   facet_grid(plate_col ~ plate_row, scales = "free") +
#   theme_light()
