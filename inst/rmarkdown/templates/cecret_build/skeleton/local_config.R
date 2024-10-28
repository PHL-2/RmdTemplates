##### Config file for Cecret #####

run_uploaded_2_basespace <- TRUE # set this to true if the run was uploaded to BaseSpace and the data was not manually transferred to a local folder

##############
# 1_select_sample_n_metadata variables
##############

create_sample_replicates <- 4 #enter a positive integer for the number of biological replicates created during concentration and extraction

# set this to TRUE to copy the platemap to the shared drive
copy_platemap <- FALSE

sample_type_acronym <- "WW" #use WW for wastewater samples

prj_description <- "COVIDSeq" #no spaces, should be the same as the R project

# number of unspecified environmental swabs to add to plate
#enviro_number <- 10

##############
# 4_run_nextflow
##############

ec2_tmp_fp <- "~/tmp_bs_dl"

update_pangolin_dataset <- TRUE

update_freyja_and_cecret_pipeline <- TRUE

cecret_version <- "master"
