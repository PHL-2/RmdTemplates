########### Config file for wastewater testing using the ddPCR ###########

##### File paths #####

## Filepath of exported ddPCR analysis file (Please add quotation marks)

dat_fp = "//Volumes/Export"

###### Volumes used in experiment ######

## Wastewater volume (mL)
ww_vol_mL = 40

## MHV added to waste water (uL)
mhv_ww_input_uL = 8

## Control matrix volume (mL)
control_vol_mL = 20



#### Concentration step ####

## Concentration: Sample volume used (uL)
init_concentration_vol_uL = 34375

## Concentration: Elution volume (uL)
end_concentration_vol_uL = 400


#### Extraction step ####

## Extraction: Sample volume used (uL)
extraction_vol_uL = 400

## Extraction: elution volume (uL)
elution_vol_uL = 50

#### ddPCR ####

## ddPCR reaction volume (uL)
rxn_vol_uL = 20

## ddPCR: sample input volume (uL)
input_vol_uL = 9


##### Stock information #####
## Interal control stock information
icInfo = "MHV"

## SC2 stock information
posCtrl = "ZeptoSC2"


## Filepath to wastewater metadata in aux folder
wwMeta_fp = file.path(dirname(here::here()), "aux_files")

## Filepath to wastewater metadata in aux folder
ww_site_meta_fp = list.files(path = file.path(dirname(here()), "aux_files"),
                             pattern = "wastewater_specific_metadata.csv",
                             full.names = TRUE,
                             recursive = TRUE)

ww_nwss_meta_fp = list.files(path = file.path(dirname(here()), "aux_files"),
                             pattern = "wastewater_data_submission_metadata.csv",
                             full.names = TRUE,
                             recursive = TRUE)

