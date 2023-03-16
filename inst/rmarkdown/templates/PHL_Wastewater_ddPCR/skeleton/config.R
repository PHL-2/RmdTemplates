

########### Config file for wastewater testing using the ddPCR ########### 

##### File paths ##### 

## Filepath of exported ddPCR analysis file (Please add quotation marks)
dat_fp = "C:/Users/christopher.gu/OneDrive - City of Philadelphia/Wastewater/build/data_csv"


###### Volumes used in experiment ######

## Wastewater volume (mL)
ww_vol_mL = 40

## MHV added to waste water (uL)
mhv_ww_input_uL = 8

## Control matrix volume (mL)
control_vol_mL = 25

## Control matrix volume date info ([year][month][day])
control_matrix_info = "202212XX Wastewater"

## MHV added to control matrix (uL)
mhv_to_control_uL = 5

## Sc2 added to control matrix (uL)
sc2_to_control_uL = 100



#### Concentration step ####

## Concentration: Sample volume used (uL)
init_concentration_vol_uL = 9600

## Concentration: Elution volume (uL)
end_concentration_vol_uL = 500



#### Extraction step ####

## Extraction: Sample volume used (uL)
extraction_vol_uL = 400

## Extraction: elution volume (uL)
elution_vol_uL = 100

#### ddPCR ####

## ddPCR reaction volume (uL)
rxn_vol_uL = 20

## ddPCR: sample input volume (uL)
input_vol_uL = 9


##### Stock information ##### 

####MHV stock information
MHV_stock_name = "UPennMHV"
MHV_concentration_copiesuL = 69620.3

#SC2 stock information
SC2_stock_name = "ZeptoSC2"
SC2_concentration_copiesuL = 200


