
#############################
### 0. Load dependencies. ###
#############################
{
  library(readxl)
  library(writexl)
  library(openxlsx)
}

############################
### 1. Declare functions ### 
############################
dcols=function(x){data.frame(colnames(x))}

###################################################
### 2. Set working directory and create folders ###
###################################################

main_wd <- getwd()
setwd(main_wd)
input_dir <- "Analysis/Inputs/2_Processed_data"
output_dir <- "Analysis/Outputs"

###########################
####### 3. Load data ######
###########################
contrast <- read.table(file.path(input_dir,"txt/contrasts_nomenclature.txt"))

feat_path = file.path(input_dir,"txt/feature_data_filtered.txt")
feature_data = read.table(feat_path)

meta_path <- file.path(input_dir,"txt/metadata_32_samples.txt")
metadata = read.table(meta_path)

MyData <-readRDS(file = file.path(input_dir,"RDS/Contrasts_stat_0.05.RDS")) 
names(MyData)                    

# Glycerol+dextrosa en Exp
# C2 --> Glycerol+dextrosa en Stat
# C5 --> Glycerol en Exp
# C6 --> Glycerol en Stat
# C7 --> Glycerol+Lipidos en Exp
# C8 --> Glycerol+Lipidos en Stat
# C11 --> Lípidos en Exp
# C12 --> Lípidos en Stat

data_to_save <- MyData[1:12]
names(data_to_save) <- c( "Interaction_C7_C8_G_LCFA", "Interaction_C1_C2_G_D","Interaction_C5_C6_G", "Interaction_C11_C12_LCFA",            
                    "C1_G_D_EXP" , "C7_G_LCFA_EXP","C5_G_EXP","C11_LCFA_EXP" ,                                   
                    "C2_G_D_STAT" , "C8_G_LCFA_STAT","C6_G_STAT","C12_LCFA_STAT" ) 
# Crear un objeto de libro Excel
wb <- createWorkbook()

# Iterar sobre la lista de data frames y agregar cada uno como una hoja
for (nombre_df in names(data_to_save)) {
  addWorksheet(wb, sheetName = nombre_df)
  writeData(wb, sheet = nombre_df, x = data_to_save[[nombre_df]])
}

### Save the data
saveWorkbook(wb, file = file.path(input_dir,"xlsx","Iron_effects.xlsx"), overwrite = TRUE)
