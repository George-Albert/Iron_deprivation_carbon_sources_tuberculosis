#############################
### 0. Load dependencies. ###
#############################

{
  library(tidyverse)
  library(dplyr)
  
}

############################
### 1. Declare functions ### 
############################
dcols <- function(df){data.frame(colnames(df))}

###################################################
### 2. Set working directory and create folders ###
###################################################


main_wd <- getwd()
setwd(main_wd)

input_dir <- "Analysis/Inputs/2_Processed_data"
output_dir <- "Analysis/Outputs"

####################
### 3. Load Data ###
####################

contrast <- read.table(file.path(input_dir,"txt/contrast_matrix.txt"))
df_name_contrast <- read.table(file.path(input_dir,"txt/contrasts_nomenclature.txt"))
myData <- readRDS(file.path(input_dir,"RDS/Contrasts_stat_0.05.RDS"))
feature_data = read.table(file.path(input_dir,"txt/feature_data_filtered.txt"))


#################################################################
##              Create LogFC, BH and lfcSE tables              ##
#################################################################

{
  
  list_name <- c( "Effect_of_Fe_in_response_to_Growth_arrest_G_L",    
                  "Effect_of_Fe_in_response_to_Growth_arrest_G_D",
                  "Effect_of_Fe_in_response_to_Growth_arrest_G",         
                  "Effect_of_Fe_in_response_to_Growth_arrest_L",             
                  "Iron_effect_in_EXP_G_D",                       
                  "Iron_effect_in_EXP_G_L",                           
                  "Iron_effect_in_EXP_G",                                
                  "Iron_effect_in_EXP_L",                                    
                  "Iron_effect_in_STAT_G_D",                      
                  "Iron_effect_in_STAT_G_L",                          
                  "Iron_effect_in_STAT_G",                               
                  "Iron_effect_in_STAT_L",
                  "Phase_effects_with_Fe_G_D",                    
                  "Phase_effects_with_Fe_G_L",                        
                  "Phase_effects_with_Fe_G",                             
                  "Phase_effects_with_Fe_L",                                 
                  "Phase_effects_without_Fe_G_D",                 
                  "Phase_effects_without_Fe_G_L",                     
                  "Phase_effects_without_Fe_G",                          
                  "Phase_effects_without_Fe_L")
  
  names(myData) <- list_name
  
  lapply(names(myData), function(x) assign(x,myData[[x]],envir = .GlobalEnv))  

length(which(rownames(feature_data)!= rownames(Iron_effect_in_EXP_G_D)))
### Create log Fold Change df
LogFC_df <- data.frame(Iron_effect_in_EXP_G_D$log2FoldChange_shrunken, Iron_effect_in_EXP_G$log2FoldChange_shrunken,
                       Iron_effect_in_EXP_G_L$log2FoldChange_shrunken, Iron_effect_in_EXP_L$log2FoldChange_shrunken,
                       Iron_effect_in_STAT_G_D$log2FoldChange_shrunken, Iron_effect_in_STAT_G$log2FoldChange_shrunken,
                       Iron_effect_in_STAT_G_L$log2FoldChange_shrunken, Iron_effect_in_STAT_L$log2FoldChange_shrunken,
                       Phase_effects_with_Fe_G_D$log2FoldChange_shrunken, Phase_effects_with_Fe_G$log2FoldChange_shrunken,
                       Phase_effects_with_Fe_G_L$log2FoldChange_shrunken, Phase_effects_with_Fe_L$log2FoldChange_shrunken,
                       Phase_effects_without_Fe_G_D$log2FoldChange_shrunken, Phase_effects_without_Fe_G$log2FoldChange_shrunken,
                       Phase_effects_without_Fe_G_L$log2FoldChange_shrunken, Phase_effects_without_Fe_L$log2FoldChange_shrunken,
                       Effect_of_Fe_in_response_to_Growth_arrest_G_D$log2FoldChange_shrunken,
                       Effect_of_Fe_in_response_to_Growth_arrest_G$log2FoldChange_shrunken,
                       Effect_of_Fe_in_response_to_Growth_arrest_G_L$log2FoldChange_shrunken,
                       Effect_of_Fe_in_response_to_Growth_arrest_L$log2FoldChange_shrunken)

rownames(LogFC_df) <- rownames(Iron_effect_in_EXP_G_D)

### Create BH df
BH_df <- data.frame(Iron_effect_in_EXP_G_D$BH, Iron_effect_in_EXP_G$BH,
                    Iron_effect_in_EXP_G_L$BH, Iron_effect_in_EXP_L$BH,
                    Iron_effect_in_STAT_G_D$BH, Iron_effect_in_STAT_G$BH,
                    Iron_effect_in_STAT_G_L$BH, Iron_effect_in_STAT_L$BH,
                    Phase_effects_with_Fe_G_D$BH, Phase_effects_with_Fe_G$BH,
                    Phase_effects_with_Fe_G_L$BH, Phase_effects_with_Fe_L$BH,
                    Phase_effects_without_Fe_G_D$BH, Phase_effects_without_Fe_G$BH,
                    Phase_effects_without_Fe_G_L$BH, Phase_effects_without_Fe_L$BH,
                    Effect_of_Fe_in_response_to_Growth_arrest_G_D$BH,
                    Effect_of_Fe_in_response_to_Growth_arrest_G$BH,
                    Effect_of_Fe_in_response_to_Growth_arrest_G_L$BH,
                    Effect_of_Fe_in_response_to_Growth_arrest_L$BH)

rownames(BH_df) <- rownames(Iron_effect_in_EXP_G_D)

### Create log Fold Change Standard Error df
lfcSE_df <- data.frame(Iron_effect_in_EXP_G_D$lfcSE_shrunken, Iron_effect_in_EXP_G$lfcSE_shrunken,
                       Iron_effect_in_EXP_G_L$lfcSE_shrunken, Iron_effect_in_EXP_L$lfcSE_shrunken,
                       Iron_effect_in_STAT_G_D$lfcSE_shrunken, Iron_effect_in_STAT_G$lfcSE_shrunken,
                       Iron_effect_in_STAT_G_L$lfcSE_shrunken, Iron_effect_in_STAT_L$lfcSE_shrunken,
                       Phase_effects_with_Fe_G_D$lfcSE_shrunken, Phase_effects_with_Fe_G$lfcSE_shrunken,
                       Phase_effects_with_Fe_G_L$lfcSE_shrunken, Phase_effects_with_Fe_L$lfcSE_shrunken,
                       Phase_effects_without_Fe_G_D$lfcSE_shrunken, Phase_effects_without_Fe_G$lfcSE_shrunken,
                       Phase_effects_without_Fe_G_L$lfcSE_shrunken, Phase_effects_without_Fe_L$lfcSE_shrunken,
                       Effect_of_Fe_in_response_to_Growth_arrest_G_D$lfcSE_shrunken,
                       Effect_of_Fe_in_response_to_Growth_arrest_G$lfcSE_shrunken,
                       Effect_of_Fe_in_response_to_Growth_arrest_G_L$lfcSE_shrunken,
                       Effect_of_Fe_in_response_to_Growth_arrest_L$lfcSE_shrunken)

rownames(lfcSE_df) <- rownames(Iron_effect_in_EXP_G_D)
}
### Check congruence
length(which(rownames(feature_data)!= rownames(LogFC_df)))
length(which(rownames(feature_data)!= rownames(BH_df)))
length(which(rownames(feature_data)!= rownames(lfcSE_df)))

### Save the data

write.table(LogFC_df, file= file.path(input_dir,"txt/LogFC_0.05.txt"))
write.table(BH_df, file= file.path(input_dir,"txt/BH_0.05.txt"))
write.table(lfcSE_df, file= file.path(input_dir,"txt/lfcSE_0.05.txt"))

