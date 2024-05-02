
########################
####### Depencies ######
########################
library(tidyverse)

#######################################
####### 1. Set working directory ######
#######################################
main_wd <- getwd()
setwd(main_wd)
### Create folder`s names ###

doc_folder <- "Docs"
analysis_folder <- "Analysis"
codes_folder <- "Codes"
inputs_folder <- "Inputs"
outputs_folder <- "Outputs"
raw_folder <- "1_Raw_data"
processed_folder <- "2_Processed_data"

### Create folders list
folder_list <- list(codes_folder,inputs_folder,outputs_folder,raw_folder,processed_folder)

dir.create(doc_folder,showWarnings = F)
### Create Codes, Inputs and Outputs folders inside Analysis ####
sapply(folder_list[1:3], FUN=function(x) dir.create(file.path(analysis_folder,x),showWarnings = F,recursive=T))
### Create Raw and Processed folders inside Inputs ####
sapply(folder_list[4:5], FUN=function(x) dir.create(file.path(analysis_folder,inputs_folder,x),showWarnings = F,recursive=T))













