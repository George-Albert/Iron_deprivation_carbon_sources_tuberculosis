########################
####### Depencies ######
########################
{
  library(stringr)
  library(tidyverse)
  library(dplyr)
  library(readxl)
  library(writexl)
  library(dplyr)
  library(writexl)
}
###################################################
### 1. Set working directory and create folders ###
###################################################

main_wd <- getwd()
setwd(main_wd)
input_dir <- "Analysis/Inputs/2_Processed_data"
output_dir <- "Analysis/Outputs"

###########################
####### 1. Functions ######
###########################
dcols=function(x){data.frame(colnames(x))}

###########################
####### 2. Load data ######
###########################
metadata_66_samples <- read.table(file.path(input_dir,"txt/metadata_66_samples.txt"),check.names = F)
reads_66_samples <- read.table(file.path(input_dir,"txt/reads_66_samples.txt"),check.names = F)
metadata_old <- read.table(file.path(input_dir,"txt/metadata_32_samples_old.txt"),check.names = F)
# reads_C11_C12 <- read.table(file.path(input_dir,"reads_C11_C12_samples.txt"),check.names = F)

### Subset the data
# metadata_66_samples$short_setup <- paste0(metadata_66_samples$Culture,"_Fe_",metadata_66_samples$Iron)
metadata <- metadata_66_samples[which(metadata_66_samples$Original_ID %in% metadata_old$Original_ID),]

### Keep C12 old (HN00166209_C11_C12_C2r_Marzo_2022)
### C12--> old; C11--> New.
met_C11_C12 <- metadata_66_samples[which(metadata_66_samples$Culture %in% c("C12")),]
met_C11_C12 <-  met_C11_C12[1:4,]

metadata[which(metadata$Culture %in% c("C12")),] <- met_C11_C12
rownames(metadata) <- metadata$Sample_ID
### to extract C12
# metadata_C11_C12 <- metadata_66_samples[which(metadata_66_samples$Culture=="C12"),]
# metadata_C11_C12 <- metadata_C11_C12[grep(c("1$|2$"),metadata_C11_C12$Sample_ID),]
# metadata[which(metadata$Culture=="C12"),] <- metadata_C11_C12
# rownames(metadata) <- metadata$Sample_ID
# metadata <- metadata[order(metadata$short_setup),]

# metadata_old <- metadata_old[order(match(metadata_old$Original_ID,metadata$Original_ID)),]
# 
# length(which(metadata$Original_ID != metadata_old$Original_ID))

reads_32 <- reads_66_samples[,rownames(metadata)]

length(which(colnames(reads_32) != rownames(metadata)))

### Save the data
write.table(metadata,file.path(input_dir,"txt/metadata_32_samples.txt"))
write.table(reads_32,file.path(input_dir,"txt/reads_32_samples.txt"))



