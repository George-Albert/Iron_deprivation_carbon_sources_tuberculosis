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
  }
###################################################
### 1. Set working directory and create folders ###
###################################################


main_wd <- getwd()
setwd(main_wd)
input_dir <- "Analysis/Inputs/1_Raw_data"
output_dir <- "Analysis/Outputs"

###########################
####### 1. Functions ######
###########################
dcols=function(x){data.frame(colnames(x))}

###########################
####### 2. Load data ######
###########################
reads_company <- read.table(file.path(main_wd,"Analysis/Inputs/1_Raw_data/reads_company.txt"),check.names = F)
reads_new <- read.table(file.path(main_wd,"Analysis/Inputs/1_Raw_data/reads_new.txt"),check.names = F)

#Exponential phase

#*C1 <- Glycerol, Dextrose, NO LCFA
#*C5 <- Glycerol, No Dextrose, NO LCFA
#*C7 <- Glycerol, Dextrose, LCFA
#*C11 <- NO Glycerol, NO Dextrose, LCFA
#*C14 <- NULL

#Stationary phase

#*C2 <- Glycerol, Dextrose, no LCFA
#*C6 <- Glycerol, No Dextrose, NO LCFA
#*C8 <- Glycerol, Dextrose, LCFA
#*C10 <- Glycerol, NO Dextrose, LCFA
#*C12 <- NO Glycerol, NO Dextrose, LCFA
#*C15 <- NULL

#long Stationary phase

#*C10 <- Glycerol, NO Dextrose, LCFA
#*C13 <- NO Glycerol, NO Dextrose, LCFA
full_path_df <- read.table(file.path(input_dir,"Mapping_raw_data/list_of_fastq_paths/fastq_files_paths.txt"))
sample_name  <- read.table(file.path(input_dir,"Mapping_raw_data/list_of_fastq_paths/raw_fastq_names.txt"))

base_name <- basename(full_path_df$fastq_paths)
base_name <- substr(base_name,1,nchar(base_name)-11)
base_name <- base_name[-53]
base_name <- unique(base_name)

sample_name <- apply(sample_name,1,function(x) {substr(x,1,nchar(x)-11)})
sample_name <- sample_name[-53]
sample_name <- unique(sample_name)

# check congruence between full_path and sample name
all.equal(base_name,sample_name)
# TRUE

full_path_df <- substr(full_path_df$fastq_paths,1,nchar(full_path_df$fastq_paths)-11)
full_path_df <- full_path_df[-53]
full_path_df <- unique(full_path_df)

reads_company <- reads_company[,sample_name]
reads_new <- reads_new[,sample_name]

length(which(rownames(reads_company)!=rownames(reads_new[rownames(reads_company),])))
length(which(colnames(reads_company)!=colnames(reads_new)))

### change the samples names from HN00153515_C1_C2_Julio_2021
# 1 10_3MMmADN_F1C1 
# 2 10_3MMmADNmF1C1 
# 3 16_3MMmADN_F2C1 
# 4 16_3MMmADNmF2C1 
# 5 23_4MMmADN_F1C2 
# 6 23_4MMmADN_F2C2 
# 7 23_4MMmADNmF1C2
# 8 23_4MMmADNmF2C2
sample_number <-c(1:8) 
real_names <- c("10_3MMmADN_F1C1", "10_3MMmADNmF1C1", "16_3MMmADN_F2C1", "16_3MMmADNmF2C1", "23_4MMmADN_F1C2",
  "23_4MMmADN_F2C2","23_4MMmADNmF1C2", "23_4MMmADNmF2C2")
real_names <-paste0(real_names,"_rep1")

reads <- reads_new
colnames(reads)[9:16] <- real_names 
### Create metadata
#Lets name reads to reads_news

metadata=data.frame(Original_ID=substr(colnames(reads),1,(nchar(colnames(reads)))),Full_path=full_path_df)

{
  #Create Iron column
  metadata$Iron <- NA
  
  for (i in (1:nrow(metadata))) {
    
    if (grepl("mF", metadata[[1]][i], fixed = TRUE)) {
      
      metadata$Iron[i]= "YES"
      
    } else{
      metadata$Iron[i]= "NO"
    }
  }
  
  #Create Dextrose column
  metadata$Dextrose <- NA
  
  for (i in (1:nrow(metadata))) {
    
    if (grepl("C11", metadata[[1]][i], fixed = TRUE)||grepl("C12", metadata[[1]][i], fixed = TRUE)||
        grepl("C14", metadata[[1]][i], fixed = TRUE)||grepl("C15", metadata[[1]][i], fixed = TRUE)||
        grepl("C7", metadata[[1]][i], fixed = TRUE)||grepl("C8", metadata[[1]][i], fixed = TRUE)) {
      
      metadata$Dextrose[i]= "NO"
      
    } else{
      if (grepl("C1", metadata[[1]][i], fixed = TRUE)||grepl("C2", metadata[[1]][i], fixed = TRUE)){
        
        metadata$Dextrose[i]= "YES"
        
      } else{
        metadata$Dextrose[i]= "NO"
      }
      
    }
  }
  
  #Create Growth column
  
  metadata$Growth <- NA
  for (i in (1:nrow(metadata))) {
    
    if (grepl("C2", metadata[i,1], fixed = TRUE)||grepl("C6", metadata[i,1], fixed = TRUE)||
        grepl("C8", metadata[i,1], fixed = TRUE)||grepl("C12", metadata[i,1], fixed = TRUE)||
        grepl("C15", metadata[i,1], fixed = TRUE)) {
      
      metadata$Growth[i]= "STAT"
      
    } else{
      
      if (grepl("C1", metadata[i,1], fixed = TRUE)||grepl("C5", metadata[i,1], fixed = TRUE)||
          grepl("C7", metadata[i,1], fixed = TRUE)||grepl("C11", metadata[i,1], fixed = TRUE)||
          grepl("C14", metadata[i,1], fixed = TRUE)) {
        
        metadata$Growth[i]= "EXP"
        
      } else{
        metadata$Growth[i]= "STAT"
      }
      
    }
  }
  
  
  #Create LCFA column
  
  metadata$LCFA <- NA
  
  for (i in (1:nrow(metadata))) {
    if (grepl("C1$|C1_",ignore.case = TRUE,metadata[i,1])||grepl("C5$", ignore.case = TRUE,metadata[i,1])||
        grepl("C2", ignore.case = TRUE,metadata[i,1])||grepl("C6$", ignore.case = TRUE,metadata[i,1])||
        grepl("C14$", ignore.case = TRUE,metadata[i,1])||grepl("C15", ignore.case = TRUE,metadata[i,1])){
      
      metadata$LCFA[i]= "NO"
      
    } else{
      
      if (grepl("C7$", ignore.case = TRUE,metadata[i,1])|grepl("C8$", ignore.case = TRUE,metadata[i,1])|
          grepl("C11$",ignore.case = TRUE,metadata[i,1])|grepl("C12$",ignore.case = TRUE,metadata[i,1])){
        
        metadata$LCFA[i]= "YES"
        
      } else{
        
        metadata$LCFA[i]= "YES"
      }
    }
  }
  
  #Create Glycerol column
  
  metadata$Glycerol <- NA
  
  for (i in (1:nrow(metadata))) {
    if (grepl("C11$|_C11R", ignore.case = TRUE,metadata[i,1])||grepl("C12", metadata[i,1], fixed = TRUE)||
        grepl("C14", metadata[i,1], fixed = TRUE)||grepl("C15", metadata[i,1], fixed = TRUE)){
      
      metadata$Glycerol[i]= "NO"
      
    } else{
      
      metadata$Glycerol[i]= "YES"
      
    }
  }
}
  #Create metadata setup
  metadata$setup=paste0("Fe.",metadata$Iron,"_Dextrose.",metadata$Dextrose,"_Growth.",
                        metadata$Growth,"_LCFA.",metadata$LCFA,"_Glycerol.",metadata$Glycerol)

  length(which(metadata$Original_ID==colnames(reads))) 
  length(which(metadata$Original_ID!=colnames(reads))) 
  
  # reads=reads[,order(metadata$setup)]  
  # metadata=metadata[order(metadata$setup),]

  name_cul <- str_extract(metadata$Original_ID, "(C\\d+)")
  metadata$Culture <- name_cul
  
  rownames(metadata)=paste0(metadata$Culture,"_Fe_",metadata$Iron, "_", 
                            ave(paste0(metadata$Culture,"_Fe_",metadata$Iron), 
                                paste0(metadata$Culture,"_Fe_",metadata$Iron), FUN = seq_along))
  
  metadata$Sample_ID <- rownames(metadata)
  
  length(which(metadata$Original_ID!=colnames(reads)))
  # duplicated(colnames(reads))
  
  colnames(reads)=rownames(metadata)
 
  metadata=metadata[,c("Iron","Dextrose","Growth","LCFA","Glycerol","Culture",
                       "Sample_ID","Original_ID","setup","Full_path")]
 
  metadata <- metadata[order(as.numeric(gsub("C","",metadata$Culture))),]
  metadata$short_setup <- paste0(metadata$Culture,"_Fe_",metadata$Iron)
  reads <- reads[,rownames(metadata)]
  
  cult <- paste0("C",c(1:2,5:6,7:8,11:12))
  
  ### Set the colors for the cultures
  #C1 <- "orchid"
  #C2 <- "darkgreen"
  #C5 <- "blue"
  #C6 <- "red"
  #C7 <- "orange"
  #C8 <- "midnightblue"
  #C11 <- "green"
  #C12 <- "darkgoldenrod2
  color_vec <- c(C1 = "#f6a1a1", C2 = "#852121", C5 = "#ffc499", C6 = "#8f4a17", 
                 C7 = "#a9dea9", C8 = "#2c682c", C11 = "#8cd9ff", C12 = "#355f91", C14 ="#df82f0",C15 = "#c149d6")
  cultures <- metadata$short_setup
  
  ### Add color and fill columns to the metadata
  metadata$color <- color_vec[match(metadata$Culture, unique(metadata$Culture))]
  metadata$fill <- "white"
  metadata$fill[metadata$Iron == "YES"] <- (color_vec[match(metadata$Culture[metadata$Iron == "YES"],
                                                            unique(metadata$Culture[metadata$Iron == "YES"]))])
  
  dir.create("Analysis/Inputs/2_Processed_data/xlsx",showWarnings = F)
  dir.create("Analysis/Inputs/2_Processed_data/txt",showWarnings = F)
  dir.create("Analysis/Inputs/2_Processed_data/RDS",showWarnings = F)
  
  ### Save data
  write_xlsx(metadata,file.path(getwd(),"Analysis/Inputs/2_Processed_data/xlsx/metadata_66_samples.xlsx"))
  write.table(metadata,file.path(getwd(),"Analysis/Inputs/2_Processed_data/txt/metadata_66_samples.txt"))
  write.table(reads,file.path(getwd(),"Analysis/Inputs/2_Processed_data/txt/reads_66_samples.txt"))  
