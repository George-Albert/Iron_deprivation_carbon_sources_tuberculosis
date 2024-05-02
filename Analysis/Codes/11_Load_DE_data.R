#############################
### 0. Load dependencies. ###
#############################
{
  library(readxl)
  library(writexl)
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
contrast <- read.table(file.path(input_dir,"txt/contrast_matrix.txt"))

feat_path = file.path(input_dir,"txt/feature_data_filtered.txt")
feature_data = read.table(feat_path)

meta_path <- file.path(input_dir,"txt/metadata_32_samples.txt")
metadata = read.table(meta_path)

df_name_contrast <- read.table(file.path(input_dir,"txt/contrasts_nomenclature.txt"))
file_names <- paste0(df_name_contrast$title,".txt")
file_names <- file_names[c(1:16,44:47)]
current_folder <- file.path(output_dir,"Data/txt/th_0.05_th_size_0")

##################################################################
####### 4. Create a full list of the DE stats per condition ######
##################################################################
### Get full names including folder path
list.of.files = list.files(current_folder, full.names = TRUE)

### Keep only the basename (file names) matching dataframe column
clean_list <- list.of.files[basename(list.of.files) %in% file_names]
data_name <- basename(clean_list)
data_name <- substr(data_name,1,nchar(data_name)-4)

### Read the data
myData <- lapply(clean_list, read.table)
myData <- setNames(object = myData, data_name)

### Save the list
saveRDS(myData, file = file.path(input_dir,"RDS/Contrasts_stat_0.05.RDS")) 
# mydf <- readRDS("Contrasts_stat.RDS")






