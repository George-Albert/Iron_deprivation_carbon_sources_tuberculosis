#############################
### 0. Load dependencies. ###
#############################
library(openxlsx)

safe_path_name <- function(x) {
  x <- gsub("->", "to", x, fixed = TRUE)
  x <- gsub("[<>:\"/\\\\|?*]", "_", x)
  x <- gsub("\\s+", "_", x)
  x
}

###################################################
### 2. Set working directory and create folders ###
###################################################
main_wd <- getwd()
setwd(main_wd)
input_dir <- "Inputs/002_Processed_data"
output_dir <- "Outputs"

###########################
####### 3. Load data ######
###########################
contrast <- read.table(file.path(input_dir,"txt/contrast_matrix.txt"))

feat_path = file.path(input_dir,"txt/feature_data_filtered.txt")
feature_data = read.table(feat_path)

meta_path <- file.path(input_dir,"txt/metadata_32_samples.txt")
metadata = read.table(meta_path)

df_name_contrast <- read.table(file.path(input_dir,"txt/contrasts_nomenclature.txt"))
file_names <- paste0(safe_path_name(df_name_contrast$title),".txt")
file_names <- file_names[c(1:16,44:47)]
current_folder <- file.path(output_dir,"Data/txt/th_0.05_th_size_0")

##################################################################
####### 4. Create a full list of the DE stats per condition ######
##################################################################
### Get full names including folder path
list.of.files = list.files(current_folder, full.names = TRUE)

### Keep only the basename (file names) matching dataframe column
clean_list <- list.of.files[basename(list.of.files) %in% file_names]
stopifnot(length(clean_list) == length(file_names))
data_name <- basename(clean_list)
data_name <- substr(data_name,1,nchar(data_name)-4)

### Read the data
myData <- lapply(clean_list, read.table)
myData <- setNames(object = myData, data_name)

### Save the list
saveRDS(myData, file = file.path(input_dir,"RDS/Contrasts_stat_0.05.RDS"))
# mydf <- readRDS("Contrasts_stat.RDS")

### I also save the data in an Excel file with short names for the contrasts because
### excel does not allow long names for the sheets.
short_names <- names(myData)

short_names <- gsub("Effect_of_Fe_in_response_to_Growth_arrest_", "FeResp_GA_", short_names)
short_names <- gsub("Iron_effect_in_", "FeEff_", short_names)
short_names <- gsub("Phase_effects_with_Fe_", "Phase_Fe+_", short_names)
short_names <- gsub("Phase_effects_without_Fe_", "Phase_Fe-_", short_names)

short_names <- gsub("GLYCEROL-DEXTROSE", "G-D", short_names)
short_names <- gsub("GLYCEROL-LCFA", "G-LCFA", short_names)
short_names <- gsub("GLYCEROL", "G", short_names)

short_names <- gsub("[()]", "", short_names)

names(myData) <- short_names
write.xlsx(myData, "Table_1_Section_2.xlsx")
