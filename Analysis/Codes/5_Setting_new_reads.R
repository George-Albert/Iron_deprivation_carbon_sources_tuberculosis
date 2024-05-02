
########################
####### Depencies ######
########################
library(stringr)
library(tidyverse)
library(dplyr)
library(xlsx)
library(purrr)
library(ggridges)
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

###########################
####### 2. Load data ######
###########################

### Read Kallisto´s files         
list_tsv <- list.files(file.path(output_dir,"2_Reports/2_kallisto_output"),pattern = ".tsv", recursive = T)
list_txt <- list.files(file.path(output_dir,"2_Reports/reads_from_tximport"),pattern = "\\.txt$", recursive = T)

length(list_tsv)
length(list_txt)

### Load each file and keep it in a list
list_of_files <- lapply(list_txt, function(x) read.table(file.path(output_dir,"2_Reports/reads_from_tximport",x), header = TRUE))
length(list_of_files)

### Obtain the sample name
sample_name <- basename(list_txt)
sample_name <- substr(sample_name,1,nchar(sample_name)-4)

data <- list_of_files
data <- setNames(data,sample_name)

### est count data frame ### 
est_count_list <- lapply(1:length(data), function(x) data[[x]][1])
est_count_list <- setNames(est_count_list,sample_name)
est_count_df <- do.call(cbind, est_count_list)
colnames(est_count_df) <- sub("^X", "", colnames(est_count_df))
est_count_df <- est_count_df %>% 
  add_column(RV_code = rownames(est_count_df),.before =est_count_df$'03_12_AN_F1C14') 
est_count_df <- est_count_df[order(est_count_df$RV_code),]
est_count_df$RV_code <- NULL

write.table(est_count_df, file = file.path(input_dir,"reads_new.txt"),col.names = T)

######################################
### Prepare Reads from the company ###
######################################
full_list_paths <- read.delim(file.path(input_dir,"Mapping_raw_data/full_path_to_NAS.txt"))

### Select the paths where the fastq are in
names(full_list_paths) <- "paths"

xlsx_path <- strsplit(full_list_paths[[1]], "\\\"")
indices <- grepl("result_RNAseq_excel", xlsx_path)
xlsx_path_filt <- xlsx_path[indices]
xlsx_path1 <- lapply(xlsx_path_filt, function(x) gsub(" ", "_", x))
xlsx_path <- lapply(xlsx_path1, function(x) gsub("-", "_", x))
xlsx_path_df <- data.frame(do.call(cbind, xlsx_path))
xlsx_path_df <- t(xlsx_path_df)
rownames(xlsx_path_df) <- c(1:8)
colnames(xlsx_path_df) <- "excel_paths"

### Create the folders name
NAS_pass_dir <- file.path(input_dir,"passNAS.txt")
my_files <- xlsx_path_df
folders_name <- read.table(file.path(input_dir,"Mapping_raw_data/list_of_fastq_paths/folders_names.txt"))
folders_name1 <- data.frame(folders_name=apply(folders_name,1,FUN = function(x) gsub(" ", "_", gsub("-", "_", x))))

for(i in c(1:(dim(my_files)[1]))){
  
  remote_file <- paste0("/volume1/Raw_data/Colab/",xlsx_path_df[i])
  ### 1_result_RNAseq_excel is the folder with the expression profiles from Macrogen 
  local_folder <- file.path(input_dir,"1_result_RNAseq_excel",folders_name1[i,])
  
  dir.create(local_folder,recursive = T,showWarnings = F)
  
  ### type your NASS credentials
  ssh_command1 <- paste0("sshpass -f ",NAS_pass_dir," scp -O -P 28 -rp JCardenas@155.210.137.80:",remote_file," ",local_folder)

  system(ssh_command1)
}

list_xlsx <- list.files(file.path(input_dir,"1_result_RNAseq_excel"),pattern = "\\.xlsx$", recursive = T)
list_of_xlsxfiles <- lapply(list_xlsx, function(x) read.xlsx(file.path(input_dir,"1_result_RNAseq_excel",x),sheetIndex = 1, header = TRUE))
length(list_of_xlsxfiles)

extract_reads_from_excel <- function(df_list) {
  # Iterate over the list of df
  lapply(df_list, function(df) {
    # Select the columns that contain "Read_Count" in the name
    df %>% select(contains("Read_Count"),contains("locus_tag"))%>% 
      rename_all(~sub("_Read_Count_", "", df))
  })
    # Combine the df of the list by columns
    # df <- bind_cols() 
}
  
  # reads_df <- extract_reads_from_excel(list_of_xlsxfiles)
unique_colnames <- unique(unlist(map(list_of_xlsxfiles, names)))

list_df <- lapply(list_of_xlsxfiles, function(df){
  df %>%
    rename_at(vars(contains("Locus_tag")), ~"locus_tag") %>%
    select(Gene_ID, locus_tag, contains("Read_Count")) %>%
    rename_with(~str_remove(., "_Read_Count"), contains("_Read_Count")) %>%
    rename_with(~str_remove(., "X"), contains("X")) %>%
    arrange(df,locus_tag)
   
  })

unique_colnames <- unique(unlist(map(list_df, names)))

for (i in seq_along(list_df)) {
  # Assign the name to each df
  df_name <- paste0("df_", i)
  
  # Assign the df to an object with an specific name in the global environment
  assign(df_name, list_df[[i]])
  # colnames(get(df_name)) <- colnames(df_list[[i]])
}

### Remove "gene-" in locus_tag column
RV_codes=substr(df_1[,"locus_tag"],6,nchar(df_1[,"locus_tag"]))
df_1[,"locus_tag"]=RV_codes

### Sort df_5 
df_5 <- df_5[order(df_5$locus_tag),]
### rename locus_tag as RV_codes
colnames(df_5)[2] <- "RV_codes"
rownames(df_5) <- df_5$RV_codes
df_5[,c("Gene_ID","RV_codes")] <- NULL

### Substitute the . by the RV_code in the df_2
df_2[which(df_2$locus_tag=="."),"locus_tag"] <- c("RV2280","RV3476c")

### Sort all df by RV_codes 
df_name <- paste0("df_", 1:8)

for (name in df_name){
  if (name == "df_5") {
    next
  }
  print(name)
  df <- get(name)
  df <- df[order(df$locus_tag),]
  ### rename locus_tag as RV_codes
  print("changing names")
  colnames(df)[2] <- "RV_codes"
  rownames(df) <- df$RV_codes
  print("taking the common values")
  common <- intersect(rownames(df), rownames(df_5))
  df <- df[common,]
  df[,c("Gene_ID","RV_codes")] <- NULL
  df <- df[,order(names(df))]
  
  assign(name,df)
}
for (name in df_name){
  dif <- length(which(rownames(name)!=rownames(df_1)))
  print(paste0("There are ",dif, " between ", name, " and df_1"))
  }

list_new_df <- lapply(df_name, function(name) get(name))
reads_company <- do.call(cbind, list_new_df)

write.table(reads_company, file = file.path(main_wd,"Analysis/Inputs/1_Raw_data/reads_company.txt"),col.names = T)

#############################################################
### Correlate the new reads with the reads of the company ###
#############################################################

### Load read news
reads_news <- read.table(file.path(input_dir,"reads_new.txt"))

reads_news <- reads_news[order(rownames(reads_news)),]
reads_company <- reads_company[order(rownames(reads_company)),]
reads_news <- reads_news[rownames(reads_company),]

length(which(rownames(reads_news)!= rownames(reads_company)))

column_names <- colnames(reads_news)
column_names <- gsub("^X","",column_names)
colnames(reads_news) <- column_names
  
reads_news <- reads_news[,order(colnames(reads_news))]
reads_company <- reads_company[,order(names(reads_company))]

length(which(colnames(reads_news)!= colnames(reads_company)))

corr_matrix <- cor(reads_news,reads_company,method = "pearson")
corr_matrix_row <- sapply(1:nrow(as.matrix(reads_news)), function(i) cor(as.matrix(reads_news)[i,],
                                                                        as.matrix(reads_company)[i,],method = "pearson"))
corr_matrix <- data.frame(corr_matrix)
diagonal <- corr_matrix[row(corr_matrix)==col(corr_matrix)]
diagonal <- data.frame(corr_value=diagonal)
rownames(diagonal) <- colnames(reads_news)

corr_matrix_row <- data.frame(corr_value=corr_matrix_row)
rownames(corr_matrix_row) <- rownames(reads_news)
corr_matrix_row[which(is.na(corr_matrix_row)),] <- 0

#                         ====================================
#                         === Density Plot: each CONDITION ===
#                         ====================================
{
  dens1 <- corr_matrix_row
  
  dens2 <- diagonal
}

# write.table(dens,file="processed_all_JAC/dens_table.txt")

density <- function(data,title){
  
  frame <- data.frame(data)
  corr_value <- frame$corr_value
  y <- factor(frame$corr_value)
  sam <- reorder(y,corr_value)
  # color <- c("#00AFBB", "#E7B800","#50486D", "#FC4E07","#FFA373")
  color <- "#00AFBB"
  
  ggplot(frame, aes(x=corr_value)) +
    geom_density(alpha = 0.6) +
    scale_fill_manual(values = color)+
    labs(x="Correlation", y = "Density",title=title)+
    theme_ridges(font_size = 8,center_axis_labels = TRUE,grid = FALSE)+
    theme(plot.title = element_text(hjust = 0.5))
  
  
}

dir.create(file.path(output_dir,"3_Correlation_reads_old_vs_reads_new"),recursive = T,showWarnings = F)
pdf(file = file.path(output_dir,"3_Correlation_reads_old_vs_reads_new","dens_plot_correlation.pdf"),width=9,height=5)

density(dens2,"Density Plot by sample")
density(dens1,"Density Plot by gene")

dev.off()
