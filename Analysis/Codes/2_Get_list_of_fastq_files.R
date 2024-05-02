
############################################################################
############################################################################
###                                                                      ###
###                    GETTING A LIST OF FASTQ FILES                     ###
###                                                                      ###
############################################################################
############################################################################

getwd()
list_dir <- "Analysis/Inputs/1_Raw_data/Mapping_raw_data/list_of_fastq_paths"
dir.create(list_dir)

full_list_paths <- read.delim(file.path(getwd(),"Analysis/Inputs/1_Raw_data/Mapping_raw_data/full_path_to_NAS.txt"))

### Select the paths where the fastq are in
names(full_list_paths) <- "paths"
# indices <- grepl("RawFASTQ", full_list_paths$paths)
# raw_fastq_path <- full_list_paths[indices]

### split the paths based on "\"
raw_fastq_path <- strsplit(full_list_paths[[1]], "\\\"")
### Get the indices
indices <- grepl("RawFASTQ", raw_fastq_path)
raw_fastq_path_filt <- raw_fastq_path[indices]
### Substitute spaces by "_"
raw_fastq_path1 <- lapply(raw_fastq_path_filt, function(x) gsub(" ", "_", x))
### Substitute "-" by "_"
raw_fastq_path <- lapply(raw_fastq_path1, function(x) gsub("-", "_", x))
raw_fastq_path_df <- data.frame(do.call(cbind, raw_fastq_path))
### Create a matrix
raw_fastq_path_df <- t(raw_fastq_path_df)
rownames(raw_fastq_path_df) <- c(1:133)
colnames(raw_fastq_path_df) <- "fastq_paths"

### Divide raw_fastq_path in parts using "/" as separators
# parts <- strsplit(raw_fastq_path[[1]], "/")
parts <-lapply(raw_fastq_path,FUN=function(x) unlist(strsplit(x, "/")))

### Convert the list of character vectors in a data frame
parts_df <- data.frame(do.call(rbind, parts))

### Remove the 1st and the 3rd column 

parts_df_rm <- parts_df[,c("X2","X4")]
colnames(parts_df_rm) <- c("folders_names","fastq_names")

### create a list with the folders names
folders_name <- unique(parts_df_rm$folders_names)
raw_fastq_names <- parts_df_rm$fastq_names

# (other way to do the same as above)Extract the filename from a complete path in a dataframe and store it 
# in a new variable named 'raw_fastq_names'
raw_fastq_names <- sub(".*/([^/]+)$", "\\1", raw_fastq_path_df)
# list_subchains <- strsplit(raw_fastq_names, "\"")
# raw_fastq_names <-t(data.frame(do.call(cbind, list_subchains)))
# rownames(raw_fastq_names) <- c(1:nrow(raw_fastq_names))
# colnames(raw_fastq_names) <- "file_name"

list_by_folders <- split(parts_df_rm,parts_df_rm$folders_names,drop = T)

### Define a function to delete 1st column of each data frame
delete_first_col <- function(df) {
  df[, -1]
}

### Apply the function to each df of the list using lapply
list_by_folders <- lapply(list_by_folders, delete_first_col)

write.table(parts_df_rm,file.path(list_dir,"folder_and_fastq_names.txt"))
write.table(folders_name,file.path(list_dir,"folders_names.txt"))
write.table(raw_fastq_names,file.path(list_dir,"raw_fastq_names.txt"))
write.table(raw_fastq_path_df,file.path(list_dir,"fastq_files_paths.txt"))

### Save a list of df
saveRDS(list_by_folders, file.path(list_dir,"list_by_folders.RDS"))



