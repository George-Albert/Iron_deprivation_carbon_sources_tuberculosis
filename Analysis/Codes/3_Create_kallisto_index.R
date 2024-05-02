########################
####### Depencies ######
########################
library(tidyverse)
library(Rgff)
library(Rsubread)
library(gread)
#######################################
####### 1. Set working directory ######
#######################################


main_wd <- getwd()
setwd(main_wd)
###########################
####### 1. Functions ######
###########################

###########################
####### 2. Load data ######
###########################
input_dir <- "Analysis/Inputs/1_Raw_data/Mapping_raw_data"
feature_ASM19595v2 <- read.delim(file.path(input_dir,"GCF_000195955.2_ASM19595v2_feature_table.txt"))
#feature_count
#Tab-delimited text file reporting counts of gene, RNA, CDS, and similar features, based on data reported in 
#the *_feature_table.txt.gz file (see below). Separate counts are provided for different sets of sequences in 
#the assembly corresponding to the primary assembly, non-nuclear assembly, all alt-loci sequences, and all 
#patch scaffolds.
feature_count <- read.delim(file.path(input_dir,"GCF_000195955.2_ASM19595v2_feature_count.txt"))

### Rearranging the feature data and filtering by gene
feature_ASM19595v2 <- feature_ASM19595v2 %>% relocate("locus_tag","GeneID")
feature_data <- feature_ASM19595v2[which(feature_ASM19595v2$X..feature == "gene"),]
feature_data <- feature_data[ , colSums(is.na(feature_data))==0] # remove nan columns
feature_data <- feature_data[ , !colSums(feature_data=="") == nrow(feature_data)] # remove empty columns
colnames(feature_data)[3] <- "INSDC_feature"
write.table(feature_data,file.path("Inputs/2_Processed_data/feature_data_genes.txt"))

### Check the amount of genes
char_to_numeric <- as.numeric(feature_count[2:nrow(feature_count),]$Unique.Ids)
sum(char_to_numeric,na.rm = T)
### 4008 genes

#########################
#### Kallisto index #####
######################### 

### Constructing index
genomic_fasta_ncbi <- file.path(file.path(input_dir,"gene.fna"))

command = paste0("kallisto index -i index_ncbi_ASM19595v2.idx ",genomic_fasta_ncbi)
system(command)

output_dir <-file.path("Analysis/Outputs")
dir.create(file.path(output_dir,"1_Kallisto/index"),recursive = T,showWarnings = F)

### locate the index file
index_file <- list.files(pattern = "idx")
index_file
### copy to de dir we create for
file.copy(from = index_file,
          to = paste0(file.path(output_dir,"1_Kallisto/index", index_file)))
### remove the copy from the working directory
file.remove(from = index_file)

### inspect index
index <- file.path(output_dir,"1_Kallisto/index", index_file)
command=paste0("kallisto inspect ",index," gtf=",gtf_path)
# inspect_index <- system(command,intern = T)
capture.output(system(command,intern = T),file = paste0(file.path(output_dir,"1_Kallisto/index/inspect_index")))


















