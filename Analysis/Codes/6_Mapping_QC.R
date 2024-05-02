
########################
####### Depencies ######
########################
library(tidyverse)
library(tximportData)
library(tximport)
library(Biostrings)

###################################################
### 1. Set working directory and create folders ###
###################################################

main_wd <- getwd()
setwd(main_wd)
input_dir <- "Analysis/Inputs/1_Raw_data/Mapping_raw_data"
output_dir <- "Analysis/Outputs"
# insert NAS password file
NAS_pass_dir <- ("Analysis/Inputs/1_Raw_data/passNAS.txt") ### Insert the name of your NASS password file
print(NAS_pass_dir)
NAS_pass_path <- readline("insert a NASS pasword file path without the quotes: ")

if (NAS_pass_path == "") {
  NAS_pass_path <- NAS_pass_dir
}

NAS_pass_path

# setwd(output_dir)
# getwd()

###########################
####### 1. Functions ######
###########################

###########################
####### 2. Load data ######
###########################
full_path_df <- read.table(file.path(input_dir,"list_of_fastq_paths/fastq_files_paths.txt"))
folders_name <- read.table(file.path(input_dir,"list_of_fastq_paths/folders_names.txt"))
raw_fastq_names <- read.table(file.path(input_dir,"list_of_fastq_paths/raw_fastq_names.txt"))
tx2gene <- read.table(file.path(input_dir,"tx2gene.txt"))

### add a column just with the file name without the extension fastq.gz
ext_rm <- nchar(".fastq.gz")
colnames(raw_fastq_names) <- "full_names"
raw_fastq_names$fastq_wout_ext <- substr(raw_fastq_names$full_names,1,nchar(raw_fastq_names$full_names)-ext_rm)
length(unique(raw_fastq_names$fastq_wout_ext))

### Remove "19_01_ADN_F1C2r_1.fastq" (row 53) from the list (it is not a fastq sample)
full_path_df <- data.frame(full_path_df[-53,])
raw_fastq_names <- raw_fastq_names[-53,]
### list of all fastqc files
my_files <- raw_fastq_names[1]

####################
### QC Mapping ####
#################### 

### Create directories
# output_dir <-file.path(getwd(),"Analysis/Outputs")
reports_dir <- "2_Reports"
{
  dir1 <- "1_trimm_report"
  dir2 <-"1_trimm_files"
  dir3 <-"1_trimm_fastqc"  
  dir4 <-"1_trimm_outputs"  
  dir5 <- "3_multiqc_reports"
  dir6 <- "2_kallisto_output"
  dir7 <- "1_trimm_zip"
  
  trimm_report <- file.path(output_dir,reports_dir,dir1)
  trimm_files <- file.path(output_dir,reports_dir,dir2)
  trimm_fastqc <- file.path(output_dir,reports_dir,dir3)
  trimm_outputs <- file.path(output_dir,reports_dir,dir4)
  multiqc_reports <- file.path(output_dir,reports_dir,dir5)
  kallisto <- file.path(output_dir,reports_dir,dir6)
  trimm_zip <- file.path(output_dir,reports_dir,dir7)
}


folder_list <- list(dir1,dir2,dir3,dir4,dir5,dir6,dir7)
sapply(folder_list, FUN=function(x) dir.create(file.path(output_dir,reports_dir,x),showWarnings = F,recursive=T))

### Create file where I`m going to download the raw fastq
raw_fastq_folder <- file.path(output_dir,"2_Reports/raw_fastqc")
dir.create(path=raw_fastq_folder,showWarnings = F)

### Path to the gtf and index files
gtf_file <- file.path(input_dir,"genomic.gtf")
index <- file.path(output_dir,"1_Kallisto/index/index_ncbi_ASM19595v2.idx")

### Inputs to Trim Galore
q_threshold=20
max_n=1
stringency_threshold=5
length_threshold=20
num_cores <- 7
# i="raw_fastqc_C1-C2_rep1"
# j=1
my_files

start_time <- Sys.time()
# NAS_password <- "Sanzlab2021*"

# for(i in (c(57:58))){
for(i in c(1:(dim(my_files)[1]/2))){
  
  remote_file1 <- paste0("/volume1/Raw_data/Colab/",full_path_df[((2*i)-1),])
  remote_file2 <- paste0("/volume1/Raw_data/Colab/",full_path_df[(2*i),])
  local_folder <- raw_fastq_folder
  
  ### Type your credentials to copy from NAS
  ssh_command1 <- paste0("sshpass -f ",NAS_pass_dir," scp -O -P 28 -rp JCardenas@155.210.137.80:",remote_file1," ",local_folder)
  ssh_command2 <- paste0("sshpass -f ",NAS_pass_dir," scp -O -P 28 -rp JCardenas@155.210.137.80:",remote_file2," ",local_folder)

  system(ssh_command1)
  system(ssh_command2)

  sample_name1 <- raw_fastq_names[((2*i)-1),2]
  sample_name2 <- raw_fastq_names[(2*i),2]
  dir <- substr(raw_fastq_names[((2*i)-1),2],1,nchar(raw_fastq_names[((2*i)-1),2])-2)
  sample_name <- dir
  output <- file.path(output_dir,"2_Reports/raw_fastqc_output",dir)
  dir.create(output,recursive = T,showWarnings = F)
  
  ### FASTQC pre trimming
  print(my_files[((2*i)-1),])
  fastqc_1 <- file.path(raw_fastq_folder,my_files[((2*i)-1),])
  print(fastqc_1)
    
  print(my_files[(2*i),])
  fastqc_2 <- file.path(raw_fastq_folder,my_files[(2*i),])
  print(fastqc_2)
    
  command=paste0("fastqc ",fastqc_1," ",fastqc_2," -t 8 -o ",output)
  system(command)
    
  #######################
  ####  TRIMM GALORE ####
  ####################### 
  forward=fastqc_1
  reverse=fastqc_2
    
  forward_name = substr(my_files[((2*i)-1),],1,(nchar(my_files[((2*i)-1),])-9))
  reverse_name = substr(my_files[(2*i),],1,(nchar(my_files[(2*i),])-9))
    
  print(forward)
  print(reverse)
   
  #" -o ",trimm_outputs,
    
  command=paste0("trim_galore -q ",q_threshold," -stringency ",stringency_threshold," --length ",length_threshold," --max_n ",max_n," --cores ",num_cores," --trim-n -paired --fastqc ",fastqc_1," ",fastqc_2)
  system(command)
    
  ### trim_galore trimming reports
  dir.create(file.path(trimm_report,sample_name),recursive=TRUE,showWarnings = F)

  mv_report_command_1=paste0("mv ",forward_name,".fastq.gz_trimming_report.txt ",trimm_report,"/",sample_name,"/",forward_name,".fastq.gz_trimming_report.txt")
  mv_report_command_2=paste0("mv ",reverse_name,".fastq.gz_trimming_report.txt ",trimm_report,"/",sample_name,"/",reverse_name,".fastq.gz_trimming_report.txt")
    
  system(mv_report_command_1)
  system(mv_report_command_2)
    
  ### trimmed files
  dir.create(file.path(trimm_files,sample_name),recursive=TRUE,showWarnings = F)
    
  mv_trimmed_file_command_1=paste0("mv ",forward_name,"_val_1.fq.gz ",trimm_files,"/",sample_name,"/",forward_name,"_val_1.fq.gz")
  mv_trimmed_file_command_2=paste0("mv ",reverse_name,"_val_2.fq.gz ",trimm_files,"/",sample_name,"/",reverse_name,"_val_2.fq.gz")
    
  system(mv_trimmed_file_command_1)
  system(mv_trimmed_file_command_2)
    
  ###  fastQC reports on trimmed data.
  dir.create(file.path(trimm_fastqc,sample_name),recursive=TRUE,showWarnings = F)
    
  mv_trimmed_file_fastqc_report_command_1=paste0("mv ",forward_name,"_val_1_fastqc.html ",trimm_fastqc,"/",sample_name,"/",forward_name,"_val_1_fastqc.html")
  mv_trimmed_file_fastqc_report_command_2=paste0("mv ",reverse_name,"_val_2_fastqc.html ",trimm_fastqc,"/",sample_name,"/",reverse_name,"_val_2_fastqc.html")
    
  system(mv_trimmed_file_fastqc_report_command_1)
  system(mv_trimmed_file_fastqc_report_command_2)
    
  ###  ZIP data.
  dir.create(file.path(trimm_zip,sample_name),recursive=TRUE,showWarnings = F)
    
  mv_zip_file_fastqc_report_command_1=paste0("mv ",forward_name,"_val_1_fastqc.zip ",trimm_zip,"/",sample_name,"/",forward_name,"_val_1_fastqc.zip")
  mv_zip_file_fastqc_report_command_2=paste0("mv ",reverse_name,"_val_2_fastqc.zip ",trimm_zip,"/",sample_name,"/",reverse_name,"_val_2_fastqc.zip")
    
  system(mv_zip_file_fastqc_report_command_1)
  system(mv_zip_file_fastqc_report_command_2)
  
  ### Remove the raw data used
  rm_fastq_command <- paste0("rm ", fastqc_1," ",fastqc_2)
  system(rm_fastq_command)
  
  ####################
  #### Kallisto  #####
  #################### 
      
  file <- list.files(file.path(trimm_files,sample_name))
  ### Select the trim files for Kallisto   
  print(file[1])
  file1=file[1]
      
  print(file[2])
  file2=file[2]
      
  name = substr(file[2],1,(nchar(file[2])-14))
      
  kallisto_abundance <- file.path(kallisto,name)
  dir.create(kallisto_abundance,recursive = T,showWarnings = F)
      
  ### command=paste0("kallisto quant -i ",index," -t 12 ",trimm_files,"/",i,"/",file1," ",trimm_files,"/",i,"/",file2," --gtf ",gtf_file," -o ",kallisto)
  command=paste0("kallisto quant -i ",index," -t 12 ",trimm_files,"/",sample_name,"/",file1," ",trimm_files,"/",sample_name,"/",file2," -o ",kallisto_abundance)
      
  system(command)
  
  ### Remove the trim files data used
  fq1 <- paste0(trimm_files,"/",sample_name,"/",file[1])
  fq2 <- paste0(trimm_files,"/",sample_name,"/",file[2])
  
  rm_trimm_command <- paste0("rm ", fq1," ",fq2)
  system(rm_trimm_command)
  
  ####################
  #### tximport  #####
  #################### 
  
  file_tximport <- file.path(kallisto_abundance,"abundance.tsv")
  
  genes_sample <- tximport(file_tximport, type = "kallisto", tx2gene = tx2gene[,c(1,4)], countsFromAbundance = "lengthScaledTPM")
  reads=data.frame(genes_sample$counts)
  colnames(reads) <- name
  
  reads_from_tximport <- paste0("reads_from_tximport/",name)
  dir.create(file.path(output_dir,"2_Reports",reads_from_tximport),showWarnings = F,recursive = T)
  write.table(reads, file.path(output_dir,"2_Reports",reads_from_tximport,paste0(colnames(reads),".txt")))
  
}
end_time <- Sys.time()
print(end_time - start_time)

####################
####  MultiQC   ####
#################### 

output_multiqc <- file.path(output_dir,"2_Reports/raw_fastqc_output")

raw <- file.path(output_multiqc)
html <- file.path(trimm_fastqc)
txt <- file.path(trimm_report)
out <- file.path(multiqc_reports,"multiqc_reports")

command=paste0("multiqc ",raw," --profile-runtime -n ",out,"_raw -o ", out)
command1=paste0("multiqc ",html," --profile-runtime -n ",out,"_trimmed -o ", out)
command2=paste0("multiqc ",txt," --profile-runtime -n ",out,"_txt -o ", out)
command3=paste0("multiqc ",main_wd," --profile-runtime -n ",out,"_txt -o ", out)

system(command)
system(command1)
system(command2)
system(command3)

