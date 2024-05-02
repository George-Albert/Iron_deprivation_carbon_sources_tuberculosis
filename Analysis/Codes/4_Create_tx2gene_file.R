
###########################################################################
###########################################################################
###                                                                     ###
###                           TXIMPORT SCRIPT                           ###
###                                                                     ###
###########################################################################
###########################################################################


##################################################################
##                         Dependencies                         ##
##################################################################
library(tximportData)
library(tximport)
library(Biostrings)

input_dir <- "Analysis/Inputs/1_Raw_data/Mapping_raw_data"
output_dir <- "Analysis/Outputs"
dir <- system.file("extdata", package = "tximportData")
list.files(dir)

#################################################################
##                     List and load files                     ##
#################################################################
feature_data <- read.table(file.path(input_dir,"feature_data_genes.txt"))
fasta_path <- file.path(input_dir,"gene.fna")

# Read sequences from a FASTA file located at 'fasta_path' 
# and store them in a DNAStringSet object named 'fasta_df'
fasta_df <- readDNAStringSet(fasta_path)
trans_name <- names(fasta_df)
trans_name <- str_split_fixed(trans_name," ",n=3)
trans_name <- data.frame(trans_name)
geneid <- trans_name[,3]
trans_name1 <- str_split_fixed(geneid,"\\] \\[",n=3)
sequence = paste(fasta_df)

gene_from_feature <- data.frame(Gene_id=feature_data$GeneID,Locus_tag=feature_data$locus_tag)
tx2gene <- data.frame(transcript_id=trans_name[,1],Gene_id=trans_name1[,2],symbol=trans_name[,2])
tx2gene$Gene_id <- as.numeric(str_remove(tx2gene$Gene_id,"GeneID="))

gene_from_feature <- gene_from_feature[order(gene_from_feature$Gene_id),]
rownames(gene_from_feature) <- gene_from_feature$Gene_id
rownames(tx2gene) <- tx2gene$Gene_id

tx2gene <- tx2gene[rownames(gene_from_feature),]

length(which(rownames(gene_from_feature)!=rownames(tx2gene)))

tx2gene$Locus_tag <- gene_from_feature$Locus_tag
tx2gene <- tx2gene[order(tx2gene$Locus_tag),]

write.table(tx2gene,file.path(input_dir,"tx2gene.txt"))






