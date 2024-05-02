#############################
### 0. Load dependencies. ###
#############################

{
  library(tidyverse)
  library(ggplot2)
  library(reshape2)
  library(dplyr)
  library(gridExtra)
  library(ComplexHeatmap)
  library(dendextend)
  library(umap)
  library(rgl)
  library(RColorBrewer)
  library(circlize)
  library(eulerr)
  library(ggrepel)
  library(cowplot)
}

############################
### 1. Declare functions ### 
############################
dcols=function(x){data.frame(colnames(x))}

###################################################
### 2. Set working directory and create folders ###
###################################################
main_wd <- getwd()
input_dir <- "Analysis/Inputs/2_Processed_data"
output_dir <- "Analysis/Outputs"

#################
##  Load Data  ##
#################
myData <- readRDS(file.path(input_dir,"RDS/Contrasts_stat_0.05.RDS"))

LogFC_df <- read.table(file.path(input_dir,"txt/LogFC_0.05.txt"))
BH_df <- read.table(file.path(input_dir,"txt/BH_0.05.txt"))
lfcSE_df <- read.table(file.path(input_dir,"txt/lfcSE_0.05.txt"))

# # Define the pattern to remove the parenthesis 
# names(myData)
# pattern <- paste(c("[(]", "[)]"), collapse = "|")

### Extract the 4 interactions
# 17 Effect_of_Fe_in_response_to_Growth_arrest_G_D.log2FoldChange_shrunken
# 18   Effect_of_Fe_in_response_to_Growth_arrest_G.log2FoldChange_shrunken
# 19 Effect_of_Fe_in_response_to_Growth_arrest_G_L.log2FoldChange_shrunken
# 20   Effect_of_Fe_in_response_to_Growth_arrest_L.log2FoldChange_shrunken

vec_int <- c(17:20)
LogFC_df_int <- LogFC_df[,vec_int]
BH_df_int <- BH_df[,vec_int]
lfcSE_df_int <- lfcSE_df[,vec_int]

colnames(LogFC_df_int) <- paste0(substring(colnames(LogFC_df_int),1,nchar(colnames(LogFC_df_int))-24),".LogFC")
# colnames(BH_df_int) <- paste0(colnames(int_BH),".BH")
colnames(lfcSE_df_int) <- paste0(substring(colnames(lfcSE_df_int),1,nchar(colnames(lfcSE_df_int))-9))

length(which(rownames(LogFC_df_int)!=rownames(BH_df_int)))
length(which(rownames(LogFC_df_int)!=rownames(lfcSE_df_int)))

### Extract the 8 responses to Growth arrest
vec_iron <- c(9:16)
LogFC_df_iron <- LogFC_df[,vec_iron]
BH_df_iron <- BH_df[,vec_iron]
lfcSE_df_iron <- lfcSE_df[,vec_iron]

colnames(LogFC_df_iron) <- paste0(substring(colnames(LogFC_df_iron),1,nchar(colnames(LogFC_df_iron))-24),".LogFC")
# colnames(BH_df_int) <- paste0(colnames(int_BH),".BH")
colnames(lfcSE_df_iron) <- paste0(substring(colnames(lfcSE_df_iron),1,nchar(colnames(lfcSE_df_iron))-9))

length(which(rownames(LogFC_df_iron)!=rownames(BH_df_iron)))
length(which(rownames(LogFC_df_iron)!=rownames(lfcSE_df_iron)))

###########################
#### 3. Check Coherence ###
###########################


threshold=0.05
threshold_int=0.05

name <- c("C1_C2_Glycerol-Dextrose","C5_C6_Glycerol","C7_C8_Glycerol-LCFA","C11_C12_LCFA")

vec <- c(9:12)

Genes_list_enhanced <- list()

for (i in vec) {
  
  indice <- i-8
  print(paste0("Calculating ",name[indice], " layer"))

  summary=data.frame(
    logFC_with=LogFC_df[,i],
    logFC_without=LogFC_df[,i+4],
    logFC_int=LogFC_df[,i+8],
    BH_with=BH_df[,i],
    BH_without=BH_df[,i+4],
    BH_int=BH_df[,i+8])
  rownames(summary)=rownames(LogFC_df)

  {
  summary$label <- "Background"  
  summary$label[which((summary$BH_with<threshold & summary$BH_without<threshold) & 
                        summary$logFC_with*summary$logFC_without>0)]="Coherent response"
  
  summary$label[which((summary$BH_with<threshold & summary$BH_without<threshold) &
                        summary$logFC_with*summary$logFC_without<0)]="Flipped response"
  
  summary=summary[order(summary$label),]
  
  test=cor.test(summary$logFC_with,summary$logFC_without)
  test$estimate
  # 0.8945709
  test$p.value
  # 0
  
  length(which(summary$label=="Coherent response"))
  # 2348
  length(which(summary$label=="Coherent response"))/(nrow(summary)-length(which(summary$label=="Background")))
  # 0.9779259
  length(which(summary$BH_int<0.05))
  # 1160
  length(which(summary$BH_without<threshold))
  # 2927
  
  #Intercept again
  length(which(summary$BH_with<threshold & summary$BH_without<threshold))
  # 2401
  }
  #############################################
  ### Absolute effect sizes  ###
  #############################################
  {
  threshold_int=0.05
  threshold=0.01
  {
  summary$label_int="Background"
  summary$label_int[which( summary$logFC_with>0 & summary$logFC_without>0)]="UP"
  summary$label_int[which( summary$logFC_with>0 & summary$logFC_without>0 & 
                             summary$BH_int<threshold_int & summary$logFC_int>0)]="UP_More_without_Fe"
  summary$label_int[which( summary$logFC_with>0 & summary$logFC_without>0 & 
                             summary$BH_int<threshold_int & summary$logFC_int<0)]="UP_More_with_Fe"
  
  summary$label_int[which( summary$logFC_with<0 & summary$logFC_without<0)]="DOWN"
  summary$label_int[which( summary$logFC_with<0 & summary$logFC_without<0 & 
                             summary$BH_int<threshold_int & summary$logFC_int>0)]="DOWN_More_with_Fe"
  summary$label_int[which( summary$logFC_with<0 & summary$logFC_without<0 & 
                             summary$BH_int<threshold_int & summary$logFC_int<0)]="DOWN_More_without_Fe"
  
  summary$label_int[which(summary$logFC_with<0 & summary$logFC_without>0 & summary$BH_int<threshold_int)]="FLIPPED_int_down_up_weak"
  summary$label_int[which(summary$logFC_with>0 & summary$logFC_without<0 & summary$BH_int<threshold_int)]="FLIPPED_int_up_down_weak"
  
  summary$label_int[which(summary$BH_with<threshold & summary$BH_without<threshold & summary$logFC_with<0 & summary$logFC_without>0 & summary$BH_int<threshold_int)]="FLIPPED_int_down_up_strong"
  summary$label_int[which(summary$BH_with<threshold & summary$BH_without<threshold & summary$logFC_with>0 & summary$logFC_without<0 & summary$BH_int<threshold_int)]="FLIPPED_int_up_down_strong"
  }
  
  summary_df <- data.frame(summary(factor(summary$label_int)))
  
  up_chunk=summary[which(summary$label_int %in% c("UP_More_with_Fe","UP_More_without_Fe")),]
  up_chunk$int_rebuilt=abs(up_chunk$logFC_without)-abs(up_chunk$logFC_with)
  up_chunk$label_int="Upregulated_genes"
  
  down_chunk=summary[which(summary$label_int %in% c("DOWN_More_with_Fe","DOWN_More_without_Fe")),]
  down_chunk$int_rebuilt=abs(down_chunk$logFC_without)-abs(down_chunk$logFC_with)
  down_chunk$label_int="Downregulated_genes"
  
  # tab=rbind(up_chunk,down_chunk)
 
  ### Save the genes enhanced_upregulated, dampened_upreg, enhanced_downregulated and dampened_downreg
  enhanced_upregulated <- rownames(summary[which(summary$label_int=="UP_More_without_Fe"),])
  damped_upreg <- rownames(summary[which(summary$label_int=="UP_More_with_Fe"),])
  
  enhanced_downregulated <- rownames(summary[which(summary$label_int=="DOWN_More_without_Fe"),])
  damped_downreg <- rownames(summary[which(summary$label_int=="DOWN_More_with_Fe"),])
  
  dir.create(file.path(output_dir,"Enhanced_and_damped_genes",name[indice]),recursive = T,showWarnings = F)
  file_name <- c("enhanced_upregulated.txt","damped_upreg.txt","enhanced_downregulated.txt","damped_downreg.txt")
  dir <- file.path(output_dir,"Enhanced_and_damped_genes",name[indice])
 
  writeLines(enhanced_upregulated, file.path(dir,file_name[1]))
  writeLines(damped_upreg, file.path(dir,file_name[2]))
  writeLines(enhanced_downregulated, file.path(dir,file_name[3]))
  writeLines(damped_downreg, file.path(dir,file_name[4]))
  }
}


