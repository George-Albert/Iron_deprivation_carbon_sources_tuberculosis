#############################
### 0. Load dependencies. ###
#############################

{ library(tidyverse)
  library(limma)
  library(ggplot2)
  library(DESeq2)
  library(edgeR)
  library(dplyr)
  library(locfdr)
  library(xlsx)
  library(writexl)
  library(stringr)
  library(ggrepel)
  library(gtools)
}

############################
### 1. Declare functions ###
############################
dcols=function(x){data.frame(colnames(x))}
# generate_samples_to_compare <- function(metadata, sample_ids) {
### Reads_data <-  data frame of reads
### met_df <- data.frame of the metadata
### samples_to_compare <-  the samples i am going to compare between with
generate_pca_plt <- function(reads_data, met_df, samples_to_compare) {

  # Filter the reads by selected cultures
  filtered_reads <- reads_data[, samples_to_compare]
  met_df <- subset(metadata, Sample_ID %in% samples_to_compare)
  # Do the PCA 2D
  y <- DGEList(counts = filtered_reads)
  y <- calcNormFactors(y)
  design <- model.matrix(~short_setup,data=met_df)
  v=voom(y,design,plot=TRUE)
  exp=v$E
  pca <-prcomp(t(exp),scale=T)
  sum_pca=data.frame(summary(pca)$importance)
  # Create a df with the meta data and the data from PCA
  pca_df <- as.data.frame(pca$x)

  # Bind the metadata and PCA data
  combined_data <- bind_cols(met_df,pca_df)

  dir.create(file.path(output_dir,"002_PCA","PCA_tables"),recursive=TRUE,showWarnings = F)
  write.table(combined_data,file=file.path(output_dir,"002_PCA","PCA_tables",paste0(tab_name,".txt")))

  fill_base <- combined_data$fill[which(!duplicated(combined_data$short_setup))]
  color_base <- combined_data$color[which(!duplicated(combined_data$short_setup))]
  # Generate the plot
  pca_plt <- ggplot(combined_data, aes(x = PC1, y = PC2, fill = short_setup, color = short_setup)) +
    geom_point(size = 3, shape = 21,stroke=1.5) +
    geom_label_repel(label = combined_data$short_setup,color="black",fill=combined_data$fill,
                     max.overlaps = 21, size = 2, nudge_x = 0.1, nudge_y = 0.1,show.legend = F) +
    theme_bw() +
    scale_fill_manual(values = fill_base,name = "Samples",labels=unique(combined_data$short_setup)) +
    scale_color_manual(values = color_base,name = "Samples",labels=unique(combined_data$short_setup)) +
    xlab(paste0("PC1:",100*round(sum_pca$PC1[2],digits=3),"% variance explained"))+
    ylab(paste0("PC2:",100*round(sum_pca$PC2[2],digits=3),"% variance explained"))+theme_bw()+
    theme(
      axis.text.y   = element_text(size=14),
      axis.text.x   = element_text(size=14),
      axis.title.y  = element_text(size=14),
      axis.title.x  = element_text(size=14),
      # panel.background = element_blank(),
      # panel.grid.major = element_blank(),
      # panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      panel.border = element_rect(colour = "black", fill=NA, linewidth=1,linetype="solid"),
      # legend.title=element_blank(),
      # legend.position="none",
      legend.text=element_text(size=14),
      legend.key.size = unit(1, 'lines'))

  pca_plt

  dir.create(file.path(output_dir,"002_PCA","PCA_figures",dir_pca),recursive=TRUE,showWarnings = F)
  pdf(file = file.path(output_dir,"002_PCA","PCA_figures",dir_pca,paste0(pca_name,".pdf")),width=6,height=5)
  print(pca_plt)
  dev.off()

  return(pca_plt)
}


###################################################
### 2. Set working directory and create folders ###
###################################################

main_wd <- getwd()
setwd(main_wd)
input_dir  <- "Inputs/002_Processed_data"
output_dir <- "Outputs"

pca_dir <- file.path(output_dir,"001_Figures_paper")
dir.create(pca_dir)

###########################
####### 3. Load data ######
###########################
metadata     <- read.table(file.path(input_dir,"txt/metadata_32_samples.txt"),check.names = F)
reads        <- read.table(file.path(input_dir,"txt/reads_32_samples.txt"),check.names = F)
feature_data <- read.table(file.path(input_dir,"txt/feature_data_filtered.txt"),check.names = F)
# DE_genes_df<- read.table(file.path(input_dir,"DE_genes_per_cond.txt"),check.names = F)

reads <- reads[rownames(feature_data),]
length(which(rownames(feature_data)!=(rownames(reads))))
length(which(colnames(reads)!=rownames(metadata)))

metadata <- metadata[order(metadata$short_setup),]
cult <- paste0("C",c(1:2,5:6,7:8,11:12))

### Set the colors for the cultures
#C1 <- "lightcoral"
#C2 <- "red2"
#C5 <- "lightgoldenrod"
#C6 <- "yellow2"
#C7 <- "skyblue"
#C8 <- "blue4"
#C11 <- "lawngreen"
#C12 <- "darkgreen

# color_vec <- c(C1 = "#f6a1a1", C2 = "#852121", C5 = "#ffc499", C6 = "#8f4a17",
#                C7 = "#a9dea9", C8 = "#2c682c", C11 = "#8cd9ff", C12 = "#355f91")

color_vec <- c(C1 = "#d3c7d9", C2 = "#5f4290", C5 = "#76c6e7", C6 = "#3c8ec8",
               C7 = "#bbdd93", C8 = "#559d3f", C11 = "#ee9f9b", C12 = "#d1382c")

cultures <- metadata$short_setup

# df_met      <- Dataframe of metadata
# culture_vec <- Vector with the name of the Cultures
# reads_df    <- Dataframe of the reads
# iron        <- YES or NO options
# tab_name    <- Name of the pca tables
# pca_name    <- Name of the PCA figures
# dir_pca     <- folder of PCA (attending to conditions or effects)

### Effect_Iron_on_response_to_growth_arrest
{
vec <- c(1:4)
dir_pca <-"Effect_Iron_on_response_to_growth_arrest"
samples <- data.frame(samples=metadata$Sample_ID)

for (i in vec) {
  # culture_vec <- cult[c(2*i,(2*i)+2)]
  culture_vec <- cult[c(2*i-1,2*i)]
  print((culture_vec))

  patterns <- paste0(culture_vec,"_.*_.*_.*")
  samples_to_compare <- samples$samples[grep(paste(patterns, collapse = "|"),samples$samples)]

  name <- unique(substr(samples_to_compare,1,nchar(samples_to_compare)-2))
  tab_name <- paste0("PCA_",culture_vec[1],"_vs_",culture_vec[2],"_table")
  print(tab_name)
  pca_name <- paste0("PCA_",culture_vec[1],"_vs_",culture_vec[2])

  pca_plt <- generate_pca_plt(reads_data=reads, metadata, samples_to_compare)
  print(pca_plt)
}
}

### Carbon_effects_with_Fe
{
samples <- data.frame(samples=metadata$Sample_ID)
vec <- c(1:5)
dir_pca <- "Carbon_effects_with_Fe"

for (i in vec) {
  if (i %% 2 == 0) {  # verify if is even
    next  # jump to the next iteration
  }
  start <- i
  end <- start + 3
  sequence <- start:end
  culture_vec <- cult[sequence]
  print(culture_vec)
  patterns <- paste0(culture_vec,"_.*_YES_.*")
  samples_to_compare <- samples$samples[grep(paste(patterns, collapse = "|"),samples$samples)]

  name <- unique(substr(samples_to_compare,1,nchar(samples_to_compare)-2))
  tab_name <- paste0("PCA_",culture_vec[1],"_",culture_vec[2],"_vs_",culture_vec[3],"_",culture_vec[4],"_table")
  print(tab_name)
  pca_name <- paste0("PCA_",culture_vec[1],"_",culture_vec[2],"_vs_",culture_vec[3],"_",culture_vec[4])

  pca_plt <- generate_pca_plt(reads_data=reads, metadata, samples_to_compare)
  print(pca_plt)
}
}

### Carbon_effects_without_Fe
{
  samples <- data.frame(samples=metadata$Sample_ID)

  vec <- c(1:5)
  dir_pca <- "Carbon_effects_without_Fe"

  for (i in vec) {
    if (i %% 2 == 0) {  # verify if is even
      next  # jump to the next iteration
    }
    start <- i
    end <- start + 3
    sequence <- start:end
    culture_vec <- cult[sequence]
    print(culture_vec)
    patterns <- paste0(culture_vec,"_.*_NO_.*")
    samples_to_compare <- samples$samples[grep(paste(patterns, collapse = "|"),samples$samples)]

    name <- unique(substr(samples_to_compare,1,nchar(samples_to_compare)-2))
    tab_name <- paste0("PCA_",culture_vec[1],"_",culture_vec[2],"_vs_",culture_vec[3],"_",culture_vec[4],"_table")
    print(tab_name)
    pca_name <- paste0("PCA_",culture_vec[1],"_",culture_vec[2],"_vs_",culture_vec[3],"_",culture_vec[4])

    pca_plt <- generate_pca_plt(reads_data=reads, metadata, samples_to_compare)
    print(pca_plt)
  }
}
### Carbon_sources_vs_Phase_effects
{
  samples <- data.frame(samples=metadata$Sample_ID)
  dir_pca <- "Carbon_sources_vs_Phase_effects"

  culture_vec <- cult[c(3:6)]
  print(culture_vec)
  patterns <- paste0(culture_vec,"_.*_NO_.*")
  samples_to_compare <- samples$samples[grep(paste(patterns, collapse = "|"),samples$samples)]

  name <- unique(substr(samples_to_compare,1,nchar(samples_to_compare)-2))
  tab_name <- paste0("PCA_",culture_vec[1],"_",culture_vec[2],"_vs_",culture_vec[3],"_",culture_vec[4],"_NO_table")
  print(tab_name)
  pca_name <- paste0("PCA_",culture_vec[1],"_",culture_vec[2],"_vs_",culture_vec[3],"_",culture_vec[4],"_NO")

  pca_plt <- generate_pca_plt(reads_data=reads, metadata, samples_to_compare)
  print(pca_plt)
}

### EXP_phase_vs_Carbon_effects
{
vec <- c(1:3)
### dir_pca <- "STAT_phase_vs_Carbon_effects"
dir_pca <-"EXP_phase_vs_Carbon_effects"
for (i in vec) {
  # culture_vec <- cult[c(2*i,(2*i)+2)]
  culture_vec <- cult[c(2*i-1,(2*i-1)+2)]
  print((culture_vec))
  patterns <- paste0(culture_vec,"_.*_.*_.*")
  samples_to_compare <- samples$samples[grep(paste(patterns, collapse = "|"),samples$samples)]

  name <- unique(substr(samples_to_compare,1,nchar(samples_to_compare)-2))
  tab_name <- paste0("PCA_",culture_vec[1],"_vs_",culture_vec[2],"_table")
  print(tab_name)
  pca_name <- paste0("PCA_",culture_vec[1],"_vs_",culture_vec[2])

  pca_plt <- generate_pca_plt(reads_data=reads, metadata, samples_to_compare)
  print(pca_plt)
}
}
### STAT_phase_vs_Carbon_effects
{
  vec <- c(1:3)
  dir_pca <- "STAT_phase_vs_Carbon_effects"
  for (i in vec) {
    culture_vec <- cult[c(2*i,(2*i)+2)]
    print((culture_vec))
    patterns <- paste0(culture_vec,"_.*_.*_.*")
    samples_to_compare <- samples$samples[grep(paste(patterns, collapse = "|"),samples$samples)]

    name <- unique(substr(samples_to_compare,1,nchar(samples_to_compare)-2))
    tab_name <- paste0("PCA_",culture_vec[1],"_vs_",culture_vec[2],"_table")
    print(tab_name)
    pca_name <- paste0("PCA_",culture_vec[1],"_vs_",culture_vec[2])

    pca_plt <- generate_pca_plt(reads_data=reads, metadata, samples_to_compare)
    print(pca_plt)
  }
}


###Iron_effects_vs Carbon_sources
{
  samples <- data.frame(samples=metadata$Sample_ID)
  dir_pca <- "Iron_effects_vs_Carbon_sources"

  culture_vec<- cult[c(2,8)]
  # culture_vec <- cult[c(1:2,7:8)]
  print(culture_vec)
  # patterns <- c(paste0(culture_vec[c(1,2)],"_.*_YES_.*"),paste0(culture_vec[c(3,4)],"_.*_NO_.*"))
  patterns <- paste0(culture_vec[c(1,2)],"_.*_.*_.*")
  samples_to_compare <- samples$samples[grep(paste(patterns, collapse = "|"),samples$samples)]

  name <- unique(substr(samples_to_compare,1,nchar(samples_to_compare)-2))
  # tab_name <- paste0("PCA_",culture_vec[1],"_",culture_vec[2],"_YES_vs_",culture_vec[3],"_",culture_vec[4],"_NO_table")
  # print(tab_name)
  # pca_name <- paste0("PCA_",culture_vec[1],"_",culture_vec[2],"_YES_vs_",culture_vec[3],"_",culture_vec[4],"_NO")
  tab_name <- paste0("PCA_",culture_vec[1],"_vs_",culture_vec[2],"_table")
  print(tab_name)
  pca_name <- paste0("PCA_",culture_vec[1],"_vs_",culture_vec[2])

  pca_plt <- generate_pca_plt(reads_data=reads, metadata, samples_to_compare)
  print(pca_plt)
}
{
  culture_vec <- cult[c(1:2,7:8)]
  print(culture_vec)
  patterns <- c(paste0(culture_vec[c(1,2)],"_.*_YES_.*"),paste0(culture_vec[c(3,4)],"_.*_NO_.*"))
  samples_to_compare <- samples$samples[grep(paste(patterns, collapse = "|"),samples$samples)]

  name <- unique(substr(samples_to_compare,1,nchar(samples_to_compare)-2))
  tab_name <- paste0("PCA_",culture_vec[1],"_",culture_vec[2],"_YES_vs_",culture_vec[3],"_",culture_vec[4],"_NO_table")
  print(tab_name)
  pca_name <- paste0("PCA_",culture_vec[1],"_",culture_vec[2],"_YES_vs_",culture_vec[3],"_",culture_vec[4],"_NO")

  pca_plt <- generate_pca_plt(reads_data=reads, metadata, samples_to_compare)
  print(pca_plt)
}

### Exp effects
{
  # dir_pca <-"Stat effects"
  dir_pca <-"Exp effects"

  # culture_vec <- cult[c(2,4,6,8)]
  culture_vec <- cult[c(1,3,5,7)]

  print((culture_vec))

  patterns <- paste0(culture_vec,"_.*_.*_.*")
  samples_to_compare <- samples$samples[grep(paste(patterns, collapse = "|"),samples$samples)]

  name <- unique(substr(samples_to_compare,1,nchar(samples_to_compare)-2))
  tab_name <- c("PCA_all_exp_table")
  print(tab_name)
  pca_name <- c("PCA_all_exp")

  pca_plt <- generate_pca_plt(reads_data=reads, metadata, samples_to_compare)
  print(pca_plt)
}
### Stat effects
{
  dir_pca <-"Stat effects"

  culture_vec <- cult[c(2,4,6,8)]

  print((culture_vec))

  patterns <- paste0(culture_vec,"_.*_.*_.*")
  samples_to_compare <- samples$samples[grep(paste(patterns, collapse = "|"),samples$samples)]

  name <- unique(substr(samples_to_compare,1,nchar(samples_to_compare)-2))
  tab_name <- c("PCA_all_stat_table")
  print(tab_name)
  pca_name <- c("PCA_all_stat")

  pca_plt <- generate_pca_plt(reads_data=reads, metadata, samples_to_compare)
  print(pca_plt)
}


##################################################################
##                  All Cultures PCA wo labels                  ##
##################################################################

{
  samples <- data.frame(samples=metadata$Sample_ID)
  dir_pca <-"All_Cultures"
  # reads_data=reads
  culture_vec <- cult
  print((culture_vec))

  patterns <- paste0(culture_vec,"_.*_.*_.*")
  samples_to_compare <- samples$samples[grep(paste(patterns, collapse = "|"),samples$samples)]
  name <- unique(substr(samples_to_compare,1,nchar(samples_to_compare)-2))
  tab_name <- c("PCA_all_C12_table")
  print(tab_name)
  pca_name <- c("PCA_C12_all")

  # Filter the reads by selected cultures
  filtered_reads <- reads[, samples_to_compare]
  met_df <- subset(metadata, Sample_ID %in% samples_to_compare)
  # Do the PCA 2D
  y <- DGEList(counts = filtered_reads)
  y <- calcNormFactors(y)
  design <- model.matrix(~short_setup,data=met_df)
  v=voom(y,design,plot=TRUE)
  exp=v$E
  exp <- exp[,mixedsort(colnames(exp))]
  pca <-prcomp(t(exp),scale=T)
  sum_pca=data.frame(summary(pca)$importance)
  # Create a df with the meta data and the data from PCA
  pca_df <- as.data.frame(pca$x)

  # Bind the metadata and PCA data
  combined_data <- bind_cols(met_df,pca_df)

  dir.create(file.path(output_dir,"002_PCA","PCA_tables"),recursive=TRUE,showWarnings = F)
  write.table(combined_data,file=file.path(output_dir,"002_PCA","PCA_tables",paste0(tab_name,".txt")))

  fill_base <- combined_data$fill[which(!duplicated(combined_data$short_setup))]
  color_base <- combined_data$color[which(!duplicated(combined_data$short_setup))]
  # Generate the plot
  pca_plt <- ggplot(combined_data, aes(x = PC1, y = PC2, fill = mixedsort(short_setup), color = mixedsort(short_setup))) +
    geom_point(size = 3, shape = 21,stroke=1.5) +
    # geom_label_repel(label = combined_data$short_setup,color="black",fill=combined_data$fill,
    #                  max.overlaps = 21, size = 2, nudge_x = 0.1, nudge_y = 0.1,show.legend = F) +
    theme_bw() +
    scale_fill_manual(values = fill_base,name = "Samples",labels=unique(combined_data$short_setup)) +
    scale_color_manual(values = color_base,name = "Samples",labels=unique(combined_data$short_setup)) +
    xlab(paste0("PC1:",100*round(sum_pca$PC1[2],digits=3),"% variance explained"))+
    ylab(paste0("PC2:",100*round(sum_pca$PC2[2],digits=3),"% variance explained"))+theme_bw()+
    theme(
      axis.text.y   = element_text(size=14),
      axis.text.x   = element_text(size=14),
      axis.title.y  = element_text(size=14),
      axis.title.x  = element_text(size=14),
      # panel.background = element_blank(),
      # panel.grid.major = element_blank(),
      # panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      panel.border = element_rect(colour = "black", fill=NA, linewidth=1,linetype="solid"),
      # legend.title=element_blank(),
      # legend.position="none",
      legend.text=element_text(size=14),
      legend.key.size = unit(1, 'lines'))

    pca_plt

  dir.create(file.path(output_dir,"002_PCA","PCA_figures",dir_pca),recursive=TRUE,showWarnings = F)
  pdf(file = file.path(output_dir,"002_PCA","PCA_figures",dir_pca,paste0(pca_name,".pdf")),width=6,height=5)
  print(pca_plt)
  dev.off()
    ### Save in 4_Figures_paper
  dir.create(file.path(output_dir,"001_Figures_paper"),recursive=TRUE,showWarnings = F)
  pdf(file = file.path(output_dir,"001_Figures_paper","1_Figure_2A_PCA.pdf"),width=6,height=5)
  print(pca_plt)
  dev.off()


  ##################################################################
  ##                  All Cultures PCA with labels                ##
  ##################################################################
  dir_pca <-"All_Cultures_with_labels"

  pca_plt <- generate_pca_plt(reads_data=reads, metadata, samples_to_compare)
  print(pca_plt)
  metadata <- metadata[order(metadata$short_setup),]


  ##==========================
  ##  Correlation patterns  ==
  ##==========================
  combined_data <- combined_data %>% add_column(culture_number=rep(0,dim(samples)[1]))
  combined_data$culture_number[which(combined_data$Growth=="STAT")]=1

  combined_data <- combined_data %>% add_column(iron_number=rep(0,dim(samples)[1]))
  combined_data$iron_number[which(combined_data$Iron == "NO")]=1

  combined_data$growth_phase <- factor(
    combined_data$culture_number,
    levels = c(0, 1),
    labels = c("Exponential", "Stationary")
  )

  t.test(PC1 ~ growth_phase, data = combined_data)
  # Welch Two Sample t-test
  #
  # data:  PC1 by growth_phase
  # t = -3.7591, df = 28.872, p-value = 0.0007703
  # alternative hypothesis: true difference in means between group Exponential and group Stationary is not equal to 0
  # 95 percent confidence interval:
  #   -73.00976 -21.55129
  # sample estimates:
  #   mean in group Exponential  mean in group Stationary
  # -23.64026                  23.64026

  combined_data$LCFA_status <- ifelse(
    combined_data$LCFA == "YES",
    "LCFA",
    "non_LCFA"
  )

  combined_data$LCFA_status <- factor(
    combined_data$LCFA_status,
    levels = c("non_LCFA", "LCFA")
  )

  t.test(PC2 ~ LCFA_status, data = combined_data)
  # Welch Two Sample t-test
  #
  # data:  PC2 by LCFA_status
  # t = -3.8542, df = 23.348, p-value = 0.0007907
  # alternative hypothesis: true difference in means between group non_LCFA and group LCFA is not equal to 0
  # 95 percent confidence interval:
  #   -44.61411 -13.46625
  # sample estimates:
  #   mean in group non_LCFA     mean in group LCFA
  # -14.52009               14.52009
  ### Print the stats
  dir.create(file.path(output_dir, "002_PCA","Stats"),recursive=TRUE)

  sink(file.path(output_dir, "002_PCA","Stats","PC1_vs_culture_correlation_all_samples.txt"))
  print("cor.test(combined_data$PC1,combined_data$culture_number)")
  cor.test(combined_data$PC1,combined_data$culture_number)
  # cor.test(combined_data$PC1[which(combined_data$Culture==c("C1","C2"))],combined_data$culture_number[which(combined_data$Culture==c("C1","C2"))])
  sink()

  #t = 13.119, df = 6, p-value = 1.21e-05
  #alternative hypothesis: true correlation is not equal to 0
  #95 percent confidence interval:
  # 0.9057638 0.9970361
  #sample estimates:
  # cor
  #0.9830122

  sink(file.path(output_dir, "002_PCA","Stats","PC2_vs_iron_correlation_stat_samples.txt"))
  print("cor.test(datos$PC2[which(datos$Culture==\"C6\")],datos$iron_number[which(datos$Culture==\"C6\")])")
  # cor.test(combined_data$PC2[which(combined_data$Growth=="STAT")],combined_data$iron_number[which(combined_data$Growth=="STAT")])
  cor.test(combined_data$PC2[which(combined_data$Culture=="C8")],combined_data$iron_number[which(combined_data$Culture=="C8")])
  sink()

}

#################################################################
##                   Dendrogram All Cultures                   ##
#################################################################

library(dendextend)
### hclust of the samples

reads_data=reads
filtered_reads <- reads_data[, samples_to_compare]
met_df <- subset(metadata, Sample_ID %in% samples_to_compare)
y <- DGEList(counts = filtered_reads)
y <- calcNormFactors(y)
design <- model.matrix(~short_setup,data=met_df)
v=voom(y,design,plot=TRUE)
exp=v$E

hclust_matrix <- exp
dist_matrix <- dist(hclust_matrix)
h_clust <- hclust(dist_matrix,method="ward.D2")
# columns_name_raw <- c("C1_Fe_NO", "C1_Fe_NO", "C1_Fe_YES",  "C1_Fe_YES", "C2_Fe_NO", "C2_Fe_NO",
#                   "C2_Fe_YES", "C2_Fe_YES",  "C5_Fe_NO", "C5_Fe_NO", "C5_Fe_YES","C5_Fe_YES",
#                   "C6_Fe_NO", "C6_Fe_NO", "C6_Fe_YES",  "C6_Fe_YES", "C7_Fe_NO", "C7_Fe_NO",
#                   "C7_Fe_YES", "C7_Fe_YES", "C8_Fe_NO", "C8_Fe_NO", "C8_Fe_YES","C8_Fe_YES",
#                   "C11_Fe_NO", "C11_Fe_NO", "C11_Fe_YES","C11_Fe_YES", "C12_Fe_NO", "C12_Fe_NO",
#                   "C12_Fe_YES", "C12_Fe_YES")

columns_name_raw <- gsub("_[0-9]+$", "", colnames(hclust_matrix))

colnames(hclust_matrix) <- columns_name_raw
col_dend <- hclust(dist(t(hclust_matrix)))

# columns_name <- colnames(hclust_matrix)
column_title = "Samples"
clust_name <- "Expression levels"

# plot(as.dendrogram(h_clust))
dend <-as.dendrogram(col_dend)
# columns_name <- columns_name_raw[order.dendrogram(dend)]

# Extract the prefix from the sample names
prefix <- sub("_.*", "", columns_name_raw)
# Assign colors based on the prefix
colors_assigned <- color_vec[prefix]
# Check the assigned colors
head(data.frame(Sample = columns_name_raw, Prefix = prefix, Color = colors_assigned), 32)

color_vec_rep <- colors_assigned
names(color_vec_rep) <- columns_name_raw
names(color_vec_rep)

# Extract the condition from the sample names attending to Fe status
condition <- ifelse(grepl("NO", columns_name_raw), "NO",
                    ifelse(grepl("YES", columns_name_raw), "YES", NA))
# Assign shapes based on the condition: 19 for "YES" and 1 for "NO"
pch_vec <- ifelse(grepl("YES", columns_name_raw), 19,
                  ifelse(grepl("NO", columns_name_raw), 1, NA))
names(pch_vec) <- columns_name_raw

# Check the assigned shapes
head(data.frame(Sample = columns_name_raw, Condition = condition, PCH = pch_vec))

# Reorder color_vec_rep and pch_vec according to the dendrogram order
order_dend    <- order.dendrogram(dend)
color_vec_rep <- color_vec_rep[order_dend]
pch_vec       <- pch_vec[order_dend]
# pch_vec <- c(19,19,1,1,1,19,1,19,1,19,1,19,1,19,1,19,rep(1,6),rep(19,8),1,1)

dir_pca <- "All_Cultures"
pdf(file = file.path(output_dir,"002_PCA","PCA_figures",dir_pca,paste0("dendrogram_all_samples.pdf")),width=8,height=5)

dend %>%
# set("labels", columns_name) %>%
  set("leaves_pch", pch_vec) %>%
  set("labels_col", F) %>%
  set("leaves_cex", 2) %>%  # node point size
  set("leaves_col", color_vec_rep) %>% # node point color
  set("branches_k_color",value = c("grey", "grey30"), k = 2)%>%
  set("branches_lwd", c(2,1,2)) %>%
  # rect.dendrogram(k=2, border = 8, lty = 5, lwd = 2)%>%
  plot(main = "Samples")

dev.off()

dir.create(file.path(output_dir,"001_Figures_paper"),recursive=TRUE,showWarnings = F)
pdf(file = file.path(output_dir,"001_Figures_paper","2_Figure_2B_Dendrogram.pdf"),width=8,height=5)
dend %>%
  # set("labels", columns_name) %>%
  set("leaves_pch", pch_vec) %>%
  set("labels_col", F) %>%
  set("leaves_cex", 2) %>%  # node point size
  set("leaves_col", color_vec_rep) %>% # node point color
  set("branches_k_color",value = c("grey", "grey30"), k = 2)%>%
  set("branches_lwd", c(2,1,2)) %>%
  # rect.dendrogram(k=2, border = 8, lty = 5, lwd = 2)%>%
  plot(main = "Samples")

dev.off()

first_cluster <- plot(dend, xlim = c(1, 16), ylim = c(0,65))
second_cluster <- plot(dend, xlim = c(17, 32), ylim = c(0,100))
