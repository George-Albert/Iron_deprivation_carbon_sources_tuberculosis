
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
input_dir  <- "Inputs/001_Raw_data"
output_dir <- "Outputs"

###########################
####### 1. Functions ######
###########################

###########################
####### 2. Load data ######
###########################
### Load ncRNA reads news
reads_ncRNA <- read.table(file.path(input_dir,"reads_bowtie.txt"))

### Load reads news
reads_new   <- read.table(file.path(input_dir,"reads_new.txt"))

#############################################################
### Correlate the new reads with the reads of the company ###
#############################################################

reads_ncRNA <- reads_ncRNA[,colnames(reads_new)]

# Check if the column names match
length(which(colnames(reads_ncRNA)!= colnames(reads_new)))

#Remove the "X" prefix from column names in reads_new
column_names <- colnames(reads_new)
column_names <- gsub("^X","",column_names)
colnames(reads_new) <- column_names
# Remove the "X" prefix from column names in reads_ncRNA
column_names <- colnames(reads_ncRNA)
column_names <- gsub("^X","",column_names)
colnames(reads_ncRNA) <- column_names

# Check if the column names match
length(which(colnames(reads_ncRNA)!= colnames(reads_new)))

# Extract the gene names
ncgenes <- rownames(reads_ncRNA)
genes <- rownames(reads_new)

# Find ncRNA genes that are also in the new reads
idx_ncRNA <- which(ncgenes %in% genes)
reads_ncRNA_filtered <- reads_ncRNA[idx_ncRNA,]
reads_ncRNA_filtered <- reads_ncRNA_filtered[order(rownames(reads_ncRNA_filtered)),]

#Check if the row names match
length(which(rownames(reads_ncRNA_filtered)!= rownames(reads_new)))

zero_sd_rows <- sapply(1:nrow(reads_new), function(i) {
      sd1 <- sd(as.numeric(reads_new[i, ]), na.rm = TRUE)
      sd2 <- sd(as.numeric(reads_ncRNA_filtered[i, ]), na.rm = TRUE)
      return(sd1 == 0 || sd2 == 0)
  })
reads_sd <- reads_new[which(zero_sd_rows),]
# Correlate the new reads with the ncRNA reads

corr_matrix <- cor(reads_new,reads_ncRNA_filtered,method = "pearson")
corr_matrix_row <- sapply(1:nrow(as.matrix(reads_new)), function(i) cor(as.matrix(reads_new)[i,],
                                                                        as.matrix(reads_ncRNA_filtered)[i,],method = "pearson"))
corr_matrix <- data.frame(corr_matrix)
diagonal <- corr_matrix[row(corr_matrix)==col(corr_matrix)]
diagonal <- data.frame(corr_value=diagonal)
rownames(diagonal) <- colnames(reads_new)

corr_matrix_row <- data.frame(corr_value=corr_matrix_row)
rownames(corr_matrix_row) <- rownames(reads_new)
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

dir.create(file.path(output_dir,"001_Correlation_reads_ncRNA_vs_reads_new"),recursive = T,showWarnings = F)
pdf(file = file.path(output_dir,"001_Correlation_reads_ncRNA_vs_reads_new","dens_plot_correlation.pdf"),width=9,height=5)

density(dens2,"Density Plot by sample")
density(dens1,"Density Plot by gene")

dev.off()

# Save the correlation matrices
dir.create(file.path(input_dir,"001_Correlation_reads_ncRNA_vs_reads_new"),recursive = T,showWarnings = F)
write.table(corr_matrix_row, file = file.path(input_dir,"001_Correlation_reads_ncRNA_vs_reads_new","corr_matrix_row.txt"))
write.table(diagonal, file = file.path(input_dir,"001_Correlation_reads_ncRNA_vs_reads_new","diagonal.txt"))

# Save summary statistics
stats_names <- c("max","quantile_75","mean","median", "quantile_25","min", "sd")
corr1 <- corr_matrix_row$corr_value
corr2 <- diagonal$corr_value

stats_corr1 <- c(
  max(corr1, na.rm = TRUE),
  quantile(corr1, 0.75, na.rm = TRUE),
  mean(corr1, na.rm = TRUE),
  median(corr1, na.rm = TRUE),
  quantile(corr1, 0.25, na.rm = TRUE),
  min(corr1, na.rm = TRUE),
  sd(corr1, na.rm = TRUE)
)

stats_corr2 <- c(
  max(corr2, na.rm = TRUE),
  quantile(corr2, 0.75, na.rm = TRUE),
  mean(corr2, na.rm = TRUE),
  median(corr2, na.rm = TRUE),
  quantile(corr2, 0.25, na.rm = TRUE),
  min(corr2, na.rm = TRUE),
  sd(corr2, na.rm = TRUE)
)

# Crear el data frame final
summary_stats_df <- data.frame(
  stat           = stats_names,
  corr_by_gene   = round(stats_corr1, 3),
  corr_by_sample = round(stats_corr2, 3)
)

print(summary_stats_df)
# Save the summary statistics to a text file

write.table(summary_stats_df, file = file.path(input_dir,"001_Correlation_reads_ncRNA_vs_reads_new","summary_statistics.txt"),
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
