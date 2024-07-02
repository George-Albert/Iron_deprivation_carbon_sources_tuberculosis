#############################
### 0. Load dependencies. ###
#############################

{
  library(tidyverse)
  library(limma)
  library(DESeq2)
  library(edgeR)
  library(dplyr)
  library(ggrepel)
  library(igraph)
  library(xlsx)
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(circlize)
  library(ggsci)
  library(ggdendro)
  library(dendextend)
  library(dendsort)
}

############################
### 1. Declare functions ### 
############################
dcols=function(x){data.frame(colnames(x))}
calc_ht_size = function(ht, unit = "inch") {
  pdf(NULL)
  ht = draw(ht)
  w = ComplexHeatmap:::width(ht)
  w = convertX(w, unit, valueOnly = TRUE)
  h = ComplexHeatmap:::height(ht)
  h = convertY(h, unit, valueOnly = TRUE)
  dev.off()
  
  c(w, h)
}
#################################
### 2. Set working directory  ###
#################################
main_wd <- getwd()
setwd(main_wd)
input_dir <- "Analysis/Inputs/2_Processed_data"
output_dir <- "Analysis/Outputs"
routes_dir <- file.path(input_dir,"Routes_GO")
#################################
### 2. Set working directory  ###
#################################

main_wd <- getwd()
setwd(main_wd)
input_dir <- "Analysis/Inputs/2_Processed_data"
output_dir <- "Analysis/Outputs"
routes_dir <- file.path(input_dir,"Routes_GO")

###################################
### 3. Load GO genes per route  ###
###################################

list_files <- list.files(routes_dir,pattern = ".*.txt")
datalist = lapply(list_files, function(x) read.table(file.path(routes_dir,x),header=F))
names(datalist) <- gsub(".txt","", list_files)
### Fix the RV codes
datalist[[1]][2,] <- "Rv1099c"  
datalist[[3]][19,] <- "Rv3043c"

metadata<- read.table(file.path(input_dir,"txt/metadata_32_samples.txt"),check.names = F)
reads <- read.table(file.path(input_dir,"txt/reads_32_samples.txt"),check.names = F)
feature_data <- read.table(file.path(input_dir,"txt/feature_data_filtered.txt"),check.names = F)

reads <- reads[rownames(feature_data),]
length(which(rownames(feature_data)!=(rownames(reads))))
length(which(colnames(reads)!=rownames(metadata)))

### Sort by Culture and then by Iron keeping the pattern of Culture
metadata$Culture <- factor(metadata$Culture,levels = unique(metadata$Culture)) 

metadata <- metadata %>%
  arrange(Culture, Iron)

reads <- reads[,rownames(metadata)]
colnames(reads)

### Voom the reads
y <- DGEList(counts = reads)
y <- calcNormFactors(y)
design <- model.matrix(~short_setup,data=metadata)
v=voom(y,design,plot=TRUE)
exp=v$E

colnames(exp) <- gsub("_(\\d+)","",colnames(exp))
colnames(exp) <- gsub("C1_","G-D_EXP_",colnames(exp))
colnames(exp) <- gsub("C2","G-D_STAT",colnames(exp))
colnames(exp) <- gsub("C5","G_EXP",colnames(exp))
colnames(exp) <- gsub("C6","G_STAT",colnames(exp))
colnames(exp) <- gsub("C7","G-L_EXP",colnames(exp))
colnames(exp) <- gsub("C8","G-L_STAT",colnames(exp))
colnames(exp) <- gsub("C11","L_EXP",colnames(exp))
colnames(exp) <- gsub("C12","L_STAT",colnames(exp))

print(colnames(exp))

exp <- data.frame(exp)

# Create cluster by metabolic route 
clusters = lapply(datalist, function(x) {
  
  vec <- unlist(x,use.names = F)
  exp_filt <- exp[vec,]
  return(exp_filt)
  
})

glycolisis <- clusters[[1]]
methylcitrate <- clusters[[2]]
oxidative_phosph <- clusters[[3]]
pentose_phosphate <- clusters[[4]]
TCA <- clusters[[5]]

columns_name <- factor(colnames(glycolisis),levels =colnames(glycolisis))## Create a column names vector

# Plot
ht_by_route <- function(df,row_title="",column_title="",filename){
  
  ### Apply hclust function
  hclust_matrix <- as.matrix(df)
  dist_matrix <- dist(hclust_matrix)
  hc_clust <- hclust(dist_matrix,method="ward.D2")
  col_dend <- hclust(dist(t(hclust_matrix)))
  colnames(hclust_matrix) <- columns_name
  # column_title = ""
  # row_title = ""
  max(abs(df))      
  # color_breaks <- c(-(max(abs(df))-0.5) , 0, max(abs(df))-0.5 )
  color_breaks <- c(3, 9)
  
  # my_palette <- c( "yellow",
  #                  "blue",
  #                  colorRampPalette(rev(brewer.pal(8, "Spectral")))(n = length(color_breaks[2:4])-1),
  #                  "red")
  col_fun = colorRamp2(color_breaks, c("white","red"))
  # col_fun = colorRamp2(color_breaks, my_palette)
  
  ht_plt <- Heatmap(hclust_matrix,
                    na_col = "grey2",
                    col = col_fun,
                    # split = k,
                    name="Log2(OR)",
                    column_order = columns_name,
                    show_column_names = T,
                    column_names_gp = gpar(fontsize = 6),
                    row_names_gp = gpar(fontsize = 6),
                    column_title = column_title,
                    column_title_side = "bottom",
                    # row_dend_reorder=T,
                    border_gp = gpar(col = "black", lty = 2),
                    # heatmap_height = unit(6, "cm"),
                    # heatmap_width = unit(8, "cm"),
                    width = unit(0.8, "snpc"), 
                    height = unit(0.2, "snpc"),
                    # show_column_dend = T,
                    # column_dend_side = "top",
                    # cluster_rows = color_branches(hc_clust),
                    # cluster_columns = color_branches(col_dend),
                    row_title = row_title,
                    row_title_gp = gpar(fontize = 2),
                    row_dend_side = "right",
                    row_names_side = "left",
                    # row_dend_width = unit(2, "cm"),
                    show_row_names = T ,
                    show_row_dend = T,
                    # cell_fun = function(j, i, x, y, width, height, fill) {
                    #   grid.text(round(hclust_matrix[i, j],digits = 2), x, y, gp = gpar(fontsize = 3))},
                    heatmap_legend_param = list(title = "LogFC",
                                                title_position = "leftcenter-rot",
                                                labels_gp = gpar(font = 3), 
                                                title_gp = gpar( fontsize = 8),
                                                at= c(0,3,10)))
  
  
  draw(ht_plt)
  
  size <- calc_ht_size(ht_plt)
  
  dir.create(file.path(output_dir,"4_Figures_paper"),recursive = T,showWarnings=F)
  
  pdf(file =filename ,width=size[1]+0.5,height=size[2]+1)
  draw(ht_plt)
  dev.off()
  
  return(draw(ht_plt))
}

heat_dir <- "Heatmap_expression_levels"
filename <- file.path(output_dir,"4_Figures_paper",paste0(heat_dir,"_glycolisis.pdf"))
glycolisis_plt <- ht_by_route (glycolisis,row_title="",column_title="Growth arrest effects",filename)

# filename <- file.path(output_dir,"4_Figures_paper","15_Figure_4C_Heatmap_methylcitrate.pdf")
# methylcitrate_plt <- ht_by_route (methylcitrate,row_title="",column_title="Interactions",filename)

filename <- file.path(output_dir,"4_Figures_paper",paste0(heat_dir,"_oxidative_phosph.pdf"))
oxidative_phosph_plt <- ht_by_route (oxidative_phosph,row_title="",column_title="Interactions",filename)

# filename <- file.path(output_dir,"4_Figures_paper","15_Figure_4C_Heatmap_pentose_phosphate.pdf")
# pentose_phosphate_plt <- ht_by_route (pentose_phosphate,row_title="",column_title="Interactions",filename)

filename <- file.path(output_dir,"4_Figures_paper",paste0(heat_dir,"_TCA.pdf"))
TCA_plt <- ht_by_route (TCA,row_title="",column_title="Interactions",filename)

