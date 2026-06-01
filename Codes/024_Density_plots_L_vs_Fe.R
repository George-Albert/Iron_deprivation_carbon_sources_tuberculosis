######################
### 0.Dependencies ###
######################
{
  library(tidyverse)
  library(shadowtext)
  library(cowplot)
}

###################################################
### 2. Set working directory and create folders ###
###################################################
main_wd <- getwd()
setwd(main_wd)
input_dir  <- "Inputs/002_Processed_data"
output_dir <- "Outputs"

###########################
####### 3. Load data ######
###########################
metadata         <- read.table(file.path(input_dir,"txt/metadata_32_samples.txt"),check.names = F)
feature_data     <- read.table(file.path(input_dir,"txt/feature_data_genes_type_added.txt"),check.names = F)
df_name_contrast <- read.table(file.path(input_dir,"txt/contrasts_nomenclature.txt"))

### Lets filter the metadata to retrieve only C7-C8 and C11-C12
metadata_filtered <- metadata %>%
  filter((Culture %in% c("C7","C8")) | (Culture %in% c("C11","C12")))

# Get the batch info
batch <- metadata_filtered$Full_path
batch <- gsub("Mtb_Latency/","",batch)
batch <- gsub("/.*","",batch)

unique(batch)

# [1] "HN00161701_C7_C8_Diciembre_2021"   "HN00170106_C14_C15_C8r_Mayo_2022"
# [3] "HN00189620_C11_C15_Abril_2023"     "HN00166209_C11_C12_C2r_Marzo_2022"

# Get the contrasts names of interest and its file names
file_names <- paste0(df_name_contrast$title,".txt")
file_names <- file_names[c(11:12,15:16,33:34,37,41)]
current_folder <- file.path(output_dir,"Data/txt/th_0.05_th_size_0")
### Get full names including folder path
list.of.files = list.files(current_folder, full.names = TRUE)
### Keep only the basename (file names) matching dataframe column
clean_list <- list.of.files[basename(list.of.files) %in% file_names]
data_name <- basename(clean_list)
data_name <- substr(data_name,1,nchar(data_name)-4)

### Read the data
myData <- lapply(clean_list, read.table)
myData <- setNames(object = myData, data_name)

CS_effect_wo_Fe <- myData[[1]][,c("log2FoldChange_shrunken","BH")]
colnames(CS_effect_wo_Fe) <- c("LogFC","BH")

CS_effect_w_Fe <- myData[[2]][,c("log2FoldChange_shrunken","BH")]
colnames(CS_effect_w_Fe) <- c("LogFC","BH")

feature_data_filtered <- feature_data[rownames(CS_effect_wo_Fe),]
length(which(rownames(CS_effect_wo_Fe)!=rownames(feature_data_filtered)))
length(which(rownames(CS_effect_w_Fe)!=rownames(feature_data_filtered)))

# CS_effect_wo_Fe$Type <- feature_data_filtered$Type
# CS_effect_w_Fe$Type   <- feature_data_filtered$Type

th <- 0.05
density_plotter=function(stat,feat=feature_data,threshold=th,culture){
  stat$Type=feat$Type
  stat=stat[which(stat$BH<threshold),]

  pc=stat$LogFC[which(stat$Type=="CDS")]
  npc=stat$LogFC[which(stat$Type=="Stable_RNA")]
  test=wilcox.test(npc,pc,alt="greater")
  print(paste0("Wilcox rank test p=",test$p.value," (alt: npc>pc)"))


  pl=ggplot(stat)+geom_density(aes(x=LogFC,fill=Type),alpha=0.8)+
    xlim(-10,10)+
    ylim(0,0.4)+
    xlab(paste0("logFC Carbon source response at ",culture))+
    theme(legend.position="none")+
    geom_vline(xintercept=0)+
    theme_classic()+
    theme(
      axis.text.y   = element_text(size=14),
      axis.text.x   = element_text(size=14),
      axis.title.y  = element_text(size=14),
      axis.title.x  = element_text(size=14),
      #panel.background = element_blank(),
      #panel.grid.major = element_blank(),
      #panel.grid.minor = element_blank(),
      #axis.line = element_line(colour = "black"),
      panel.border = element_rect(colour = "black", fill=NA, linewidth = 1,linetype="solid"),
      #legend.title=element_blank(),
      legend.position="top",
      #legend.text=element_text(size=14),
      legend.key.size = unit(1, 'lines'))

  return(pl)
}


effect_name <- "Carbon_source_effects_in_absence_and_presence_of_Fe_GLYCEROL_LCFA-LCFA"
print(paste("Processing:", effect_name))

# Create dataframes for each condition
stat_noFe <- CS_effect_wo_Fe
stat_Fe   <- CS_effect_w_Fe

pl_noFe <- density_plotter(stat_noFe,feat = feature_data_filtered, culture = "Fe/-")
pl_Fe   <- density_plotter(stat_Fe,  feat = feature_data_filtered, culture = "Fe/+")

combined_pl <- plot_grid(pl_noFe, pl_Fe, nrow = 2, align = "hv")
suplementary_dir <- file.path(output_dir, "001_Figures_paper","Suplementary_Figures")
dir.create(suplementary_dir, showWarnings = FALSE)
pdf(file = file.path(suplementary_dir, paste0("Density_", effect_name, "_Fe_effects.pdf")),
    width = 6, height = 6)
print(combined_pl)
dev.off()



#########################################################
##  Same analysis but now with glycerol and LCFA only  ##
#########################################################

### Lets filter the metadata to retrieve only C7-C8 and C11-C12
metadata_filtered <- metadata %>%
  filter((Culture %in% c("C5","C6")) | (Culture %in% c("C11","C12")))

# Get the batch info
batch <- metadata_filtered$Full_path
batch <- gsub("Mtb_Latency/","",batch)
batch <- gsub("/.*","",batch)

unique(batch)

#  "HN00150037_C5_C6       Junio_2021"       "HN00189620_C11_C15  Abril_2023"
#  "HN00166209_C11_C12_C2r Marzo_2022"

# Get the contrasts names of interest and its file names
file_names <- paste0(df_name_contrast$title,".txt")
file_names <- file_names[50:51]
current_folder <- file.path(output_dir,"Data/txt/th_0.05_th_size_0")
### Get full names including folder path
list.of.files = list.files(current_folder, full.names = TRUE)
### Keep only the basename (file names) matching dataframe column
clean_list <- list.of.files[basename(list.of.files) %in% file_names]
data_name <- basename(clean_list)
data_name <- substr(data_name,1,nchar(data_name)-4)

### Read the data
myData <- lapply(clean_list, read.table)
myData <- setNames(object = myData, data_name)

CS_effect_wo_Fe <- myData[[1]][,c("log2FoldChange_shrunken","BH")]
colnames(CS_effect_wo_Fe) <- c("LogFC","BH")

CS_effect_w_Fe <- myData[[2]][,c("log2FoldChange_shrunken","BH")]
colnames(CS_effect_w_Fe) <- c("LogFC","BH")

length(which(rownames(CS_effect_wo_Fe)!=rownames(feature_data_filtered)))
length(which(rownames(CS_effect_w_Fe)!=rownames(feature_data_filtered)))

# CS_effect_wo_Fe$Type <- feature_data_filtered$Type
# CS_effect_w_Fe$Type   <- feature_data_filtered$Type

th <- 0.05
density_plotter=function(stat,feat=feature_data,threshold=th,culture){
  stat$Type=feat$Type
  stat=stat[which(stat$BH<threshold),]

  pc=stat$LogFC[which(stat$Type=="CDS")]
  npc=stat$LogFC[which(stat$Type=="Stable_RNA")]
  test=wilcox.test(npc,pc,alt="greater")
  print(paste0("Wilcox rank test p=",test$p.value," (alt: npc>pc)"))


  pl=ggplot(stat)+geom_density(aes(x=LogFC,fill=Type),alpha=0.8)+
    xlim(-10,10)+
    ylim(0,0.4)+
    xlab(paste0("logFC Carbon source response at ",culture))+
    theme(legend.position="none")+
    geom_vline(xintercept=0)+
    theme_classic()+
    theme(
      axis.text.y   = element_text(size=14),
      axis.text.x   = element_text(size=14),
      axis.title.y  = element_text(size=14),
      axis.title.x  = element_text(size=14),
      #panel.background = element_blank(),
      #panel.grid.major = element_blank(),
      #panel.grid.minor = element_blank(),
      #axis.line = element_line(colour = "black"),
      panel.border = element_rect(colour = "black", fill=NA, linewidth = 1,linetype="solid"),
      #legend.title=element_blank(),
      legend.position="top",
      #legend.text=element_text(size=14),
      legend.key.size = unit(1, 'lines'))

  return(pl)
}


effect_name <- "Carbon_source_effects_in_absence_and_presence_of_Fe_GLYCEROL-LCFA"
print(paste("Processing:", effect_name))

# Create dataframes for each condition
stat_noFe <- CS_effect_wo_Fe
stat_Fe   <- CS_effect_w_Fe

pl_noFe <- density_plotter(stat_noFe,feat = feature_data_filtered, culture = "Fe/-")
pl_Fe   <- density_plotter(stat_Fe,  feat = feature_data_filtered, culture = "Fe/+")

combined_pl <- plot_grid(pl_noFe, pl_Fe, nrow = 2, align = "hv")
suplementary_dir <- file.path(output_dir, "001_Figures_paper","Suplementary_Figures")
dir.create(suplementary_dir, showWarnings = FALSE)
pdf(file = file.path(suplementary_dir, paste0("Density_", effect_name, "_Fe_effects.pdf")),
    width = 6, height = 6)
print(combined_pl)
dev.off()
