#############################
### 0. Load dependencies. ###
#############################
{
  library(tidyverse)
  library(readxl)
  library(writexl)
}

############################
### 1. Declare functions ### 
############################
dcols=function(x){data.frame(colnames(x))}

###################################################
### 2. Set working directory and create folders ###
###################################################

main_wd <- getwd()
setwd(main_wd)
input_dir <- "Analysis/Inputs/2_Processed_data"
output_dir <- "Analysis/Outputs"

###########################
####### 3. Load data ######
###########################

### Matrix of contrasts
contrast <- read.table(file.path(input_dir,"txt/contrast_matrix.txt"))

### feature data
feat_path = file.path(input_dir,"txt/feature_data_filtered.txt")
feature_data = read.table(feat_path)

### metadata
meta_path <- file.path(input_dir,"txt/metadata_32_samples.txt")
metadata = read.table(meta_path)

### Contrast with the nomenclature used
df_name_contrast <- read.table(file.path(input_dir,"txt/contrasts_nomenclature.txt"))

### load the list
MyData <- readRDS(file = file.path(input_dir,"RDS/Contrasts_stat_0.05.RDS")) 
# mydf <- readRDS("Contrasts_stat.RDS")

### load the logFc,BH and lfcSE dfs
LogFC_df <- read.table(file.path(input_dir,"txt/LogFC_0.05.txt"))
BH_df <- read.table(file.path(input_dir,"txt/BH_0.05.txt"))
lfcSE_df <- read.table(file.path(input_dir,"txt/lfcSE_0.05.txt"))

### Extract the 8 Iron effects @  EXP and STAT
# Iron_effect_in_EXP_G_D.log2FoldChange_shrunken
# Iron_effect_in_EXP_G.log2FoldChange_shrunken
# Iron_effect_in_EXP_G_L.log2FoldChange_shrunken
# Iron_effect_in_EXP_L.log2FoldChange_shrunken
# Iron_effect_in_STAT_G_D.log2FoldChange_shrunken
# Iron_effect_in_STAT_G.log2FoldChange_shrunken
# Iron_effect_in_STAT_G_L.log2FoldChange_shrunken
# Iron_effect_in_STAT_L.log2FoldChange_shrunken

vec_iron <- c(5:8)
LogFC_df_iron <- LogFC_df[,vec_iron]
BH_df_iron <- BH_df[,vec_iron]
lfcSE_df_iron <- lfcSE_df[,vec_iron]

### Change the columns names
colnames(LogFC_df_iron) <- paste0(substring(colnames(LogFC_df_iron),1,nchar(colnames(LogFC_df_iron))-24),".LogFC")
# colnames(BH_df_iron) <- paste0(colnames(BH_df_iron),".BH")
colnames(lfcSE_df_iron) <- paste0(substring(colnames(lfcSE_df_iron),1,nchar(colnames(lfcSE_df_iron))-9))

length(which(rownames(LogFC_df_iron)!=rownames(BH_df_iron)))
length(which(rownames(LogFC_df_iron)!=rownames(lfcSE_df_iron)))

###Create Iron effects
iron_effects <- data.frame(LogFC_df_iron,BH_df_iron)

rv_non_coding <- iron_effects[grep("RVnc", rownames(iron_effects)), ]

{
  rv_non_coding$Pattern_G_D_stat   = "Background"
  rv_non_coding$Pattern_G_stat     = "Background"
  rv_non_coding$Pattern_G_L_stat   = "Background"
  rv_non_coding$Pattern_L_stat     = "Background"
  threshold <- 0.05
}

col_logFC <- c(1:4)
for (i in col_logFC){
  
  #Upregulated per condition
  rv_non_coding[which((rv_non_coding[(i+4)] < threshold) & 
                       (rv_non_coding[i] > 0)),(i+8)] <- "UP"
  #Downregulated per condition
  rv_non_coding[which(rv_non_coding[(i+4)] < threshold & 
                       rv_non_coding[i] < 0),(i+8)] <- "DOWN"
  
}

DE_genes_vec <- c()

for (i in c(9:12)){
  
  print(paste("There are", length(which(rv_non_coding[,i]=="UP" | 
                                          rv_non_coding[,i]=="DOWN")), "genes in",colnames(rv_non_coding)[i]))
  
  DE_genes_vec[i-8] <- length(which(rv_non_coding[,i]=="UP" | rv_non_coding[,i]=="DOWN"))
}

names_samples <- gsub('Pattern_','',colnames(rv_non_coding[9:12]))
# col <- c("lightcoral","red2","lightgoldenrod","yellow3","skyblue","blue4",
#          "lawngreen","darkgreen")
col <- c( "#852121", "#8f4a17","#2c682c","#355f91")


DE_genes_df <- data.frame(DE=DE_genes_vec,Contrast=names_samples,color=col)

### Save the list of DE genes in Lipids
write.table(rv_non_coding,file.path(input_dir,"txt/DE_genes_at_stat_RVnc_th<0.05.txt"))

rv_non_coding$RV <- rownames(rv_non_coding)
rv_non_coding <- rv_non_coding %>% relocate(RV, .before = Iron_effect_in_STAT_G_D.LogFC )
write_xlsx(rv_non_coding,file.path(input_dir,"xlsx/DE_genes_at_stat_RVnc_th<0.05.xlsx"))

df <- DE_genes_df
df$Contrast <- factor(df$Contrast,levels = df$Contrast)
breaks_values <- c(pretty(df$DE))

bar_plt <- ggplot(df, aes(y=DE, x=Contrast, fill= Contrast)) + 
  geom_bar(stat="identity",width = 0.8)+
  scale_fill_manual(values = col)+
  scale_y_continuous(expand = c(0,0), limits = c(0,10))+
  # scale_color_manual(values = col)+
  # geom_col_pattern(data=df[5:8,],pattern="stripe",
  #                  pattern_fill = color_vec[5:8],
  #                  pattern_angle = 45,
  #                  pattern_spacing = 0.05)+
  # scale_pattern_manual(values = color_vec[5:8])+
  # coord_flip()+
  theme_classic()+
  geom_hline(yintercept=0)+
  # theme(aspect.ratio = .9)+
  labs(x="Contrasts", y = "Number of Genes")+
  theme(
    axis.text.x   = element_text(size=12),
    axis.text.y   = element_text(size=12),
    axis.title.x  = element_text(size=12),
    axis.title.y  = element_text(size=12),
    # axis.ticks.y = element_blank(),
    #panel.background = element_blank(),
    #panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    #axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill=NA, linewidth=1,linetype="solid"),
    legend.title=element_blank(),
    legend.position="top",
    #legend.text=element_text(size=14),
    legend.key.size = unit(1, 'lines'))
bar_plt

bar_pl_1 <- bar_plt+
  geom_shadowtext(
    data = df,
    aes(x=Contrast, y=DE, label = DE),
    hjust = 0.5,vjust = -0.5,
    position = "stack",
    # nudge_x = -0.3,
    colour = "black",
    bg.colour = "white",
    bg.r = 0.2,
    # family = "Econ Sans Cnd",
    size = 5)

bar_pl_1

pdf(file=file.path(output_dir,"4_Figures_paper","barplot_DE_noncod_genes_stat_effects.pdf"),width=7,height=7) 
print(bar_pl_1)
dev.off()
