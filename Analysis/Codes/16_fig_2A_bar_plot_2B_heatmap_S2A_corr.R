#############################
### 0. Load dependencies. ###
#############################
{
  library(readxl)
  library(writexl)
  library(tidyverse)
  library(cowplot)
  library(ComplexHeatmap)
  library(shadowtext)
  library(circlize)
  library(corrplot)
  library(dendextend)
  library(ggthemes)
  library(plotly)
  library(ggdendro)
}

############################
### 1. Declare functions ### 
############################
dcols=function(x){data.frame(colnames(x))}
### Compute the size of the heatmap
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
# iron_response_to_stat <- read.table(file.path(input_dir,"txt/iron_response_to_stat.txt"))

### Extract the 8 Iron effects @  EXP and STAT
# Iron_effect_in_EXP_G_D.log2FoldChange_shrunken
# Iron_effect_in_EXP_G.log2FoldChange_shrunken
# Iron_effect_in_EXP_G_L.log2FoldChange_shrunken
# Iron_effect_in_EXP_L.log2FoldChange_shrunken
# Iron_effect_in_STAT_G_D.log2FoldChange_shrunken
# Iron_effect_in_STAT_G.log2FoldChange_shrunken
# Iron_effect_in_STAT_G_L.log2FoldChange_shrunken
# Iron_effect_in_STAT_L.log2FoldChange_shrunken

vec_iron <- c(1:8)
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
write.table(iron_effects,file.path(input_dir,"txt/iron_effects.txt"))

{
  iron_effects$Pattern_G_D_exp    = "Background"
  iron_effects$Pattern_G_exp      = "Background"
  iron_effects$Pattern_G_L_exp    = "Background"
  iron_effects$Pattern_L_exp      = "Background"
  iron_effects$Pattern_G_D_stat   = "Background"
  iron_effects$Pattern_G_stat     = "Background"
  iron_effects$Pattern_G_L_stat   = "Background"
  iron_effects$Pattern_L_stat     = "Background"
  threshold <- 0.05
}

### Separate by direction (UP or DOWN)
col_logFC <- c(1:8)
for (i in col_logFC){
  
  #Upregulated per condition
  iron_effects[which((iron_effects[(i+8)] < threshold) & 
                                (iron_effects[i] > 0)),(i+16)] <- "UP"
  #Downregulated per condition
  iron_effects[which(iron_effects[(i+8)] < threshold & 
                       iron_effects[i] < 0),(i+16)] <- "DOWN"
  
}

### DE genes count
DE_genes_vec <- c()

# length(which(iron_effects[,"Pattern_G_L_exp"]=="UP"))
# length(which(iron_effects[,"Pattern_G_L_exp"]=="DOWN"))
# 
# length(which(iron_effects[,"Pattern_G_exp"]=="UP"))
# length(which(iron_effects[,"Pattern_G_L_exp"]=="DOWN"))
# 
# length(which(iron_effects[,"Pattern_L_exp"]=="UP"))
# length(which(iron_effects[,"Pattern_L_exp"]=="DOWN"))
# 
# # iron_effects["Rv2383c",c(9:12)]
# iron_effects[which(iron_effects[,"Pattern_L_exp"]=="UP"),]


for (i in c(17:24)){
  
  print(paste("There are", length(which(iron_effects[,i]=="UP" | 
                                          iron_effects[,i]=="DOWN")), "genes in",colnames(iron_effects)[i]))
  
  DE_genes_vec[i-16] <- length(which(iron_effects[,i]=="UP" | iron_effects[,i]=="DOWN"))
  }


#################################################################
##           Figure_2A_Barplot_DE_genes_iron_effects           ##
#################################################################

names_samples <- gsub('Pattern_','',colnames(iron_effects[17:24]))
# col <- c("lightcoral","red2","lightgoldenrod","yellow3","skyblue","blue4",
#          "lawngreen","darkgreen")
col <- c( "#f6a1a1", "#ffc499","#a9dea9","#8cd9ff",
          "#852121", "#8f4a17","#2c682c","#355f91")


DE_genes_df <- data.frame(DE=DE_genes_vec,Contrast=names_samples,color=col)
### Extract the DE genes in LCFA
iron_effects_in_lipids_name <- rownames(iron_effects[which(iron_effects[,"Pattern_L_stat"]=="UP" | iron_effects[,"Pattern_L_stat"]=="DOWN"),])
iron_effects_in_lipids <- iron_effects[which(iron_effects[,"Pattern_L_stat"]=="UP" | iron_effects[,"Pattern_L_stat"]=="DOWN"),c(8,16,24)]
iron_effects_in_lipids$RV <- rownames(iron_effects_in_lipids)
### Save the list of DE genes in Lipids
write.table(iron_effects_in_lipids_name,file.path(input_dir,"txt/DE_genes_iron_effects_in_lipids_th<0.05.txt"))
write_xlsx(data.frame(iron_effects_in_lipids_name),file.path(input_dir,"xlsx/DE_genes_iron_effects_in_lipids_th<0.05.xlsx"))

write.table(iron_effects_in_lipids,file.path(input_dir,"txt/DE_lfc_BH_iron_effects_in_lipids_th<0.05.txt"))
write_xlsx(iron_effects_in_lipids,file.path(input_dir,"xlsx/DE_lfc_BH_iron_effects_in_lipids_th<0.05.xlsx"))

write.table(DE_genes_df,file = file.path(input_dir, paste0("txt/DE_genes_iron_effects_",threshold,".txt")))

df <- DE_genes_df
df$Contrast <- factor(df$Contrast,levels = df$Contrast)
breaks_values <- c(pretty(df$DE))

bar_plt <- ggplot(df, aes(y=DE, x=Contrast, fill= Contrast)) + 
  geom_bar(stat="identity",width = 0.8)+
  scale_fill_manual(values = col)+
  scale_y_continuous(expand = c(0,0), limits = c(0,2500))+
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

pdf(file=file.path(output_dir,"4_Figures_paper","6_Figure_2A_Barplot_DE_genes_iron_effects.pdf"),width=7,height=7) 
print(bar_pl_1)
dev.off()

#################################################################
##                     Figure 2B Heatmap                       ##
#################################################################

dcols(iron_effects)

iron_effects_at_exp <- iron_effects[,c(1:4,9:12,17:20)]
iron_effects_at_exp_filt <- iron_effects_at_exp[which(iron_effects_at_exp$Pattern_G_D_exp   == "UP" | iron_effects_at_exp$Pattern_G_D_exp == "DOWN" |
                                                        iron_effects_at_exp$Pattern_G_exp   == "UP" | iron_effects_at_exp$Pattern_G_exp   == "DOWN" |
                                                        iron_effects_at_exp$Pattern_G_L_exp == "UP" | iron_effects_at_exp$Pattern_G_L_exp == "DOWN" |
                                                        iron_effects_at_exp$Pattern_L_exp   == "UP" | iron_effects_at_exp$Pattern_L_exp   == "DOWN" ),]
feature_data_filt <- feature_data[rownames(iron_effects_at_exp_filt),]
feature_data_filt["Rv2377c",]$symbol <- "mbtH"

### Create the vector of genes names and fill the empty spaces with the RV codes
symbol <- feature_data_filt$symbol
symbol[which(symbol=="")] <- NA
RV <- rownames(feature_data_filt) 
symbol_fill <- coalesce(symbol,RV)

rownames(iron_effects_at_exp_filt) <- symbol_fill
BH_exp <- iron_effects_at_exp_filt[,5:8]

### Construct the dendrogram
# Run clustering
hc_matrix <- as.matrix(iron_effects_at_exp_filt[,1:4])
hc_dendro <- as.dendrogram(hclust(d = dist(x = hc_matrix),method = "ward.D2"))

# Create dendro
dendro_plot <- ggdendrogram(data = hc_dendro, rotate = TRUE,labels = F)+
  scale_x_continuous(expand = c(0, 0.1)) +
  scale_y_continuous(expand = c(0, 0))+
  labs(x = "")

# Preview the plot
print(dendro_plot)

### Construct the heatmap
data <- iron_effects_at_exp_filt[,1:4]
columns_name <- c("G_D_EXP","G_EXP","G_LCFA_EXP","LCFA_EXP")
colnames(data) <- columns_name 
data$RV <- rownames(data)
dendro_order <- order.dendrogram(hc_dendro)
data$RV <- factor(x=data$RV,levels = data$RV[dendro_order],ordered = T)

#Long format logfc
data_long <- pivot_longer(data = data, cols = !RV, values_to = "LogFC" ,names_to = "Conditions")
data_long$Conditions <- factor(x=data_long$Conditions,levels = unique(data_long$Conditions))

data_BH <- BH_exp
columns_name <- c("G_D_EXP","G_EXP","G_LCFA_EXP","LCFA_EXP")
colnames(data_BH) <- columns_name 
data_BH$RV <- rownames(data_BH)
length(which(rownames(data)!=rownames(data_BH)))
data_BH$RV <- factor(x=data_BH$RV,levels = data_BH$RV[dendro_order],ordered = T)

#Long format BH
data_long_BH <- pivot_longer(data = data_BH, cols = !RV, values_to = "BH" ,names_to = "Conditions")
data_long_BH$Conditions <- factor(x=data_long_BH$Conditions,levels = unique(data_long_BH$Conditions))

length(which(data_long$RV != data_long_BH$RV))
data_long$BH <- data_long_BH$BH

base_size <- 9

ht_plt <- ggplot(data=data_long, aes(x=Conditions, y=RV, col= LogFC)) + 
  geom_tile(col="black", fill="grey99",linewidth=1.2)+
  geom_point(aes(size = -log10(BH)), shape=15)+
  # geom_text(data=data_long, aes(x=Conditions, y=RV, label = round(LogFC,2)), color = "black", size = 2)+
  coord_equal()+
  scale_size(range = c(1, 5),limits=c(0,6),breaks = c(0,2,4,6)) +
  scale_color_gradient2(mid="#f5d7ed",low="#0C6291",high="red3", limits=c(-1.6,1.6))+
  # scale_fill_distiller(palette = "RdBu",limit = c(-1.6,1.6)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0,0))+
  theme_grey(base_size = base_size)+
  labs(x=NULL, y=NULL) +
  theme_tufte(base_family="Helvetica")
ht_plt

dendro_plot1 <- dendro_plot + 
  theme(axis.text = element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        plot.margin = unit(c(t = 0, r = 5, b = 1.5, l = 0),unit = "cm"))

dendro_plot1

ht_plt1 <- ht_plt + theme(legend.position = "left", 
                          axis.ticks = element_blank(),
                          axis.text.x = element_text(size = base_size * 0.8, angle = 270,
                                                     hjust = 0, colour = "grey50"),
                          # text=element_text(family="Roboto"),
                          plot.margin = unit(c(t = 0, r = -0.5, b = 0, l = 0),unit = "cm"))

ht_plt1

grid_plt <- plot_grid(ht_plt1,NULL,dendro_plot1,ncol = 3, rel_widths = c(1,0.01,1))
grid_plt

dir.create(file.path(output_dir,"4_Figures_paper"))
filename=file.path(output_dir,"4_Figures_paper","7_Fig_2B_Heatmap_iron_effects_at_exp.pdf")

pdf(file =filename )
print(grid_plt)
dev.off()

filename=file.path(output_dir,"4_Figures_paper","7_Fig_2B_Heatmap_iron_effects_at_exp_wo_dend.pdf")

pdf(file =filename )
print(ht_plt1)
dev.off()

##################################################################
##           Corrrelation interactions vs iron effects          ##
##################################################################

lfc_int_vs_iron_effects <- LogFC_df[,c(5:8,17:20)]
BH_int_vs_iron_effects <- BH_df[,c(5:8,17:20)]
### Change the columns names
colnames(lfc_int_vs_iron_effects) <- paste0(substring(colnames(lfc_int_vs_iron_effects),1,nchar(colnames(lfc_int_vs_iron_effects))-24),".LogFC")

length(which(rownames(lfc_int_vs_iron_effects)!=rownames(BH_df_iron)))

int_vs_iron_effects <- data.frame(lfc_int_vs_iron_effects,BH_int_vs_iron_effects)

### Define colors
color <- c("#852121", "#8f4a17","#2c682c",  "#355f91") 

### Create an empty cor matrix
cor.matrix <- matrix(NA, nrow = 4, ncol = 1)
pvalues.matrix <- matrix(NA, nrow = 4, ncol = 1)

names_cor_matrix <- c()

### number of interaction to compare with
vector_to_loop <- 1:4
th_1 <- 0.05
th_2 <- 0.05

for (col in vector_to_loop) {
  
  print(paste(colnames(int_vs_iron_effects)[col],"vs",colnames(int_vs_iron_effects)[col+4]))

  df <- data.frame(LogFC1 = int_vs_iron_effects[,col],
                   LogFC2=int_vs_iron_effects[,col+4],
                   BH1=int_vs_iron_effects[,col+8],
                   BH2= int_vs_iron_effects[,col+12],
                   row.names = rownames(int_vs_iron_effects))


  df_filtered <-df[which(df$BH1 <th_1 & df$BH2 < th_2),] 
  
  ### Create dir
  df_to_corr_dir <- paste0(colnames(int_vs_iron_effects)[col],"_vs_",colnames(int_vs_iron_effects)[col+4],".txt")
  dir.create(file.path(input_dir,"txt/Correlogram_data_FigS2A"),showWarnings = F)
  write.table(df_filtered,file.path(input_dir,"txt/Correlogram_data_FigS2A",df_to_corr_dir))
  
  ####################
  ###Scatter plots ###
  ####################
  ### Define axis
  y <- df_filtered$LogFC2
  x <- df_filtered$LogFC1
  ### Define axis title
  x_title <- colnames(int_vs_iron_effects)[col]
  y_title <- colnames(int_vs_iron_effects)[col+4]
  
  ### Pearson Correlation
  test=cor.test(y,x,method = "pearson")
  print(test)
  ### Save estimate and pvalue
  names_cor_matrix[col] <- paste0(colnames(int_vs_iron_effects)[col],"_vs_",colnames(int_vs_iron_effects)[col+4],".txt")
  cor.matrix[col,] <- test$estimate
  pvalues.matrix[col,] <- -log10(test$p.value)
  
  ### plot
  pl_scatter=ggplot(df_filtered)+
    geom_point(aes(x=LogFC1,y=LogFC2),color=color[col],alpha=0.8)+
    # geom_abline(slope=0,intercept=0,color="black")+
    xlab(x_title)+
    ylab(y_title)+
    theme(legend.position=c(0.8,0.2))
  
  pl_scatter_1 <- pl_scatter+
    geom_smooth(aes(x=LogFC1,y=LogFC2),color="black",method=lm,formula = "y ~ x",show.legend = FALSE)+
    # geom_rug(aes(x=LogFC1,y=LogFC2),color=color)+
    geom_vline(xintercept = 0)+
    geom_hline(yintercept = 0)+
    ggtitle(paste0("R = ",round(test$estimate,digits=2)," Num of common genes = ",dim(df_filtered)[1]))+
    theme_minimal()+
    theme(plot.title = element_text(size = 16, face = "bold"))
  
  print(pl_scatter_1)
  scatter_dir <- "4_Figures_paper/9_Figure_S2A_Scatter_plot_correlation_stat_vs_interactions"
  dir.create(file.path(output_dir,scatter_dir),recursive = T,showWarnings = F)
  
  pdf(file.path(output_dir,scatter_dir,paste0(x_title,"_vs_",y_title,".pdf")),width=6,height=5)
  print(pl_scatter_1)
  dev.off()
}




