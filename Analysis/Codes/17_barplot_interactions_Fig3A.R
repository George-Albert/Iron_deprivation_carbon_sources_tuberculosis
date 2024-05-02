#############################
### 0. Load dependencies. ###
#############################

{ 
  library(tidyverse)
  library(readxl)
  library(writexl)
  library(eulerr)
  library(shadowtext)
}

############################
### 1. Declare functions ### 
############################
dcols=function(x){data.frame(colnames(x))}
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

### Load DE genes.txt of interest
DE_genes_count <- read.table(file=file.path(input_dir,"txt/DE_genes_of_interest_th<0.05.txt"))

### load the logFc,BH and lfcSE dfs
LogFC_df <- read.table(file.path(input_dir,"txt/LogFC_0.05.txt"))
BH_df <- read.table(file.path(input_dir,"txt/BH_0.05.txt"))
lfcSE_df <- read.table(file.path(input_dir,"txt/lfcSE_0.05.txt"))

col <- c("#852121","#8f4a17","#2c682c","#355f91")

# rownames(DE_genes_count[9:12,])<- factor(rownames(DE_genes_count[9:12,]),levels=rownames(DE_genes_count[9:12,]))

bar_plt <- ggplot(DE_genes_count[9:12,], aes(x = factor(rownames(DE_genes_count[9:12,]),levels=rownames(DE_genes_count[9:12,])), y = DE, fill=factor(rownames(DE_genes_count[9:12,]),levels=rownames(DE_genes_count[9:12,])))) + 
  geom_bar(stat="identity",position = "dodge",width =2/3)+
  scale_fill_manual(values = col)+
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
    axis.text.x   = element_text(size=12,angle=45,vjust = 0.6),
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
    data = DE_genes_count[9:12,],
    aes(x = factor(rownames(DE_genes_count[9:12,]),levels=rownames(DE_genes_count[9:12,])), y = DE,
        label=DE),
    hjust = 0.5,vjust=-0.5,
    position = "stack",
    # nudge_x = -0.3,
    colour = "black",
    bg.colour = "white",
    bg.r = 0.2,
    # family = "Econ Sans Cnd",
    size = 5)

bar_pl_1

pdf(file = file.path(output_dir,"4_Figures_paper","10_Figure_3A_Barplot_DE_genes_interaction.pdf"),width=6,height=7.5)
print(bar_pl_1)
dev.off()
