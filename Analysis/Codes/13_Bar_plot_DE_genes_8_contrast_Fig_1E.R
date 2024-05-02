######################
### 0.Dependencies ###
######################
{
  library(tidyverse)
  library(shadowtext)
}

###########################
### 1.Declare Functions ###
###########################

dcols=function(x){data.frame(colnames(x))}

#################################
### 2. Set working directory  ###
#################################


main_wd <- getwd()
setwd(main_wd)
input_dir <- "Analysis/Inputs/2_Processed_data"
output_dir <- "Analysis/Outputs"

#####################
### 3. Load data  ###
#####################

DE_genes_per_cond <- read.table(file = file.path(input_dir,"txt/DE_genes_per_cond_th<0.05.txt"))

df <- DE_genes_per_cond[c(9:16),]
df$Conditions <- rownames(df)
df$Conditions <- factor(df$Conditions,levels = df$Conditions)

color_vec <- rep(c("#852121", "#8f4a17", "#2c682c", "#355f91"),2)
fill_vec <- c("#852121", "#8f4a17", "#2c682c", "#355f91",rep("white",4))
breaks_values <- c(pretty(df$DE))

bar_plt <- ggplot(df, aes(y=DE, x=Conditions, color = Conditions,fill=Conditions)) + 
  geom_bar(stat="identity",width = 0.8,linewidth=2)+
  # coord_flip()+
  # scale_y_continuous(breaks = breaks_values,labels = abs(breaks_values),limits=c(0,3500))+
  scale_fill_manual(values = fill_vec)+
  scale_color_manual(values = color_vec)+
  theme_classic()+
  geom_hline(yintercept=0)+
  # theme(aspect.ratio = .9)+
  labs(x="Contrasts", y = "Number of Genes")+
  theme(
    axis.text.x   = element_text(size=12,angle = 45,vjust=1.1,hjust = 1.1),
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

bar_pl_1 <- bar_plt+
  geom_shadowtext(
    data = df,
    aes(x=Conditions, y=DE, label = DE),
    hjust = 0.5,vjust=-0.5,
    position = "stack",
    # nudge_x = -0.3,
    colour = "black",
    bg.colour = "white",
    bg.r = 0.2,
    # family = "Econ Sans Cnd",
    size = 5)

bar_pl_1

pdf(file = file.path(output_dir,"4_Figures_paper","3_Figure_S1A_Bar_plot_DE_genes.pdf"),width=9,height=7)
print(bar_pl_1)
dev.off()
