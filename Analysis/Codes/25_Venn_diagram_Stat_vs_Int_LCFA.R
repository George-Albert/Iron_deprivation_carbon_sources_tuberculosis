#############################
### 0. Load dependencies. ###
#############################

library(eulerr)
library(tidyverse)

############################
### 1. Declare functions ### 
############################

dcols=function(x){data.frame(colnames(x))}

#################################
### 2. Set working directory  ###
#################################

getwd()
main_wd    <-getwd()
setwd(main_wd)
input_dir  <- "Analysis/Inputs/2_Processed_data"
output_dir <- "Analysis/Outputs"

###########################
####### 3. Load data ######
###########################

### Load genes name of the LCFA interaction
L_int_genes <- read.table(file.path(input_dir,"txt/DE_genes_in_Lipids_th<0.05.txt"))
L_int_genes[order(L_int_genes$genes),]
### Load genes name of the LCFA iron effects at STAT
L_stat_genes <- read.table(file.path(input_dir,"txt/DE_genes_iron_effects_in_lipids_th<0.05.txt"))
L_stat_genes[order(L_stat_genes$x),]
Genes_list <- list("Iron effects at STAT (LCFA)"=L_stat_genes$x,
                   "Interaction (LCFA)" = L_int_genes$genes)




# plot(venn(Genes_list))
venn2 <- plot(euler(Genes_list),quantities=list(type = c("counts", "percent"),col="black", font=4, round=2, cex=1.0),
              fills = list(fill = c("red", "grey"), alpha = 0.6),
              labels = list(col = "black", fontsize = 10),
              col="black",
              lty = 2:1,
              lwd=3,
              legend = list(labels = c("Iron effects at STAT (LCFA)", "Interaction (LCFA)")))

venn2

pdf(file.path(output_dir,"4_Figures_paper","Venn_Stat_vs_Int_LCFA.pdf"),width=12,height=7)
print(venn2)
dev.off()

















