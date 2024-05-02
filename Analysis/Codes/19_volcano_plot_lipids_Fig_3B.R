#############################
### 0. Load dependencies. ###
#############################

{
  library(tidyverse)
  library(ggplot2)
  library(EnhancedVolcano)
  library(writexl)
}

############################
### 1. Declare functions ###
############################
dcols=function(x){data.frame(colnames(x))}
volcan_plot <- function(data,x,y,xintercept,th,ymax,title){
  
  datos <- data.frame(data)
  y <- -(log10(y))
  y[is.infinite(y)] <- NA
  
  sizes <- c("UP" = 2, "DOWN" = 2, "Not Sig." = 0.5) 
  alphas <- c("UP" = 1, "DOWN" = 1, "Not Sig." = 0.5)
  #FFA373
  ##50486D
  mycolors <- c("UP" ="green","Not Sig." ="grey","DOWN" ="magenta")
  
  ggplot(datos,aes(x=x,y=y,
                   size=Association,
                   alpha=Association,
                   fill=Association))+
    geom_point(shape = 21, colour = "black", size = 3, alpha = 0.8)+
    ylab("-log10(FDR)")+xlab("Log2FC")+theme_minimal()+
    scale_fill_manual(values=mycolors)+
    scale_size_manual(values = sizes,guide="none")+
    scale_alpha_manual(values = alphas, guide="none")+
    geom_vline(xintercept=c(-xintercept, xintercept),col="black",linetype = "dashed")+
    geom_hline(yintercept=-log10(th), col="black",linetype = "dashed")+
    scale_x_continuous(limits = c(-(max(x)), max(x)))+
    scale_y_continuous(limits = c(-1, ymax),expand = expansion(0))+
    ggtitle(title)+
    theme(plot.title = element_text(hjust = 0.5),
          text = element_text(size = 16),
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 14))
}

#################################
### 2. Set working directory  ###
#################################

main_wd <- getwd()
setwd(main_wd)
input_dir <- "Analysis/Inputs/2_Processed_data"
output_dir <- "Analysis/Outputs"

####################
### 2. Load data ###
####################

### feature data
feature_data <- read.table(file.path(input_dir,"txt/feature_data_filtered.txt"))
feature_data["Rv2377c",]$symbol <- "mbtH"

### metadata
meta_path <- file.path(input_dir,"txt/metadata_32_samples.txt")
metadata = read.table(meta_path)

### Contrast with the nomenclature used
df_name_contrast <- read.table(file.path(input_dir,"txt/contrasts_nomenclature.txt"))

### load the list of statistics from DeSEQ2
MyData <- readRDS(file = file.path(input_dir,"RDS/Contrasts_stat_0.05.RDS")) 

### load the logFc,BH and lfcSE dfs
LogFC_df <- read.table(file.path(input_dir,"txt/LogFC_0.05.txt"))
BH_df <- read.table(file.path(input_dir,"txt/BH_0.05.txt"))
lfcSE_df <- read.table(file.path(input_dir,"txt/lfcSE_0.05.txt"))

### vec_int is the columns of the interactions we are going to use
vec_int <- 20

int_LogFC <- LogFC_df[,vec_int]
int_BH    <- BH_df[,vec_int]
int_lfcSE <- lfcSE_df[,vec_int]

lipid_df <- data.frame(LogFC=int_LogFC,BH=int_BH,row.names = rownames(LogFC_df))

### Set the threshold
th=0.05

# lipid_df_filt <- lipid_df[which(lipid_df$BH < th),]
lipid_df_filt <- lipid_df
lipid_df_filt$Association <- "Background"
lipid_df_filt[which(lipid_df_filt$LogFC>0),]$Association <- "UP"
lipid_df_filt[which(lipid_df_filt$LogFC<0),]$Association <- "DOWN"
lipid_df_filt[which(lipid_df$BH > th),"Association"] <- "Not Sig."

#######################
### 3. Volcano plot ###
#######################

name <- "LCFA"
genes <- rownames(lipid_df_filt[which(lipid_df_filt$Association=="UP" | lipid_df_filt$Association=="DOWN"),])
feature_data_filt <- feature_data[genes,]

df_genes <- lipid_df_filt[which(lipid_df_filt$Association=="UP" | lipid_df_filt$Association=="DOWN"),]
length(which(rownames(df_genes)!=rownames(feature_data_filt)))

### Create the vector of genes names and fill the empty spaces with the RV codes
symbol <- feature_data_filt$symbol
symbol[which(symbol=="")] <- NA
RV <- rownames(feature_data_filt) 
symbol_fill <- coalesce(symbol,RV)

df_genes$RV <- symbol_fill
### Save the list of DE genes in Lipids
write.table(data.frame(genes),file.path(input_dir,"txt/DE_genes_in_Lipids_th<0.05.txt"))
write_xlsx(data.frame(genes),file.path(input_dir,"xlsx/DE_genes_in_Lipids_th<0.05.xlsx"))

write.table(df_genes,file.path(input_dir,"txt/DE_lfc_BH_in_Lipids_th<0.05.txt"))
write_xlsx(df_genes,file.path(input_dir,"xlsx/DE_lfc_BH_in_Lipids_th<0.05.xlsx"))

### Create the vector of genes names and fill the empty spaces with the RV codes
lipid_df_filt$symbol <- feature_data$symbol
symbol <- lipid_df_filt$symbol
symbol[which(symbol=="")] <- NA
RV <- rownames(lipid_df_filt) 
symbol_fill <- coalesce(symbol,RV)
rownames(lipid_df_filt) <- symbol_fill

ymax <- max(-(log10(lipid_df_filt$BH)))

plt <- volcan_plot(lipid_df_filt,x=lipid_df_filt$LogFC ,y=lipid_df_filt$BH, xintercept = 0.0,th=th,ymax=50,title = name )

# pdf(file = file.path(output_dir,"Figures_paper","volcano_plot_DE_lipids_Fig_3B.pdf"),width=10,height=12)
print(plt)
# dev.off()


plt1 <-EnhancedVolcano(lipid_df_filt,x="LogFC",y="BH",lab = rownames(lipid_df_filt),
                       title = "Volcano plot",
                       subtitle = bquote(italic(LCFA ~interaction)),
                       ylab = bquote(~-log[10]~ FDR),
                       legendLabels = c("NS", expression(Log[2] ~ FC), "FDR", expression(~-log[10]~ FDR)),
                       pCutoff = 0.05, FCcutoff = 0,
                       col = c( "grey30","grey30","royalblue","red"),
                       pointSize = 3,
                       colAlpha = 0.9,
                       drawConnectors = T,
                       arrowheads = F,
                       max.overlaps=20)
plt1 <- plt1 +
  # geom_point(size = 3,color="red")+
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5),
        plot.subtitle = element_text(size = 18, hjust = 0.5),
        axis.title = element_text(size = 24),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 2)
  )

plt1

pdf(file = file.path(output_dir,"4_Figures_paper","11_Figure_3B_volcano_plot_DE_lipids.pdf"),width=10,height=12)
print(plt1)
dev.off()
                




