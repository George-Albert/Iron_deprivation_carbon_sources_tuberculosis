#############################
### 0. Load dependencies. ###
#############################
{
  library(readxl)
  library(writexl)
  library(tidyverse)
  library(shadowtext)
  library(xlsx)
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

### Extract the 4 interactions
# 17 Effect_of_Fe_in_response_to_Growth_arrest_G_D.log2FoldChange_shrunken
# 18   Effect_of_Fe_in_response_to_Growth_arrest_G.log2FoldChange_shrunken
# 19 Effect_of_Fe_in_response_to_Growth_arrest_G_L.log2FoldChange_shrunken
# 20   Effect_of_Fe_in_response_to_Growth_arrest_L.log2FoldChange_shrunken

vec_int <- c(17:20)
LogFC_df_int <- LogFC_df[,vec_int]
BH_df_int <- BH_df[,vec_int]
lfcSE_df_int <- lfcSE_df[,vec_int]

colnames(LogFC_df_int) <- paste0(substring(colnames(LogFC_df_int),1,nchar(colnames(LogFC_df_int))-24),".LogFC")
# colnames(BH_df_int) <- paste0(colnames(int_BH),".BH")
colnames(lfcSE_df_int) <- paste0(substring(colnames(lfcSE_df_int),1,nchar(colnames(lfcSE_df_int))-9))

length(which(rownames(LogFC_df_int)!=rownames(BH_df_int)))
length(which(rownames(LogFC_df_int)!=rownames(lfcSE_df_int)))

### Extract the 8 responses to Growth arrest
vec_iron <- c(9:16)
LogFC_df_iron <- LogFC_df[,vec_iron]
BH_df_iron <- BH_df[,vec_iron]
lfcSE_df_iron <- lfcSE_df[,vec_iron]

colnames(LogFC_df_iron) <- paste0(substring(colnames(LogFC_df_iron),1,nchar(colnames(LogFC_df_iron))-24),".LogFC")
# colnames(BH_df_int) <- paste0(colnames(int_BH),".BH")
colnames(lfcSE_df_iron) <- paste0(substring(colnames(lfcSE_df_iron),1,nchar(colnames(lfcSE_df_iron))-9))

LogFC <- data.frame(LogFC_df_iron,LogFC_df_int)
BH <- data.frame(BH_df_iron,BH_df_int)

### Coherence_response just with coherent responses
gene_rows <- c(1:nrow(LogFC))
threshold <- 0.05
threshold_int=0.05
name <- c("C1_C2_Glycerol_Dextrose","C5_C6_Glycerol","C7_C8_Glycerol_LCFA","C11_C12_LCFA")
file_name <- c("enhanced_upregulated.txt","damped_upreg.txt","enhanced_downregulated.txt","damped_downreg.txt")

#################################
#### 4. Absolute Effect sizes ###
#################################
int_summary <- list()
fraction <- c()
fraction_list <- list()
vec <- c(1:4)

for (i in vec) {
  
  indice <- i
  print(paste0("Calculating ",name[indice], " layer"))
  
  summary=data.frame(
    logFC_with=LogFC[,i],
    logFC_without=LogFC[,i+4],
    logFC_int=LogFC[,i+8],
    BH_with=BH[,i],
    BH_without=BH[,i+4],
    BH_int=BH[,i+8])
  rownames(summary)=rownames(LogFC)

  ### Coherent Down direction:
  summary$label_int <- "Background"
  
  summary$label_int[which(summary$logFC_with<0 & summary$logFC_without<0)]="DOWN"
  summary$label_int[which(summary$logFC_with>0 & summary$logFC_without<0 &
                            summary$BH_int<threshold_int)]="Non_Coherent"
  summary$label_int[which(summary$logFC_with<0 & summary$logFC_without>0 &
                            summary$BH_int<threshold_int)]="Non_Coherent"
  
  down_genes=summary[which(summary$label_int %in% c("DOWN")),]
  non_coherent=summary[which(summary$label_int %in% c("Non_Coherent")),]
  
  print(paste("there are", length(down_genes$label_int), "genes downregulated"))
  print(paste("there are", length(non_coherent$label_int), "genes Non_Coherent"))  
  
  # fraction[1] <- 100*nrow(summary[which(summary$label_int %in% c("DOWN")),])/nrow(summary)
  
  ### Enhanced_DOWN:   
  summary$label_int[which( summary$logFC_with<0 & summary$logFC_without<0 & 
                             summary$BH_int<threshold_int & summary$logFC_int<0)]="Enhanced_DOWN"
  
  enhanced_down_genes <- summary[which(summary$label_int %in% c("Enhanced_DOWN")),]
  print(paste("there are", length(enhanced_down_genes$label_int), "induced in the down direction")) 
  
  ### Damped_DOWN:
  summary$label_int[which( summary$logFC_with<0 & summary$logFC_without<0 & 
                             summary$BH_int<threshold_int & summary$logFC_int>0)]="Damped_DOWN"
  damped_down_genes <- summary[which(summary$label_int %in% c("Damped_DOWN")),]
  print(paste("there are", length(damped_down_genes$label_int), "damped in the down direction")) 
  
  down_chunk=summary[which(summary$label_int %in% c("Enhanced_DOWN","Damped_DOWN")),]
  
  down_chunk$int_rebuilt=abs(down_chunk$logFC_without)-abs(down_chunk$logFC_with)
  down_chunk$label_int="Downregulated_genes"

  ### Coherent UP direction:
  summary$label_int[which( summary$logFC_with>0 & summary$logFC_without>0)]="UP"
  up_genes=summary[which(summary$label_int %in% c("UP")),]
  print(paste("there are", length(up_genes$label_int), "genes upregulated"))  
  
  ### Enhanced_UP:   
  summary$label_int[which( summary$logFC_with>0 & summary$logFC_without>0 & 
                             summary$BH_int<threshold_int & summary$logFC_int>0)]="Enhanced_UP"
  enhanced_up_genes <- summary[which(summary$label_int %in% c("Enhanced_UP")),]
  print(paste("there are", length(enhanced_up_genes$label_int), "induced in the up direction")) 
  
  ### Damped_UP:
  summary$label_int[which( summary$logFC_with>0 & summary$logFC_without>0 & 
                             summary$BH_int<threshold_int & summary$logFC_int<0)]="Damped_UP"
  damped_up_genes <- summary[which(summary$label_int %in% c("Damped_UP")),]
  print(paste("there are", length(damped_up_genes$label_int), "damped in the up direction")) 
  
  # summary$label_int[which(summary$logFC_with>0 & summary$logFC_without<0 & summary$BH_int<threshold_int)]="FLIPPED_int_up_down_weak"
  up_chunk=summary[which(summary$label_int %in% c("Enhanced_UP","Damped_UP")),]
  # print(paste0("Enhanced_UP percent is ", fraction[4]))
  
  up_chunk$int_rebuilt=abs(up_chunk$logFC_without)-abs(up_chunk$logFC_with)
  up_chunk$label_int="Upregulated_genes"
  
  ### Create a summary table with the information of how many genes are induced and in what direction  
  summary_df <- data.frame(summary=summary(factor(summary$label_int,levels = c("DOWN","UP",
                                                                               "Damped_DOWN","Damped_UP",
                                                                               "Enhanced_DOWN","Enhanced_UP",
                                                                               "Non_Coherent","Background"))))
  summary_df["Coherent",] <- sum(summary_df[3:6,])
  
  ### Compute the fraction of coherent vs non coherent
  fraction[1] <- 100*nrow(summary[which(summary$label_int %in% c("Non_Coherent")),])/sum(summary_df[3:7,])
  ### Compute each fraction of coherent responses
  fraction[2] <- 100*nrow(summary[which(summary$label_int %in% c("Enhanced_DOWN")),])/sum(summary_df[3:7,])
  fraction[3] <- 100*nrow(summary[which(summary$label_int %in% c("Enhanced_UP")),])/sum(summary_df[3:7,])
  fraction[4] <- 100*nrow(summary[which(summary$label_int %in% c("Damped_DOWN")),])/sum(summary_df[3:7,])
  fraction[5] <- 100*nrow(summary[which(summary$label_int %in% c("Damped_UP")),])/sum(summary_df[3:7,])
  fraction[6] <- sum(fraction[2:5])
  print(paste("the fraction of non coherent genes vs coherent is",fraction[1]))
  print(paste("the fraction of Enhanced_DOWN is",fraction[2]))
  print(paste("the fraction of Enhanced_UP is",fraction[3]))
  print(paste("the fraction of Damped_DOWN is",fraction[4]))
  print(paste("the fraction of Damped_UP is",fraction[5]))
  print(paste("the total fraction of coherent responses is",fraction[6]))
  
  fraction_df <- data.frame(Description=c("non_coh_vs_coh","Enhanced_DOWN","Enhanced_UP","Damped_DOWN",
                                          "Damped_UP","total fraction of coherent responses"),Fraction=fraction)
  fraction_list[[indice]] <- fraction_df
  write.table(fraction_df,file.path(input_dir,paste0("txt/fraction_genes_Fig_3_",name[indice],".txt")))
  
  nombres <- rownames(summary_df)
  summary_df <- data.frame(names=nombres,summary=summary_df)
  colnames(summary_df) <- c(name[indice],paste0(name[indice],"_summary"))
  
  
  write.table(summary_df,file.path(input_dir,paste0("txt/summary_enhanced_genes_",name[indice],".txt")))
  
  int_summary[[indice]] <- summary_df 
}
  
final_summary <- data.frame(int_summary[[1]][2],int_summary[[2]][2],int_summary[[3]][2],int_summary[[4]][2])

write.table(final_summary,file.path(input_dir,paste0("txt/final_summary_enhanced_genes.txt")))
write.table(fraction_list,file.path(input_dir,paste0("txt/final_summary_fraction_enhanced_genes.txt")))
names(fraction_list) <- c("G_D","G","G_L","L")
write_xlsx(fraction_list,file.path(input_dir,paste0("xlsx/final_summary_fraction_enhanced_genes.xlsx")))


# sum_total_genes <- final_summary[1,1]+final_summary[2,1]
# 
# fraction_EU <- 100*final_summary["Enhanced_UP",1]/sum_total_genes
# fraction_ED <- 100*final_summary["Enhanced_DOWN",1]/sum_total_genes
# fraction_DU <- 100*final_summary["Damped_UP",1]/sum_total_genes
# fraction_DD <- 100*final_summary["Damped_DOWN",1]/sum_total_genes
# fraction_non_coherent <- 100*final_summary["Non_Coherent",1]/sum_total_genes



#################################################################
##                       stacked barplot                       ##
#################################################################

categorie <- factor(rownames(final_summary[3:7,]),levels = c("Non_Coherent","Damped_DOWN","Damped_UP",
                                                             "Enhanced_DOWN","Enhanced_UP"))
categorie <- rep(categorie,4)
values <- c(final_summary$C1_C2_Glycerol_Dextrose_summary[3:7],final_summary$C5_C6_Glycerol_summary[3:7],
            final_summary$C7_C8_Glycerol_LCFA_summary[3:7],final_summary$C11_C12_LCFA_summary[3:7])
interactions <- rep(colnames(final_summary),each=5)
interactions <- gsub("^.*?_.*?_", "", interactions)
interactions <- gsub("_summary", "", interactions)
interactions <- factor(interactions,levels = c("LCFA","Glycerol_LCFA","Glycerol","Glycerol_Dextrose"))

df  <- data.frame(interactions,categorie,values)
breaks_values <- c(pretty(final_summary$C11_C12_LCFA_summary))

color <- c("grey","#A59FCB","#FFEFD9","#50486D","#FFA373")

bar_plt <- ggplot(df, aes(fill=categorie, y=values, x=interactions)) + 
  geom_bar(position="stack", stat="identity",width = 0.7, color = "black")+
  coord_flip()+
  scale_y_continuous(expand=c(0,0),breaks = breaks_values,labels = abs(breaks_values),limits=c(0,1850))+
  scale_fill_manual(values = color)+
  theme_classic()+
  geom_hline(yintercept=0)+
  # theme(aspect.ratio = .9)+
  labs(x="Interactions", y = "Number of Genes")+
  theme(
    axis.text.x   = element_text(size=12),
    axis.text.y   = element_text(size=12),
    axis.title.x  = element_text(size=12),
    axis.title.y  = element_text(size=12),
    axis.ticks.y = element_blank(),
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
    data = df[-c(16:19),],
    aes(x=interactions, y=values, label = values),
    hjust = 2,
    position = "stack",
    # nudge_x = -0.3,
    colour = "black",
    bg.colour = "white",
    bg.r = 0.2,
    # family = "Econ Sans Cnd",
    size = 3)

print(bar_pl_1)

pdf(file = file.path(output_dir,"4_Figures_paper","12_Figure_3C_Stacked_Bar_plot_interactions_enhanced.pdf"),width=8,height=5)
print(bar_pl_1)
dev.off()

