#############################
### 0. Load dependencies. ###
#############################

{
  library(edgeR)
  library(qvalue)
  library(limma)
  library(ggplot2)
  library(cowplot)
  library(reshape2)
  library(DESeq2)
  library(fdrtool)
  # library(ashr)
  library(dplyr)
  library(locfdr)
  library(xlsx)
  library(writexl)
  library(stringr)
}

############################
### 1. Declare functions ### 
############################
{
  dcols=function(x){data.frame(colnames(x))}
  
  my_name <- function(v1) {
    deparse(substitute(v1))
  }
  
  histogram <- function(data,pvalue,title){
    frame <- data.frame(data)
    hist <- ggplot(frame)+geom_histogram(aes(x=pvalue),color="black",
                                         fill="lightblue",bins = (n_genes/(n_genes/50)),boundary=0)+
      scale_x_continuous(limits = c(0, 1)) +
      theme_minimal()+ ggtitle(wrapper(title, width = 20))+labs(x="PValues", y = "Counts")+
      theme(plot.title = element_text(hjust = 0.5))
    return(hist)
  }
  
  process_contrast_DEseq=function(dds,contrast,th=th,res,title){
    
    # name <-deparse(substitute(res))
    
    result_a <- results(dds, contrast=contrast, pAdjustMethod = "BH")
    colnames(result_a)[c(2,3)]=paste0(colnames(result_a)[c(2,3)],"_raw")
    result_b <- lfcShrink(dds, contrast=contrast, type="ashr")
    colnames(result_b)[c(2,3)]=paste0(colnames(result_b)[c(2,3)],"_shrunken")
    
    result=cbind(result_a,result_b[,c("log2FoldChange_shrunken", "lfcSE_shrunken")])
    colnames(result)[colnames(result) == "padj"] <- "BH"
    result$BH <- p.adjust(result$pvalue, method = "BH")
    result=result[,c("baseMean","log2FoldChange_raw","lfcSE_raw","log2FoldChange_shrunken",
                     "lfcSE_shrunken","stat","pvalue","BH")]
    print(paste(res," With BH <", th," there are ",length(which(result$BH<th))," genes"))
    
    # sum(result$BH < 0.1, na.rm=TRUE)
    nbins=50
    pl=histogram(result,result$pvalue,title) 
    # pl=ggplot(result)+geom_histogram(aes(x=pvalue),bins=nbins,center = 1/(2*nbins))+ggtitle(name)
    output=list(data=result,figure=pl)
    return(output)
  }
  
  write_results=function(tab,name,fdata=feature_data,th=th,th_size=0){
    
    if(length(which(rownames(fdata)!=rownames(tab)))>0){
      print("tablas no congruentes")
    } else {
      
      dir.create(paste0("Analysis/Outputs/Data/xlsx/th_",th,"_th_size_",th_size,"/",name),showWarnings = FALSE,recursive=TRUE)
      dir.create(paste0("Analysis/Outputs/Data/txt/th_",th,"_th_size_",th_size),showWarnings = FALSE,recursive=TRUE)
      
      #tab$t=tab$log2FoldChange_raw/tab$lfcSE_raw
      tab$RV=rownames(tab)
      tab$Gene=fdata$symbol
      tab=tab[order(tolower(tab$RV)),c(9,10,1:8)]
      write.table(tab,paste0("Analysis/Outputs/Data/txt/th_",th,"_th_size_",th_size,"/",name,".txt"))
      
      tab=tab[order(-tab$stat),]
      
      upreg=tab[which(tab$BH<th & tab$log2FoldChange_shrunken>th_size),]
      downreg=tab[which(tab$BH<th & tab$log2FoldChange_shrunken<(-th_size)),]
      downreg=downreg[order(downreg$BH),]
      
      write_xlsx(tab,paste0("Analysis/Outputs/Data/xlsx/th_",th,"_th_size_",th_size,"/",name,"/all.xlsx"))
      write_xlsx(upreg,paste0("Analysis/Outputs/Data/xlsx/th_",th,"_th_size_",th_size,"/",name,"/upreg.xlsx"))
      write_xlsx(downreg,paste0("Analysis/Outputs/Data/xlsx/th_",th,"_th_size_",th_size,"/",name,"/downreg.xlsx"))
      
      output=list(up=upreg,down=downreg)
      
      return(output)}
    
  }
  
 
   target_genesets=function(tab,name){
    for(th in c(0.05,0.01))
    {
      for(th_size in c(0,0.2,0.5,1))
      {
        dir.create(paste0("Analysis/Outputs/Data/target_sets/th_",th,"_th_size_",th_size,"/",name),showWarnings = FALSE,recursive=TRUE)
        
        set=(which(abs(tab$log2FoldChange_shrunken)>th_size & tab$BH<th))
        print(paste("All: th=",th,"th_size=",th_size,":",length(set)))
        target=rownames(tab)[set]
        sink(paste0("Analysis/Outputs/Data/target_sets/th_",th,"_th_size_",th_size,"/",name,"/all.txt"))
        for(gen in target)
          cat(gen,"\n")
        sink()
        
        set=(which((tab$log2FoldChange_shrunken)>th_size & tab$BH<th))
        print(paste("UP: th=",th,"th_size=",th_size,":",length(set)))
        target=rownames(tab)[set]
        sink(paste0("Analysis/Outputs/Data/target_sets/th_",th,"_th_size_",th_size,"/",name,"/upreg.txt"))
        for(gen in target)
          cat(gen,"\n")
        sink()
        
        set=(which((tab$log2FoldChange_shrunken)<(-th_size) & tab$BH<th))
        print(paste("DOWN: th=",th,"th_size=",th_size,":",length(set)))
        target=rownames(tab)[set]
        sink(paste0("Analysis/Outputs/Data/target_sets/th_",th,"_th_size_",th_size,"/",name,"/downreg.txt"))
        for(gen in target)
          cat(gen,"\n")
        sink()
      }
    }
  }
  
  wrapper <- function(x, ...) 
  {
    paste(strwrap(x, ...), collapse = "\n")
  }
  
  get_voom_like_normalized_data=function(desq){
    cuentas=counts(desq)
    normalization_coeficients=dds@colData$sizeFactor
    depths=apply(cuentas,2,sum)+1
    cuentas=cuentas+0.5
    
    normalized_counts=cuentas
    for(i in 1:ncol(normalized_counts))
      normalized_counts[,i]=normalized_counts[,i]/(depths[i]*normalization_coeficients[i]/1E6)
    
    normalized_counts=log2(normalized_counts)
    return(normalized_counts)
  }
  
}

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
metadata<- read.table(file.path(input_dir,"txt/metadata_32_samples.txt"),check.names = F)
reads <- read.table(file.path(input_dir,"txt/reads_32_samples.txt"),check.names = F)
feature_data <- read.table(file.path(input_dir,"txt/feature_data_genes.txt"),check.names = F)

rownames(feature_data) <- feature_data$locus_tag
feature_data <- feature_data[rownames(reads),]

length(which(rownames(feature_data)!=(rownames(reads))))
length(which(colnames(reads)!=rownames(metadata)))

### Extract the rRNA genes
feature_data <- feature_data[which(feature_data$class!="rRNA"),]
reads <- reads[rownames(feature_data),]
length(which(rownames(feature_data)!=(rownames(reads))))

#############################
### 4. Select sampleset:  ### 
#############################

### number of samples and genes
n_samples <- ncol(reads)
n_genes <- nrow (reads)
n_Fe_NO <- nrow(metadata[which(metadata$Iron=="NO"),])
n_Fe_Yes <- n_samples-n_Fe_NO
th <- 0.05

group <- factor(metadata$short_setup)

### Design matrix 
design=model.matrix(~0+group,data=factor(metadata$short_setup))
colnames(design) <- levels(group)
contrast_vec <- colnames(design)

### Contrast matrix
{
  ##############################################################################
  ############################### 5. Contrasts #################################
  ##############################################################################
  
  #=========================== IRON Effects in EXP with different carbon sources ===============================
  #Glycerol, Dextrose
  C1vsC1 <- makeContrasts((C1_Fe_NO-C1_Fe_YES), levels=design)
  #Glycerol 
  C5vsC5 <- makeContrasts((C5_Fe_NO-C5_Fe_YES), levels=design)
  #Glycerol, LCFA
  C7vsC7 <- makeContrasts((C7_Fe_NO-C7_Fe_YES), levels=design)
  #LCFA
  C11vsC11 <- makeContrasts((C11_Fe_NO-C11_Fe_YES), levels=design)
  
  iron_contrast <- makeContrasts(C1vsC1,C5vsC5,C7vsC7,C11vsC11, levels=design)
  
  #=========================== IRON Effects in STAT with different carbon sources ===============================
  #Glycerol, Dextrose
  C2vsC2 <- makeContrasts((C2_Fe_NO-C2_Fe_YES), levels=design)
  #Glycerol 
  C6vsC6 <- makeContrasts((C6_Fe_NO-C6_Fe_YES), levels=design)
  #Glycerol, LCFA
  C8vsC8 <- makeContrasts((C8_Fe_NO-C8_Fe_YES), levels=design)
  #LCFA
  C12vsC12 <- makeContrasts((C12_Fe_NO-C12_Fe_YES), levels=design)
  
  iron_contrast_stat <- makeContrasts(C2vsC2,C6vsC6,C8vsC8,C12vsC12, levels=design)
  
  #=========================== PHASE Effects +Fe for different carbon sources ===============================
  
  #Glycerol, Dextrose
  C1vsC2 <- makeContrasts((C2_Fe_YES-C1_Fe_YES), levels=design)
  #Glycerol 
  C5vsC6 <- makeContrasts((C6_Fe_YES-C5_Fe_YES), levels=design)
  #Glycerol, LCFA
  C7vsC8 <- makeContrasts((C8_Fe_YES-C7_Fe_YES), levels=design)
  #LCFA
  C11vsC12 <- makeContrasts((C12_Fe_YES-C11_Fe_YES), levels=design)
  
  phase_contrast_iron <- makeContrasts(C1vsC2,C5vsC6,C7vsC8,C11vsC12, levels=design)
  
  #=========================== PHASE Effects -Fe for different carbon sources ===============================
  
  #Glycerol, Dextrose
  C1vsC2_Fe_NO <- makeContrasts((C2_Fe_NO-C1_Fe_NO), levels=design)
  #Glycerol 
  C5vsC6_Fe_NO  <- makeContrasts((C6_Fe_NO-C5_Fe_NO), levels=design)
  #Glycerol, LCFA
  C7vsC8_Fe_NO  <- makeContrasts((C8_Fe_NO-C7_Fe_NO), levels=design)
  #LCFA
  C11vsC12_Fe_NO  <- makeContrasts((C12_Fe_NO-C11_Fe_NO), levels=design)
  
  phase_contrast_NO <- makeContrasts(C1vsC2_Fe_NO, C5vsC6_Fe_NO, C7vsC8_Fe_NO, C11vsC12_Fe_NO, levels=design)
  
  #=========================== CARBON SOURCES Effects +Fe in EXP ===============================
  
  #Glycerol, Dextrose ---> Glycerol
  C1vsC5 <- makeContrasts((C5_Fe_YES-C1_Fe_YES), levels=design)
  #Glycerol ---> Glycerol, LCFA
  C5vsC7  <- makeContrasts((C7_Fe_YES-C5_Fe_YES), levels=design)
  #Glycerol, LCFA---> LCFA
  C7vsC11  <- makeContrasts((C11_Fe_YES-C7_Fe_YES), levels=design)
  #Glycerol, Dextrose---> LCFA
  C1vsC11 <- makeContrasts((C11_Fe_YES-C1_Fe_YES), levels=design)
  
  
  carb_contrast_exp <- makeContrasts(C1vsC5,C5vsC7,C7vsC11,C1vsC11, levels=design)
  
  #=========================== CARBON SOURCES Effects -Fe in EXP ===============================
  
  #Glycerol, Dextrose ---> Glycerol
  C1vsC5_Fe_NO <- makeContrasts((C5_Fe_NO-C1_Fe_NO), levels=design)
  #Glycerol ---> Glycerol, LCFA
  C5vsC7_Fe_NO  <- makeContrasts((C7_Fe_NO-C5_Fe_NO), levels=design)
  #Glycerol, LCFA---> LCFA
  C7vsC11_Fe_NO  <- makeContrasts((C11_Fe_NO-C7_Fe_NO), levels=design)
  
  
  carb_contrast_NO_exp <- makeContrasts(C1vsC5_Fe_NO, C5vsC7_Fe_NO, C7vsC11_Fe_NO, levels=design)
  
  #=========================== CARBON SOURCES Effects +Fe in STAT ===============================
  
  #Glycerol, Dextrose ---> Glycerol
  C2vsC6 <- makeContrasts((C6_Fe_YES-C2_Fe_YES), levels=design)
  #Glycerol ---> Glycerol, LCFA
  C6vsC8  <- makeContrasts((C8_Fe_YES-C6_Fe_YES), levels=design)
  #Glycerol, LCFA---> LCFA
  C8vsC12  <- makeContrasts((C12_Fe_YES-C8_Fe_YES), levels=design)
  #Glycerol, Dextrose---> LCFA
  C2vsC12 <- makeContrasts((C12_Fe_YES-C2_Fe_YES), levels=design)
  
  carb_contrast_stat <- makeContrasts(C2vsC6, C6vsC8, C8vsC12,C2vsC12, levels=design)
  
  #=========================== CARBON SOURCES Effects -Fe in STAT ===============================
  
  #Glycerol, Dextrose ---> Glycerol
  C2vsC6_Fe_NO <- makeContrasts((C2_Fe_NO-C6_Fe_NO), levels=design)
  #Glycerol ---> Glycerol, LCFA
  C6vsC8_Fe_NO  <- makeContrasts((C6_Fe_NO-C8_Fe_NO), levels=design)
  #Glycerol, LCFA---> LCFA
  C8vsC12_Fe_NO  <- makeContrasts((C8_Fe_NO-C12_Fe_NO), levels=design)
  
  
  carb_contrast_NO_stat <- makeContrasts(C2vsC6_Fe_NO, C6vsC8_Fe_NO, C8vsC12_Fe_NO, levels=design)
  
  #=========================== 2nd order interactions ===============================
  
  # =========================== Iron Response Effects ===============================
  # Glycerol, Dextrose
  C1C2vsC1C2 <- makeContrasts(C1vsC2_Fe_NO - C1vsC2, levels=design)
  
  # Glycerol
  C5C6vsC5C6 <- makeContrasts(C5vsC6_Fe_NO - C5vsC6, levels=design)
  
  # Glycerol,LCFA
  C7C8vsC7C8 <- makeContrasts(C7vsC8_Fe_NO - C7vsC8, levels=design)
  
  # LCFA
  C11C12vsC11C12 <- makeContrasts(C11vsC12_Fe_NO - C11vsC12, levels=design)
  # =========================== Carbon Source Effects in presence of Iron (Fe+) ===============================
  #Glycerol, Dextrose ---> Glycerol, Fe+
  C1C2vsC5C6 <- makeContrasts(C5vsC6-C1vsC2,levels=design)
  
  #Glycerol ---> Glycerol,LCFA Fe+
  C5C6vsC7C8 <- makeContrasts(C7vsC8-C5vsC6,levels=design)
  
  #Glycerol, LCFA ---> LCFA Fe+
  C7C8vsC11C12 <- makeContrasts(C11vsC12-C7vsC8,levels=design)
  
  #Glycerol, Dextrose ---> LCFA, Fe+ (Growth arrest effect)
  C1C2vsC11C12 <- makeContrasts(C11vsC12-C1vsC2,levels=design)
  # =========================== Carbon Source Effects in absence of Iron (Fe-) ===============================
  #Glycerol, Dextrose ---> Glycerol, Fe+
  C1C2vsC5C6_Fe_NO <- makeContrasts(C5vsC6_Fe_NO-C1vsC2_Fe_NO,levels=design)
  
  #Glycerol ---> Glycerol,LCFA Fe+
  C5C6vsC7C8_Fe_NO <- makeContrasts(C7vsC8_Fe_NO-C5vsC6_Fe_NO,levels=design)
  
  #Glycerol, LCFA ---> LCFA Fe+
  C7C8vsC11C12_Fe_NO <- makeContrasts(C11vsC12_Fe_NO-C7vsC8_Fe_NO,levels=design)
  
  #Glycerol, Dextrose ---> LCFA, Fe+ (Growth arrest effect)
  C1C2_Fe_NOvsC11C12_Fe_NO <- makeContrasts(C11vsC12_Fe_NO-C1vsC2_Fe_NO,levels=design)
  
  # (Growth arrest effect diagonal)
  C1C2vsC11C12_Fe_NO <- makeContrasts(C11vsC12_Fe_NO-C1vsC2,levels=design)
  
  second_order_contrast <- makeContrasts(C1C2vsC1C2,C5C6vsC5C6,C7C8vsC7C8,C11C12vsC11C12,C1C2vsC5C6,
                                         C5C6vsC7C8,C7C8vsC11C12,C1C2vsC11C12,C1C2vsC5C6_Fe_NO,C5C6vsC7C8_Fe_NO,
                                         C7C8vsC11C12_Fe_NO,C1C2_Fe_NOvsC11C12_Fe_NO,C1C2vsC11C12_Fe_NO, levels=design)
  # =========================== Effect of Fe in response to growth arrest ===============================

  C1C1vsC2C2 <- makeContrasts(C2vsC2-C1vsC1,levels=design)
  C5C5vsC6C6 <- makeContrasts(C6vsC6-C5vsC5,levels=design)
  C7C7vsC8C8 <- makeContrasts(C8vsC8-C7vsC7,levels=design)
  C11C11vsC12C12 <- makeContrasts(C12vsC12-C11vsC11,levels=design)
 
  second_order_contrast <- makeContrasts(C1C2vsC1C2,C5C6vsC5C6,C7C8vsC7C8,C11C12vsC11C12,C1C2vsC5C6,
                                         C5C6vsC7C8,C7C8vsC11C12,C1C2vsC11C12,C1C2vsC5C6_Fe_NO,C5C6vsC7C8_Fe_NO,
                                         C7C8vsC11C12_Fe_NO,C1C2_Fe_NOvsC11C12_Fe_NO,C1C2vsC11C12_Fe_NO,C1C1vsC2C2,
                                         C5C5vsC6C6,C7C7vsC8C8,C11C11vsC12C12,levels=design)
  
   # =========================== Diagonal Effects (Fe-) ===============================
  #Glycerol, Dextrose ---> Glycerol, Fe+
  C1vsC11_Fe_NO <- makeContrasts(C11_Fe_NO-C1_Fe_YES,levels=design)
  
  #Glycerol ---> Glycerol,LCFA Fe+
  C2vsC12_Fe_NO <- makeContrasts(C12_Fe_NO-C2_Fe_YES,levels=design)
  
  diagonal_contrast <- makeContrasts(C1vsC11_Fe_NO,C2vsC12_Fe_NO,levels=design)
}

contrast_matrix <- cbind(iron_contrast,iron_contrast_stat,phase_contrast_iron,phase_contrast_NO,carb_contrast_exp,carb_contrast_NO_exp,
                         carb_contrast_stat,carb_contrast_NO_stat,second_order_contrast,diagonal_contrast)

write.table(contrast_matrix,file.path(input_dir,"txt/contrast_matrix.txt"))

##################
### 6. DESeq2  ### 
##################


# design = ~ 0 + short_setup 
dds <- DESeqDataSetFromMatrix(countData = round(reads),colData = metadata,design = ~ 0 + short_setup)
dds <- DESeq(dds,test = "Wald",full = design, betaPrior = F)
print(attr(dds, "modelMatrix"))
# keep <- filterByExpr(dds,design = design, group = group)
keep <- filterByExpr(dds)
filtered_dds <- dds[keep, ]
resultsNames(filtered_dds)   

feature_data <- feature_data[rownames(filtered_dds),]
length(which(rownames(feature_data)!=rownames(filtered_dds)))
write.table(feature_data,file.path(input_dir,"txt/feature_data_filtered.txt"))

##############################################
# Obtain the normalized counts for the heatmap
##############################################

normalized_dds=get_voom_like_normalized_data(filtered_dds)

medias=normalized_dds[,1:2]

set_no=which(metadata$Iron=="NO")
set_yes=which(metadata$Iron=="YES")

medias[,1]=apply(normalized_dds[,set_no],1,mean)
medias[,2]=apply(normalized_dds[,set_yes],1,mean)

for(i in set_no){
  normalized_dds[,i]=normalized_dds[,i]-medias[,1]
}
for(i in set_yes){
  normalized_dds[,i]=normalized_dds[,i]-medias[,2]
}

##############################################
contrast_res <- colnames(contrast_matrix)
titles <- c("Iron_effect_in_EXP_(GLYCEROL-DEXTROSE)","Iron_effect_in_EXP_(GLYCEROL)",
            "Iron_effect_in_EXP_(GLYCEROL-LCFA)","Iron_effect_in_EXP_(LCFA)",
            
            "Iron_effect_in_STAT_(GLYCEROL-DEXTROSE)","Iron_effect_in_STAT_(GLYCEROL)",
            "Iron_effect_in_STAT_(GLYCEROL-LCFA)","Iron_effect_in_STAT_(LCFA)",
            
            "Phase_effects_with_Fe_(GLYCEROL-DEXTROSE)","Phase_effects_with_Fe_(GLYCEROL)",
            "Phase_effects_with_Fe_(GLYCEROL-LCFA)","Phase_effects_with_Fe_(LCFA)",
            
            "Phase_effects_without_Fe_(GLYCEROL-DEXTROSE)","Phase_effects_without_Fe_(GLYCEROL)",
            "Phase_effects_without_Fe_(GLYCEROL-LCFA)","Phase_effects_without_Fe_(LCFA)",
            
            "Carbon_source_effects_with_Fe_in_EXP_(GLYCEROL-DEXTROSE -> GLYCEROL)",
            "Carbon_source_effects_with_Fe_in_EXP_(GLYCEROL -> GLYCEROL,LCFA)",
            "Carbon_source_effects_with_Fe_in_EXP_(GLYCEROL,LCFA -> LCFA)",
            "Carbon_source_effects_with_Fe_in_EXP_(GLYCEROL-DEXTROSE -> LCFA)",
            
            "Carbon_source_effects_without_Fe_in_EXP_(GLYCEROL-DEXTROSE-> GLYCEROL)",
            "Carbon_source_effects_without_Fe_in_EXP_(GLYCEROL -> GLYCEROL,LCFA)",
            "Carbon_source_effects_without_Fe_in_EXP_(GLYCEROL,LCFA -> LCFA)",
            
            "Carbon_source_effects_with_Fe_in_STAT_(GLYCEROL-DEXTROSE -> GLYCEROL)",
            "Carbon_source_effects_with_Fe_in_STAT_(GLYCEROL -> GLYCEROL,LCFA)",
            "Carbon_source_effects_with_Fe_in_STAT_(GLYCEROL,LCFA -> LCFA)",
            "Carbon_source_effects_with_Fe_in_STAT_(GLYCEROL-DEXTROSE -> LCFA)",
            
            "Carbon_source_effects_without_Fe_in_STAT_(GLYCEROL-DEXTROSE -> GLYCEROL)",
            "Carbon_source_effects_without_Fe_in_STAT_(GLYCEROL-> GLYCEROL,LCFA)",
            "Carbon_source_effects_without_Fe_in_STAT_(GLYCEROL,LCFA -> LCFA)",
            
            "Iron_Response_(GLYCEROL-DEXTROSE)","Iron_Response_(GLYCEROL)",
            "Iron_Response_(GLYCEROL-LCFA)","Iron_Response_(LCFA)",
            
            "Carbon_source_effects_in_presence_of_Fe_(GLYCEROL-DEXTROSE -> GLYCEROL)",
            "Carbon_source_effects_in_presence_of_Fe_(GLYCEROL-> GLYCEROL,LCFA)",
            "Carbon_source_effects_in_presence_of_Fe_(GLYCEROL,LCFA -> LCFA)",
            "Carbon_source_effects_in_presence_of_Fe_(GLYCEROL-DEXTROSE -> LCFA)",
            
            "Carbon_source_effects_in_absence_of_Fe_(GLYCEROL-DEXTROSE -> GLYCEROL)",
            "Carbon_source_effects_in_absence_of_Fe_(GLYCEROL-> GLYCEROL,LCFA)",
            "Carbon_source_effects_in_absence_of_Fe_(GLYCEROL,LCFA -> LCFA)",
            "Carbon_source_effects_in_absence_of_Fe_(GLYCEROL-DEXTROSE -> LCFA)",
            
            "Growth_arrest_diagonal",
            "Effect_of_Fe_in_response_to_Growth_arrest_(GLYCEROL-DEXTROSE)",
            "Effect_of_Fe_in_response_to_Growth_arrest_(GLYCEROL)",
            "Effect_of_Fe_in_response_to_Growth_arrest_(GLYCEROL_LCFA)",
            "Effect_of_Fe_in_response_to_Growth_arrest_(LCFA)",
            
            "Iron_x_carbon_source_effects_in EXP",
            "Iron_x_carbon_source_effects_in STAT")

df_name_contrast <- data.frame(title=titles,contrasts=contrast_res)
write.table(df_name_contrast,file.path(input_dir,"txt/contrasts_nomenclature.txt"))
write_xlsx(df_name_contrast,file.path(input_dir,"xlsx/contrasts_nomenclature.xlsx"))

y <- vector("list", length(contrast_res))
up <- vector("list", length(contrast_res))
down <- vector("list", length(contrast_res))

y <- setNames(y, contrast_res)
up<- setNames(up, contrast_res)
down <- setNames(down, contrast_res)
th=0.05

for (res in contrast_res) {
  
  print(res)
  # name <-deparse(substitute(contrast))
  title <- df_name_contrast[which(df_name_contrast$contrasts==res),1]
  list_res=process_contrast_DEseq(filtered_dds,th=th,contrast=contrast_matrix[,res],
                                  res=res,title = title)
  result <- list_res$data
  # feature_data <- feature_data[rownames(result),]
  # length(which(rownames(feature_data)!=(rownames(result))))
  
  x <-   print(paste(res," With BH <", th," there are ",length(which(result$BH<th))," genes"))
  y[res] <- x
  
  pl_hist <- list_res$figure
  print(pl_hist)
  
  dir <- "4_Figures/1_signal/histograms"
  dir.create(file.path(output_dir,dir),recursive=TRUE,showWarnings = F)
  pdf(file.path(file.path(output_dir,dir),paste0(title,".pdf")),width=6,height=5)
  print(pl_hist)
  dev.off()
  
  print("before write results everything ok")
  list_DE <- write_results(tab=data.frame(result),name=title,fdata=feature_data,th=th,th_size=0)
  print("after write results everything ok")
  
  up[res] <- dim(list_DE$up)[1]
  down[res] <- dim(list_DE$down)[1]
  # up <- do.call("cbind",list_DE$up)
  # down <- do.call("cbind",list_DE$down)
  
  target_genesets(list_res$data,name=title)
  
}

y <- do.call("rbind",y)
up <- do.call("rbind",up)
down <- do.call("rbind",down)

y <- data.frame(y)
up <- data.frame(up)
down <- data.frame(down)

colnames(y) <- paste0("genes with BH <",th)
# str_extract(y[1], "(?<=are).*(?= genes)")

for (i in c(1:dim(y)[1])){
  
  y[i,1] <- str_extract(y[i,1], "(?<=are).*(?= genes)")
  
}

reg <- merge(up,down,by=0,sort=F)
rownames(reg) <- reg$Row.names
reg$Row.names <- NULL
y <- merge(reg,y, by=0,sort=F)
rownames(y) <- y$Row.names
y <-y[c(paste0("genes with BH <",th),"up","down")] 
colnames(y) <- c("DE","upreg","downreg")
DE_genes <- y

write.table(DE_genes,file.path(input_dir,"txt/DE_genes_per_cond_th<0.05.txt"))
write_xlsx(DE_genes,file.path(input_dir,"xlsx/DE_genes_per_cond_th<0.05.xlsx"))

# DE_genes_of_interest <- DE_genes[c(1:16,44:47),]
DE_genes_of_interest <- DE_genes[c(9:16,44:47),]

write.table(DE_genes_of_interest,file.path(input_dir,"txt/DE_genes_of_interest_th<0.05.txt"))
write_xlsx(DE_genes_of_interest,file.path(input_dir,"xlsx/DE_genes_of_interest_th<0.05.xlsx"))

# "Contrast" "hits"            
# "1" "Fe_at_C1" 0
# "2" "Fe_at_C2" 1741
# "3" "Fe_at_C5" 0
# "4" "Fe_at_C6" 1386
# "5" "Fe_at_C7" 8
# "6" "Fe_at_C8" 2050
# "7" "Fe_at_C11" 8
# "8" "Fe_at_C12" 109
# "15" "Int_g" 1101
# "16" "Int_h" 1028
# "17" "Int_i" 1569
# "18" "Int_j" 83
# 
#                     DE upreg downreg
# C1vsC1              0      0       0
# C5vsC5              0      0       0
# C7vsC7              7      6       1
# C11vsC11            9      9       0
# C2vsC2           1370    753     617
# C6vsC6           1234    612     622
# C8vsC8           1790    808     982
# C12vsC12          115     70      45
# C1C1vsC2C2        769    458     311
# C5C5vsC6C6        847    440     407
# C7C7vsC8C8       1274    557     717
# C11C11vsC12C12     81     59      22
