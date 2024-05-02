
#############################
### 0. Load dependencies. ###
#############################
{ library(tidyverse)
  library(readxl)
  library(writexl)
  library(shadowtext)
}

############################
### 1. Declare functions ### 
############################
dcols=function(x){data.frame(colnames(x))}

###################################################
### 2. Set working directory and create folders ###
###################################################

main_wd <- getwd()
input_dir <- "Analysis/Inputs/2_Processed_data"
output_dir <- "Analysis/Outputs"

######################
#### 3. Load data ####
######################
myData <- readRDS(file.path(input_dir,"txt/Contrasts_stat.RDS"))

LogFC_df <- read.table(file.path(input_dir,"txt/LogFC_0.05.txt"))
BH_df <- read.table(file.path(input_dir,"txt/BH_0.05.txt"))
lfcSE_df <- read.table(file.path(input_dir,"txt/lfcSE_0.05.txt"))

# # Define the pattern to remove the parenthesis 
# names(myData)
# pattern <- paste(c("[(]", "[)]"), collapse = "|")

### Extract the 8 responses to Growth arrest
vec_iron <- c(9:16)
LogFC_df_iron <- LogFC_df[,vec_iron]
BH_df_iron <- BH_df[,vec_iron]
lfcSE_df_iron <- lfcSE_df[,vec_iron]

colnames(LogFC_df_iron) <- paste0(str_extract(colnames(LogFC_df_iron),"(?<=Fe_).*?(?=\\.)"),
                                  rep(c("_with","_without"),each=4))
colnames(BH_df_iron) <- paste0(colnames(LogFC_df_iron),".BH")
colnames(lfcSE_df_iron) <- paste0(colnames(LogFC_df_iron),".lfcSE")

iron_response_to_stat <- data.frame(LogFC_df_iron,BH_df_iron)

write.table(iron_response_to_stat,file.path(input_dir,"txt/iron_response_to_stat.txt"))

th=0.05 # threshold of the fdr

### pair combinations to correlate
combination_names_with <- combn(colnames(iron_response_to_stat)[1:4],m=2)
combination_names_without <- combn(colnames(iron_response_to_stat)[5:8],m=2)

### select the samples with or without iron
vector_to_loop <- seq_along(combination_names_with[1:ncol(combination_names_with)])

############################################################################
############################################################################
###                                                                      ###
###           4. GROWTH ARREST EFFECTS VS GROWTH ARREST C5VSC6           ###
###                                                                      ###
############################################################################
############################################################################

name_of_contrasts <- df_name_contrast[9:16,1]
df_phase_effects <- LogFC_df_iron
df_iron_effects <- iron_response_to_stat
colnames(df_iron_effects) <- c(name_of_contrasts,paste0(name_of_contrasts,".BH"))

color <- c("red","mediumblue","orange")

vec=c(1:3)

### pair combinations to correlate
combination_names_with <- combn(colnames(iron_response_to_stat)[1:4],m=1)
combination_names_without <- combn(colnames(iron_response_to_stat)[5:8],m=2)

# Select the element to combine with 
element_to_combine <- colnames(iron_response_to_stat)[2]
element_to_combine_wo <- colnames(iron_response_to_stat)[6]
# Create a new vector excluding the element to combine with
iron_response_to_stat_wo_element <- colnames(iron_response_to_stat)[-2]
iron_response_to_stat_wo_element_2 <- colnames(iron_response_to_stat)[-6]
# Generate all the combinations with the element we selected
combination_names_with <- expand.grid(element_to_combine, iron_response_to_stat_wo_element[1:3])
combination_names_with <- as.matrix(t(combination_names_with))

combination_names_without <- expand.grid(element_to_combine_wo, iron_response_to_stat_wo_element_2[5:7])
combination_names_without <- as.matrix(t(combination_names_without))

# Show the result
print(combination_names_with)

combination_names <- combination_names_without
### select the samples with or without iron
vector_to_loop <- seq_along(combination_names[1:ncol(combination_names)])

for (col in vector_to_loop) {
  
  print(paste(combination_names[1,col],combination_names[2,col],
              combination_names[1,col],".BH",
              combination_names[2,col],".BH"))
  
  df <- data.frame(LogFC1=iron_response_to_stat[,combination_names[1,col]],
                   LogFC2=iron_response_to_stat[,combination_names[2,col]],
                   BH1=iron_response_to_stat[,paste0(combination_names[1,col],".BH")],
                   BH2=iron_response_to_stat[,paste0(combination_names[2,col],".BH")],
                   row.names = rownames(iron_response_to_stat))
  
  
  df_filtered <-df[which(df$BH1 < th & df$BH2 < th),] 
  
  ###########################################################################
  ###########################################################################
  ###                                                                     ###
  ###                            SCATTER PLOTS                            ###
  ###                                                                     ###
  ###########################################################################
  ###########################################################################
  
  ### Define axis
  y <- df_filtered$LogFC2
  x <- df_filtered$LogFC1

  ### Define axis title
  x_title <- paste0("Growth Arrest Effects (" ,combination_names[2,col]," Fe)")
  y_title <- paste0("Growth Arrest Effects (" ,combination_names[1,col]," Fe)")
  
  ### Pearson Correlation
  test=cor.test(y,x,method = "pearson")
  print(test)
  
  ### plot
  pl_scatter=ggplot(df_filtered)+
    geom_point(aes(x=LogFC1,y=LogFC2),color=color[col],pch=21,alpha=0.8)+
    # scale_fill_manual(values = color)+
    # geom_abline(slope=0,intercept=0,color="black")+
    xlab(x_title)+
    ylab(y_title)+
    theme(legend.position.inside=c(0.8,0.2),
          panel.background = element_rect(fill="grey90"),
          panel.grid.major = element_line(colour = "grey95"),
          panel.grid.minor = element_line(colour = "grey95")
    )
  
  pl_scatter_1 <- pl_scatter+
    geom_smooth(aes(x=LogFC1,y=LogFC2),color="black",method=lm,formula = "y ~ x",show.legend = FALSE)+
    # geom_rug(aes(x=LogFC1,y=LogFC2),color=color[i])+
    geom_vline(xintercept = 0)+
    geom_hline(yintercept = 0)+
    ggtitle(paste0("R=",round(test$estimate,digits=2)," Num_of_common genes=",dim(df_filtered)[1]))+
    theme_minimal()
  
    print(pl_scatter_1)
   
    fig_repo <- "4_Figures_paper"
    dir.create(file.path(output_dir,fig_repo,"4_Figure_S1B_Scatter_plot_correlogram"),recursive = T,showWarnings = F)
  
    pdf(file.path(output_dir,fig_repo,"4_Figure_S1B_Scatter_plot_correlogram",paste0(x_title,"_vs_",y_title,".pdf")),width=9,height=7)
    print(pl_scatter_1)
    dev.off()
}

#################################################################
##                        Same C-Source                        ##
#################################################################

combination_names <-colnames(iron_response_to_stat[1:8])
color <- c("red","mediumblue","orange","cyan")
### select the samples with or without iron
vector_to_loop <- c(1:4)

for (col in vector_to_loop) {
  
  print(paste(combination_names[col],combination_names[col+4],
              combination_names[col],".BH",
              combination_names[col+4],".BH"))
  
  df <- data.frame(LogFC1=iron_response_to_stat[,combination_names[col]],
                   LogFC2=iron_response_to_stat[,combination_names[col+4]],
                   BH1=iron_response_to_stat[,paste0(combination_names[col],".BH")],
                   BH2=iron_response_to_stat[,paste0(combination_names[col+4],".BH")],
                   row.names = rownames(iron_response_to_stat))
  
  
  df_filtered <-df[which(df$BH1 < th & df$BH2 < th),] 
  
  ####################
  ###Scatter plots ###
  ####################
  ### Define axis
  y <- df_filtered$LogFC2
  x <- df_filtered$LogFC1
  
  ### Define axis title
  x_title <- paste0("Growth Arrest Effects (" ,combination_names[col]," Fe)")
  y_title <- paste0("Growth Arrest Effects (" ,combination_names[col+4]," Fe)")
  
  ### Pearson Correlation
  test=cor.test(y,x,method = "pearson")
  print(test)
  
  ### plot
  pl_scatter=ggplot(df_filtered)+
    geom_point(aes(x=LogFC1,y=LogFC2),color=color[col],pch=21,alpha=0.8)+
    # scale_fill_manual(values = color)+
    # geom_abline(slope=0,intercept=0,color="black")+
    xlab(x_title)+
    ylab(y_title)+
    theme(legend.position.inside = c(0.8,0.2),
          panel.background = element_rect(fill="grey90"),
          panel.grid.major = element_line(colour = "grey95"),
          panel.grid.minor = element_line(colour = "grey95")
    )
  
  pl_scatter_1 <- pl_scatter+
    geom_smooth(aes(x=LogFC1,y=LogFC2),color="black",method=lm,formula = "y ~ x",show.legend = FALSE)+
    # geom_rug(aes(x=LogFC1,y=LogFC2),color=color[i])+
    geom_vline(xintercept = 0)+
    geom_hline(yintercept = 0)+
    ggtitle(paste0("R=",round(test$estimate,digits=2)," Num_of_common genes=",dim(df_filtered)[1]))+
    theme_minimal()
  
  print(pl_scatter_1)
  fig_repo <- "4_Figures_paper"
  dir.create(file.path(output_dir,fig_repo,"5_Figure_S1C_Scatter_plot_correlogram_same_C_source"),recursive = T,showWarnings = F)
  
  pdf(file.path(output_dir,fig_repo,"5_Figure_S1C_Scatter_plot_correlogram_same_C_source",paste0(x_title,"_vs_",y_title,".pdf")),width=9,height=7)
  print(pl_scatter_1)
  dev.off()
}









