#############################
### 0. Load dependencies. ###
#############################

{
  library(tidyverse)
  library(ggplot2)
  library(dplyr)
  library(corrplot)
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
output_dir <- "analysis/Outputs"

### load Iron response to stat table

iron_response_to_stat <- read.table(file.path(input_dir,"txt/iron_response_to_stat.txt"))

th=0.05 # threshold of the fdr

### pair combinations to correlate
combination_names_with <- combn(colnames(iron_response_to_stat)[1:4],m=2)
combination_names_without <- combn(colnames(iron_response_to_stat)[5:8],m=2)

### select the samples with or without iron
vector_to_loop <- seq_along(combination_names_with[1:ncol(combination_names_with)])

# Show the vector to loop through
print(vector_to_loop)

### Define colors
color_with <- "grey"
color_without <- "white"

### Select data and colors
# color <- color_with
# combination_names <- combination_names_with

# Ask the user to choose between combination_names_without and combination_names_with
data_to_select <- readline("Choose between 'combination_names_without' or 'combination_names_with': ")

# Check if the chosen variable is valid
if (data_to_select == "combination_names_without") {
  combination_names <- combination_names_without
  color <- color_without
} else if (data_to_select == "combination_names_with") {
  combination_names <- combination_names_with
  color <- color_with
} else {
  stop("The entered option is not valid. Please choose 'combination_names_without' or 'combination_names_with'.")
}

### Create an empty cor matrix
cor.matrix <- matrix(NA, nrow = 6, ncol = 1)
pvalues.matrix <- matrix(NA, nrow = 6, ncol = 1)

names_cor_matrix <- c()

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
  
  ### Create dir
  df_to_corr_dir <- paste0(combination_names[1,col],"_vs_",combination_names[2,col],"_Fig_S1B.txt")
  dir.create(file.path(input_dir,"txt/Correlogram_data_Fig_S1B"),showWarnings = F)
  write.table(df_filtered,file.path(input_dir,"txt/Correlogram_data_Fig_S1B",df_to_corr_dir))
  
  ####################
  ###Scatter plots ###
  ####################
  ### Define axis
  y <- df_filtered$LogFC2
  x <- df_filtered$LogFC1
  ### Define axis title
  x_title <- combination_names[2,col]
  y_title <- combination_names[1,col]
  
  ### Pearson Correlation
  test=cor.test(y,x,method = "pearson")
  print(test)
  ### Save estimate and pvalue
  names_cor_matrix[col] <- paste0(combination_names[1,col],"_vs_",combination_names[2,col])
  cor.matrix[col,] <- test$estimate
  pvalues.matrix[col,] <- -log10(test$p.value)
  # next
  ### plot
  pl_scatter=ggplot(df_filtered)+
    geom_point(aes(x=LogFC1,y=LogFC2),color="black",fill=color,pch=21,alpha=0.8)+
    # scale_fill_manual(values = color)+
    # geom_abline(slope=0,intercept=0,color="black")+
    xlab(x_title)+
    ylab(y_title)+
    theme(legend.position=c(0.8,0.2),
          panel.background = element_rect(fill="grey90"),
          panel.grid.major = element_line(colour = "grey95"),
          panel.grid.minor = element_line(colour = "grey95")
         )
  
  pl_scatter_1 <- pl_scatter+
    geom_smooth(aes(x=LogFC1,y=LogFC2),color="black",method=lm,formula = "y ~ x",show.legend = FALSE)+
    # geom_rug(aes(x=LogFC1,y=LogFC2),color=color)+
    geom_vline(xintercept = 0)+
    geom_hline(yintercept = 0)+
    ggtitle(paste0("R=",round(test$estimate,digits=2)," Num_of_common genes=",dim(df_filtered)[1]))
    # theme_minimal()
  
  print(pl_scatter_1)
  fig_repo <- "4_Figures_paper"
  scatter_dir <- "4_Figure_S1B_Scatter_plot_correlogram"
  dir.create(file.path(output_dir,fig_repo,scatter_dir),recursive = T,showWarnings = F)
  
  pdf(file.path(output_dir,fig_repo,scatter_dir,paste0(x_title,"_vs_",y_title,".pdf")),width=9,height=7)
  print(pl_scatter_1)
  dev.off()
}










