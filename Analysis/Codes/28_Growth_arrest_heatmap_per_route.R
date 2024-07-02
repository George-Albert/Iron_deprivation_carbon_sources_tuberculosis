######################
### 0.Dependencies ###
######################
{
  library(igraph)
  library(xlsx)
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(circlize)
  library(ggsci)
  library(ggdendro)
  library(dendextend)
  library(dendsort)
}

###########################
### 1.Declare Functions ###
###########################

dcols=function(x){data.frame(colnames(x))}
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
routes_dir <- file.path(input_dir,"Routes_GO")

###################################
### 3. Load GO genes per route  ###
###################################

list_files <- list.files(routes_dir,pattern = ".*.txt")
datalist = lapply(list_files, function(x) read.table(file.path(routes_dir,x),header=F))
names(datalist) <- gsub(".txt","", list_files)
### Fix the RV codes
datalist[[1]][2,] <- "Rv1099c"  
datalist[[3]][19,] <- "Rv3043c"

### load LogFC
LogFC_df <- read.table(file.path(input_dir,"txt/LogFC_0.05.txt"))
LogFC_df <- LogFC_df[,9:16]

colnames(LogFC_df) <- c("G_D","G", "G_L")
colnames(LogFC_df) <-gsub(".log2FoldChange_shrunken","",colnames(LogFC_df))

# Create cluster by metabolic route 
clusters = lapply(datalist, function(x) {
  
  vec <- unlist(x,use.names = F)
  LogFC_df_filt <- LogFC_df[vec,]
  return(LogFC_df_filt)
  
})

order_columns <- c("Phase_effects_with_Fe_G_D", "Phase_effects_without_Fe_G_D",
                   "Phase_effects_with_Fe_G", "Phase_effects_without_Fe_G",
                   "Phase_effects_with_Fe_G_L", "Phase_effects_without_Fe_G_L",
                   "Phase_effects_with_Fe_L", "Phase_effects_without_Fe_L")

glycolisis <- clusters[[1]]
glycolisis <- glycolisis[,order_columns]

methylcitrate <- clusters[[2]]
methylcitrate <- methylcitrate[,order_columns]

oxidative_phosph <- clusters[[3]]
oxidative_phosph <- oxidative_phosph[,order_columns]

pentose_phosphate <- clusters[[4]]
pentose_phosphate <- pentose_phosphate[,order_columns]

TCA <- clusters[[5]]
TCA <- TCA[,order_columns]

# # Create the base components
base_elements <- c("G-D", "G", "G-L", "L")
suffixes <- c("_with_Fe", "_without_Fe")
# Use expand.grid to create all combinations and then paste them together
combinations <- expand.grid(base_elements, suffixes)
combinations <- combinations[order(combinations$Var1), ]
# Combine the elements to create the desired names
vector_with_without_Fe <- paste0(combinations$Var1, combinations$Var2)
columns_name <- factor(vector_with_without_Fe,levels =vector_with_without_Fe)## Create a column names vector
 
# Plot
ht_by_route <- function(df,row_title="",column_title="",filename){
  
  ### Apply hclust function
  hclust_matrix <- as.matrix(df)
  dist_matrix <- dist(hclust_matrix)
  hc_clust <- hclust(dist_matrix,method="ward.D2")
  col_dend <- hclust(dist(t(hclust_matrix)))
  colnames(hclust_matrix) <- columns_name
  # column_title = ""
  # row_title = ""
  max(abs(df))      
  # color_breaks <- c(-(max(abs(df))-0.5) , 0, max(abs(df))-0.5 )
  color_breaks <- c(-3 , 0, 3 )
  
  # my_palette <- c( "yellow",
  #                  "blue",
  #                  colorRampPalette(rev(brewer.pal(8, "Spectral")))(n = length(color_breaks[2:4])-1),
  #                  "red")
  col_fun = colorRamp2(color_breaks, c("blue","white","red"))
  # col_fun = colorRamp2(color_breaks, my_palette)
  
  ht_plt <- Heatmap(hclust_matrix,
                    na_col = "grey2",
                    col = col_fun,
                    # split = k,
                    name="Log2(OR)",
                    column_order = columns_name,
                    show_column_names = T,
                    column_names_gp = gpar(fontsize = 6),
                    row_names_gp = gpar(fontsize = 6),
                    column_title = column_title,
                    column_title_side = "bottom",
                    row_dend_reorder=T,
                    border_gp = gpar(col = "black", lty = 2),
                    # heatmap_height = unit(6, "cm"),
                    # heatmap_width = unit(8, "cm"),
                    width = unit(0.2, "snpc"), 
                    height = unit(0.4, "snpc"),
                    # show_column_dend = T,
                    # column_dend_side = "top",
                    cluster_rows = color_branches(hc_clust),
                    # cluster_columns = color_branches(col_dend),
                    row_title = row_title,
                    row_title_gp = gpar(fontize = 2),
                    row_dend_side = "right",
                    row_names_side = "left",
                    # row_dend_width = unit(2, "cm"),
                    show_row_names = T ,
                    show_row_dend = T,
                    cell_fun = function(j, i, x, y, width, height, fill) {
                      grid.text(round(hclust_matrix[i, j],digits = 2), x, y, gp = gpar(fontsize = 3))},
                    heatmap_legend_param = list(title = "LogFC",
                                                title_position = "leftcenter-rot",
                                                labels_gp = gpar(font = 3), 
                                                title_gp = gpar( fontsize = 8),
                                                at= c(-3, 0, 3)))
              
  
  draw(ht_plt)
  
  size <- calc_ht_size(ht_plt)
  
  dir.create(file.path(output_dir,"Figures_paper"),recursive = T,showWarnings=F)
  
  pdf(file =filename ,width=size[1]+0.5,height=size[2]+1)
  draw(ht_plt)
  dev.off()
  
  return(draw(ht_plt))
}

heat_dir <- "Heatmap_growth_arrest"
filename <- file.path(output_dir,"4_Figures_paper",paste0(heat_dir,"_glycolisis.pdf"))
glycolisis_plt <- ht_by_route (glycolisis,row_title="",column_title="Growth arrest effects",filename)

# filename <- file.path(output_dir,"4_Figures_paper","15_Figure_4C_Heatmap_methylcitrate.pdf")
# methylcitrate_plt <- ht_by_route (methylcitrate,row_title="",column_title="Interactions",filename)

filename <- file.path(output_dir,"4_Figures_paper",paste0(heat_dir,"_oxidative_phosph.pdf"))
oxidative_phosph_plt <- ht_by_route (oxidative_phosph,row_title="",column_title="Interactions",filename)

# filename <- file.path(output_dir,"4_Figures_paper","15_Figure_4C_Heatmap_pentose_phosphate.pdf")
# pentose_phosphate_plt <- ht_by_route (pentose_phosphate,row_title="",column_title="Interactions",filename)

filename <- file.path(output_dir,"4_Figures_paper",paste0(heat_dir,"_TCA.pdf"))
TCA_plt <- ht_by_route (TCA,row_title="",column_title="Interactions",filename)

#################################
### Trends by clusters Fig. ###
#################################
# 
# trend_function <- function(df,name_of_path){
#   
# # Compute the mean by cluster
# result_mean <- df %>%
#   summarize(
#     Mean_G_D = mean(G_D),
#     Mean_G = mean(G),
#     Mean_G_LCFA = mean(G_L),
#     sd_G_D = sd(G_D),
#     sd_G = sd(G),
#     sd_G_LCFA = sd(G_L)
#   )
# 
# ### inspect the result
# print(result_mean)
# write.table(result_mean,file = file.path(input_dir,paste0(name_of_path,"_mean_per_cluster.txt")))
# 
# result_sd <- result_mean[4:6]
# result_mean <- result_mean[1:3]
# 
# data_mean <- pivot_longer(result_mean,cols = c(Mean_G_D,Mean_G,Mean_G_LCFA),names_to = "Conditions",values_to = "Mean")
# data_mean$Conditions <- factor(data_mean$Conditions,levels = c("Mean_G_D", "Mean_G", "Mean_G_LCFA") )
# 
# data_sd <- pivot_longer(result_sd,cols = c(sd_G_D,sd_G,sd_G_LCFA),names_to = "Conditions",values_to = "sd")
# data_sd$Conditions <- factor(data_sd$Conditions,levels = c("Mean_G_D", "Mean_G", "Mean_G_LCFA") )
# 
# data <- cbind(data_mean, data_sd[,2])
# 
# myarrow=arrow(angle = 15, ends = "last", type = "closed")
# col <- c("green", "red", "blue")
# 
# plot_and_save <- function(dataframe, filename) {
#   
#   mean_plt_fun <- ggplot(data=dataframe,aes(x=Conditions,y=Mean, colour=Conditions,group=1))+
#     geom_line(color="black",arrow=myarrow)+
#     geom_point(size=4)+
#     scale_color_manual(name="Condition",values = col)+
#     ylab("Mean")+
#     ggtitle("")+
#     theme_classic()
#   
#   print(mean_plt_fun)
#   
#   pdf(file =filename)
#   print(mean_plt_fun)
#   dev.off()
# }
# 
# filename <- file.path(output_dir,"4_Figures_paper","Trends_growth_arrest",paste0(name_of_path,"_plot_mean.pdf"))
# plot_and_save(data, filename)
# return(data)
# }
# 
# # Glycolisis
# Gly_mean <- trend_function(df=glycolisis,name_of_path="Glycolisis")
# 
# # Methylcitrate
# trend_function(df=methylcitrate,name_of_path="Methylcitrate")
# 
# # Oxidative_phosphorilation
# OP_mean <-trend_function(df=oxidative_phosph,name_of_path="Oxidative_phosphorilation")
# 
# # Pentose_phosphate
# trend_function(df=pentose_phosphate,name_of_path="Pentose_phosphate")
# 
# # TCA
# TCA_mean <-trend_function(df=TCA,name_of_path="TCA")
# 
# dataframe <- rbind(Gly_mean,OP_mean,TCA_mean)
# dataframe$Routes <- rep(c("Glycolysis","Oxidative_phosphorylation","TCA"),each=3)
# 
# myarrow=arrow(angle = 15, ends = "last", type = "closed")
# col <- c("red", "green", "blue")
# x <- dataframe$Conditions
# 
# # dataframe$Cond_route <- paste(dataframe$Conditions,dataframe$Routes,sep="_")
# ### Plot eerything in a same plot
# 
# mean_plt_fun <- ggplot(data = dataframe,aes(x = x,y = abs(Mean), color = Routes))+
#   geom_line(aes(group = Routes),arrow=myarrow,lwd=2)+
#   geom_point(aes(color = Routes),size=4)+
#   scale_color_manual(name = "Routes",values = col)+
#   geom_errorbar(aes(ymin = abs(Mean)-sd, ymax=abs(Mean)+sd),
#                 linewidth = 0.6, width = 0.1, position = position_dodge(0.2),alpha=0.5)+
#   ylab("Mean")+
#   ggtitle("")+
#   theme_classic()
# 
# print(mean_plt_fun)
# 
# filename <- file.path(output_dir,"4_Figures_paper","16_Figure_4C_Trends_per_route_error_bar.pdf")
# 
# pdf(file =filename,width = 9)
# print(mean_plt_fun)
# dev.off()
# # getwd()
# 






