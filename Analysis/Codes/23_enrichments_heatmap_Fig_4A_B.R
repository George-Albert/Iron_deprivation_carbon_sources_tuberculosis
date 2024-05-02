######################
### 0.Dependencies ###
######################
{
  library(qvalue)
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
input_folder <- "Analysis/Inputs/2_Processed_data/GO_terms"
output_dir <- "Analysis/Outputs"
enrichment_dir <- file.path(output_dir,"6_Enrichment_GO")
dir_of_gene_list <- "Analysis/Outputs/Enhanced_and_damped_genes"

#########################
### 3. Load GO lists  ###
#########################

lista_ED <- readRDS(file=file.path(input_folder,paste0("lista_down_GO_enrichment_OR_4_fdr_0.05.rds")))
lista_EU <- readRDS(file=file.path(input_folder,paste0("lista_up_GO_enrichment_OR_4_fdr_0.05.rds")))

OR_list_up <- lapply(lista_EU, function(x)
{
  # x[["OR"]][is.infinite(x[["OR"]])] <- 1e100
  x <- x[c("description","OR","Min_value_CI","Max_value_CI","fdr")]
  # x <- x[,1:2]
  # as.matrix(x)
  return(x)
}) 
OR_list_down <- lapply(lista_ED, function(x)
{
  # x[["OR"]][is.infinite(x[["OR"]])] <- 1e100
  x <- x[c("description","OR","Min_value_CI","Max_value_CI","fdr")]
  # x <- x[,1:2]
  # as.matrix(x)
  return(x)
}) 

df <- cbind(OR_list_up[[1]],OR_list_up[[2]][2:5],OR_list_up[[3]][2:5])
df_down <- cbind(OR_list_down [[1]],OR_list_down [[2]][2:5],OR_list_down [[3]][2:5])

colnames(df) <- c("description","OR_GD","lower_CI_GD","upper_CI_GD","FDR_GD","OR_G","lower_CI_G",
                  "upper_CI_G","FDR_G","OR_GL","lower_CI_GL","upper_CI_GL","FDR_GL")
colnames(df_down) <- colnames(df)
rownames(df)      <- df$description
rownames(df_down) <- df_down$description

### Save the max and min values distinct from Inf
max_values_df <- apply(df[,2:ncol(df)], 2, function(x) max(x[is.finite(x)], na.rm = TRUE))
min_values_df <- apply(df[,2:ncol(df)], 2, function(x) min(x[is.finite(x)], na.rm = TRUE))

max_values_df_down <- apply(df_down[,2:ncol(df_down)], 2, function(x) max(x[is.finite(x)], na.rm = TRUE))
min_values_df_down <- apply(df_down[,2:ncol(df_down)], 2, function(x) min(x[is.finite(x)], na.rm = TRUE))

for (col in colnames(df[,2:ncol(df)])) {
  
  replace_value <- max_values_df[col] + 2
  df[is.infinite(df[, col]), col] <- replace_value
  
  replace_value <- max_values_df_down[col] + 2
  df_down[is.infinite(df_down[, col]), col] <- replace_value
  
}

###########################
### Enhanced UP heatmap ###
###########################

### Apply hclust function
hclust_matrix <- as.matrix(df[,c(2,6,10)])
# hclust_matrix[is.infinite(hclust_matrix)] <- 1e100
hclust_matrix <- apply(hclust_matrix, MARGIN=c(1, 2), FUN = log2)
max_values_hclust <- apply(hclust_matrix, 2, function(x) max(x[is.finite(x)]))
hclust_matrix[(is.infinite(hclust_matrix))] <- -1
dist_matrix <- dist(hclust_matrix)
OR_clust <- hclust(dist_matrix,method="ward.D2")
col_dend <- hclust(dist(t(hclust_matrix)))
columns_name <- c("G_D","G","G_LCFA")
colnames(hclust_matrix) <- columns_name
columns_name <- colnames(hclust_matrix)
column_title = "Conditions"
clust_name <- "GO"

color_breaks <- c(-1, 0, 3)
# my_palette <- c( "yellow",
#                  "blue",
#                  colorRampPalette(rev(brewer.pal(8, "Spectral")))(n = length(color_breaks[2:4])-1),
#                  "red")
col_fun = colorRamp2(color_breaks, c("grey","white","red"))
# col_fun = colorRamp2(color_breaks, my_palette)
k=6
ht_plt <- Heatmap(hclust_matrix,
                  na_col = "grey2",
                  col = col_fun,
                  split = k,
                  name="Log2(OR)",
                  column_order = columns_name,
                  show_column_names = T,
                  column_names_gp = gpar(fontsize = 6),
                  row_names_gp = gpar(fontsize = 6),
                  column_title = column_title,
                  column_title_side = "bottom",
                  # row_dend_reorder=T,
                  border_gp = gpar(col = "black", lty = 2),
                  # heatmap_height = unit(6, "cm"),
                  # heatmap_width = unit(8, "cm"),
                  width=unit(2, "cm"),
                  # show_column_dend = T,
                  # column_dend_side = "top",
                  cluster_rows = color_branches(OR_clust,k=k),
                  # cluster_columns = color_branches(col_dend),
                  row_title = clust_name,
                  row_title_gp = gpar(fontize = 2),
                  row_dend_side = "right",
                  row_names_side = "left",
                  # row_dend_width = unit(2, "cm"),
                  show_row_names = T ,
                  show_row_dend = T,
                  cell_fun = function(j, i, x, y, width, height, fill) {
                    grid.text(round(hclust_matrix[i, j],digits = 2), x, y, gp = gpar(fontsize = 3))},
                  heatmap_legend_param = list(title = "Log2(OR)",
                                              title_position = "leftcenter-rot",
                                              labels_gp = gpar(font = 3),
                                              title_gp = gpar( fontsize = 8)))

draw(ht_plt)
ht_plt_up <- ht_plt

size <- calc_ht_size(ht_plt)
cluster <- dendextend:::cutree(OR_clust, k=k,order_clusters_as_data = F)
cluster_df <- data.frame(cluster)
idx <- match(rownames(hclust_matrix), names(cluster))
hcl_matrix <- cbind(hclust_matrix,cluster=cluster[idx])
hclust_df <- data.frame(hcl_matrix)
hclust_df <- hclust_df[order(hclust_df$cluster),]
hclust_df_up <- hclust_df

OR <- 4
fdr <- 0.05
write.table(hclust_df_up,file = file.path(input_folder,paste0("Cluster_OR_",OR,"_fdr_",fdr,"_up_k=",k,".txt")))

dir.create(file.path(output_dir,"4_Figures_paper"),recursive = T,showWarnings=F)
filename <- file.path(output_dir,"4_Figures_paper","13_Figure_4A_Heatmap_GO_up.pdf")

pdf(file =filename ,width=size[1]+0.5,height=size[2]+1)
draw(ht_plt)
dev.off()
  
#############################
### Enhanced DOWN heatmap ###
#############################

### Apply hclust function
hclust_matrix <- as.matrix(df_down[,c(2,6,10)])
# hclust_matrix[is.infinite(hclust_matrix)] <- 1e100
hclust_matrix <- apply(hclust_matrix, MARGIN=c(1, 2), FUN = log2)
max_values_hclust <- apply(hclust_matrix, 2, function(x) max(x[is.finite(x)]))
hclust_matrix[(is.infinite(hclust_matrix))] <- -1
dist_matrix <- dist(hclust_matrix)
OR_clust <- hclust(dist_matrix,method="ward.D2")
# OR_clust <- as.dendrogram(OR_clust)
col_dend <- hclust(dist(t(hclust_matrix)))
columns_name <- c("G_D","G","G_LCFA")
colnames(hclust_matrix) <- columns_name
columns_name <- colnames(hclust_matrix)
column_title = "Conditions"
clust_name <- "GO"

color_breaks <- c(-1, 0, 4)
# my_palette <- c( "yellow",
#                  "blue", 
#                  colorRampPalette(rev(brewer.pal(8, "Spectral")))(n = length(color_breaks[2:4])-1),
#                  "red")
col_fun = colorRamp2(color_breaks, c("grey","white","blue"))
# col_fun = colorRamp2(color_breaks, my_palette)
k=4
# plot(color_branches(OR_clust,k=k,groupLabels =T))

ht_plt1 <- Heatmap(hclust_matrix,
                  na_col = "grey2",
                  col = col_fun,
                  split = k,
                  name="Log2(OR)",
                  column_order = columns_name,
                  show_column_names = T,
                  column_names_gp = gpar(fontsize = 6),
                  row_names_gp = gpar(fontsize = 6),
                  column_title = column_title,
                  column_title_side = "bottom",
                  # row_dend_reorder=T,
                  border_gp = gpar(col = "black", lty = 2),
                  # heatmap_height = unit(6, "cm"),
                  # heatmap_width = unit(8, "cm"),
                  width=unit(2, "cm"),
                  # show_column_dend = T,
                  # column_dend_side = "top",
                  cluster_rows = color_branches(OR_clust,k=k,groupLabels =T),
                  # cluster_row_slices=T,
                  # cluster_columns = color_branches(col_dend),
                  row_title = clust_name,
                  # row_title_rot = 0,
                  row_title_gp = gpar(fontize = 2),
                  row_dend_side = "right",
                  row_names_side = "left",
                  # row_dend_reorder = F,
                  # row_dend_width = unit(2, "cm"),
                  show_row_names = T ,
                  show_row_dend = T,
                  cell_fun = function(j, i, x, y, width, height, fill) {
                    grid.text(round(hclust_matrix[i, j],digits = 2), x, y, gp = gpar(fontsize = 3))},
                  heatmap_legend_param = list(title = "Log2(OR)",
                                              title_position = "leftcenter-rot",
                                              labels_gp = gpar(font = 3),
                                              title_gp = gpar( fontsize = 8)))


draw(ht_plt1)
ht_plt_down <- ht_plt1

size <- calc_ht_size(ht_plt1)
cluster <- dendextend:::cutree(OR_clust, k=k,order_clusters_as_data = F)
cluster_df <- data.frame(cluster)
idx <- match(rownames(hclust_matrix), names(cluster))
hcl_matrix <- cbind(hclust_matrix,cluster=cluster[idx])
hclust_df <- data.frame(hcl_matrix)
hclust_df <- hclust_df[order(hclust_df$cluster),]

hclust_df_down <- hclust_df

write.table(hclust_df_down,file = file.path(input_folder,paste0("Cluster_OR_",OR,"_fdr_",fdr,"_down_k=",k,".txt")))

filename <- file.path(output_dir,"4_Figures_paper","13_Figure_4A_Heatmap_GO_down.pdf")

pdf(file =filename ,width=size[1]+0.5,height=size[2]+1)
draw(ht_plt_down)
dev.off()

  # ha = HeatmapAnnotation(pt = anno_points(1:61))
  # ha = rowAnnotation(pt = anno_points(1:61))
  
# ht_plt + rowAnnotation(foo = anno_block(
#   panel_fun = function(index, levels) {
#     grid.rect(gp = gpar(col = "black"))
#     ggplot(aes)
#   },
#   width = unit(2, "cm")))

#################################
### Trends by clusters Fig.3G ###
#################################

hclust_df_down <- read.table(file = file.path(input_folder,"Cluster_OR_4_fdr_0.05_down_k=4.txt"))

hclust_df_up <- read.table(file = file.path(input_folder,"Cluster_OR_4_fdr_0.05_up_k=6.txt"))

# Compute the mean by cluster
result_down <- hclust_df_down %>%
  group_by(cluster) %>%  # Group by column "cluster"
  summarize(
    Mean_G_D = mean(G_D),
    Mean_G = mean(G),
    Mean_G_LCFA = mean(G_LCFA)
  )

### inspect the result
print(result_down)

# Compute the mean by cluster
result_up <- hclust_df_up %>%
  group_by(cluster) %>%  # Group by column "cluster"
  summarize(
    Mean_G_D = mean(G_D),
    Mean_G = mean(G),
    Mean_G_LCFA = mean(G_LCFA)
  )

### inspect the result
print(result_up)

write.table(result_down,file = file.path(input_folder,"table_down_mean_per_cluster.txt"))
write.table(result_up,file = file.path(input_folder,"table_up_mean_per_cluster.txt"))


data_up <- pivot_longer(result_up,cols = !cluster,names_to = "Conditions",values_to = "Mean")
data_up$Conditions <- factor(data_up$Conditions,levels = c("Mean_G_D", "Mean_G", "Mean_G_LCFA") )

data_down <- pivot_longer(result_down,cols = !cluster,names_to = "Conditions",values_to = "Mean")
data_down$Conditions <- factor(data_down$Conditions,levels = c("Mean_G_D", "Mean_G", "Mean_G_LCFA") )

myarrow=arrow(angle = 15, ends = "last", type = "closed")
col <- c("#852121", "#8f4a17", "#2c682c")

list_cluster_down <- split(data_down,data_down$cluster)
list_cluster_up <- split(data_up,data_up$cluster)


plot_and_save <- function(dataframe, filename) {
  
  mean_plt_fun <- ggplot(data=dataframe,aes(x=Conditions,y=Mean, colour=Conditions,group=1))+
    geom_line(color="black",arrow=myarrow)+
    geom_point(size=4)+
    scale_color_manual(name="Condition",values = col)+
    ylab("Mean")+
    ggtitle("")+
    theme_classic()
    
  print(mean_plt_fun)
    
  pdf(file =filename)
  print(mean_plt_fun)
  dev.off()
}
dir.create(file.path(output_dir,"4_Figures_paper/Trends"),recursive = T,showWarnings = F)
for (i in 1:length(list_cluster_down)) {
  filename <- file.path(output_dir,"4_Figures_paper/Trends",paste("plot_mean_k_", i,"_down.pdf"))
  plot_and_save(list_cluster_down[[i]], filename)
}

for (i in 1:length(list_cluster_up)) {
  filename <- file.path(output_dir,"4_Figures_paper/Trends",paste("plot_mean_k_", i,"_up.pdf"))
  plot_and_save(list_cluster_up[[i]], filename)
}


####################################################################
### fractions of most enriched & least enriched C-source  Fig 4B ###
####################################################################


up <- df[,c(2,6,10)]
down <- df_down[,c(2,6,10)]
Conditions <- factor(c("G_D","G","G_L"),levels = c("G_D","G","G_L"))

### UP ###
{
### Find the max by row
idx_max_x_fila <- apply(up, 1, function(row) which.max(row))

a = length(which(idx_max_x_fila==1))
b = length(which(idx_max_x_fila==2))
c = length(which(idx_max_x_fila==3))
total <- sum(a,b,c)

most_enrich <- c(a,b,c)
most_enrich_frac <- c(a/total*100,b/total*100,c/total*100)

### Find the min by row
idx_min_x_fila <- apply(up, 1, function(row) which.min(row))

d = length(which(idx_min_x_fila==1))
e = length(which(idx_min_x_fila==2))
f = length(which(idx_min_x_fila==3))
total <- sum(d,e,f)

least_enrich <- c(d,e,f)
least_enrich_frac <- c(d/total*100,e/total*100,f/total*100)
rest <- total - (most_enrich + least_enrich)   
rest_frac <- 100 - (most_enrich_frac + least_enrich_frac)   

up_tab <- data.frame(Conditions, Most_Enriched = most_enrich_frac, Least_Enriched = least_enrich_frac,
                       REST = rest_frac)
up_tab_pivot <- pivot_longer(up_tab,cols = !Conditions, names_to = "Categorie",values_to = "values")
up_tab_pivot$Categorie <- factor(up_tab_pivot$Categorie,levels = c("Most_Enriched","Least_Enriched","REST"))

up_tab_pivot_sub <- subset(up_tab_pivot, Categorie != "REST")
### DOWN ###

### Find the max by row
idx_max_x_fila <- apply(down, 1, function(row) which.max(row))

a = length(which(idx_max_x_fila==1))
b = length(which(idx_max_x_fila==2))
c = length(which(idx_max_x_fila==3))
total <- sum(a,b,c)

most_enrich <- c(a,b,c)
most_enrich_frac <- c(a/total*100,b/total*100,c/total*100)

### Find the min by row
idx_min_x_fila <- apply(down, 1, function(row) which.min(row))

d = length(which(idx_min_x_fila==1))
e = length(which(idx_min_x_fila==2))
f = length(which(idx_min_x_fila==3))
total <- sum(d,e,f)

least_enrich <- c(d,e,f)
least_enrich_frac <- c(d/total*100,e/total*100,f/total*100)
rest <- total - (most_enrich + least_enrich)   
rest_frac <- 100 - (most_enrich_frac + least_enrich_frac)   

down_tab <- data.frame(Conditions, Most_Enriched = most_enrich_frac, Least_Enriched = least_enrich_frac,
                       REST = rest_frac)
down_tab_pivot <- pivot_longer(down_tab,cols = !Conditions, names_to = "Categorie",values_to = "values")
down_tab_pivot$Categorie <- factor(down_tab_pivot$Categorie,levels = c("Most_Enriched","Least_Enriched","REST"))
down_tab_pivot_sub <- subset(down_tab_pivot, Categorie != "REST")

}

#############################
### Bar plot of fractions ###
#############################

color <- c("#852121", "#8f4a17", "#2c682c")

### UP bar fraction plot ###

bar_plt <- ggplot(up_tab_pivot_sub, aes(fill=Conditions, y=values, x=Categorie)) + 
  geom_bar(position=position_fill(reverse = TRUE), stat="identity",width = 1/2, color = "black")+
  coord_flip()+
  scale_y_continuous(labels = scales::percent_format())+
  scale_fill_manual(values = color)+
  theme_classic()+
  geom_hline(yintercept=0)+
  # theme(aspect.ratio = .9)+
  labs(x="Interactions", y = "Fractions")+
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
    data = up_tab_pivot_sub,
    aes(x=Categorie, y=values, label = sprintf("%0.1f", round(values, digits = 1))),
    hjust = 1,vjust=-4,
    position = position_fill(reverse = TRUE),
    check_overlap = T,
    # nudge_x = -0.3,
    colour = "black",
    bg.colour = "white",
    bg.r = 0.2,
    # family = "Econ Sans Cnd",
    size = 3)

bar_pl_1

pdf(file = file.path(output_dir,"4_Figures_paper","14_Figure_4B_Most_Least_enrich_UP_wo_REST.pdf"),width=7,height=4)
print(bar_pl_1)
dev.off()

### DOWN fraction plot ###

bar_plt2 <- ggplot(down_tab_pivot_sub, aes(fill=Conditions, y=values, x=Categorie)) + 
  geom_bar(position=position_fill(reverse = TRUE), stat="identity",width = 1/2, color = "black")+
  coord_flip()+
  scale_y_continuous(labels = scales::percent_format())+
  scale_fill_manual(values = color)+
  theme_classic()+
  geom_hline(yintercept=0)+
  # theme(aspect.ratio = .9)+
  labs(x="Interactions", y = "Fractions")+
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

bar_plt2


bar_pl_3 <- bar_plt2+
  geom_shadowtext(
    data = down_tab_pivot_sub,
    aes(x=Categorie, y=values, label = sprintf("%0.1f", round(values, digits = 1))),
    hjust = 2,vjust=-4,
    position = position_fill(reverse = TRUE),
    # nudge_x = -0.3,
    colour = "black",
    bg.colour = "white",
    bg.r = 0.2,
    # family = "Econ Sans Cnd",
    size = 3)
  
bar_pl_3
  
pdf(file = file.path(output_dir,"4_Figures_paper","14_Figure_4B_Most_Least_enrich_DOWN_wo_Rest.pdf"),width=7,height=4)
print(bar_pl_3)
dev.off()
