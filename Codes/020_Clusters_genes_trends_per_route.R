######################
### 0.Dependencies ###
######################
{
  library(igraph)
  # library(xlsx)
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(circlize)
  # library(ggsci)
  library(ggdendro)
  library(dendextend)
  # library(dendsort)
  library(tidyverse)
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
input_dir <- "Inputs/002_Processed_data"
output_dir <- "Outputs"
routes_dir <- file.path(input_dir,"Routes_GO")
go_terms_dir <- file.path(input_dir, "GO_terms")

###################################
### 3. Load GO genes per route  ###
###################################

split_gene_string <- function(x) {
  genes <- trimws(unlist(strsplit(paste(x, collapse = ", "), ",")))
  unique(genes[nzchar(genes)])
}

route_from_cluster_terms <- function(cluster_file, term_patterns) {
  if (!file.exists(cluster_file)) {
    return(character(0))
  }

  go_clusters <- read.table(cluster_file, check.names = FALSE)
  term_hits <- Reduce(`|`, lapply(term_patterns, grepl, x = rownames(go_clusters), ignore.case = TRUE))
  if (!any(term_hits)) {
    return(character(0))
  }

  split_gene_string(unlist(go_clusters[term_hits, c("query_G_D", "query_G", "query_G_L")]))
}

route_from_final_heatmaps <- function(final_gene_sets_file, heatmap_name) {
  if (!file.exists(final_gene_sets_file)) {
    return(character(0))
  }

  gene_sets <- read.table(final_gene_sets_file, header = TRUE, sep = "\t")
  unique(gene_sets$gene[gene_sets$heatmap == heatmap_name])
}

load_route_gene_lists <- function(routes_dir, go_terms_dir, input_dir) {
  route_files <- list.files(routes_dir, pattern = "\\.txt$", full.names = TRUE)

  if (length(route_files) > 0) {
    route_lists <- lapply(route_files, function(x) unique(read.table(x, header = FALSE)[[1]]))
    names(route_lists) <- tools::file_path_sans_ext(basename(route_files))
    return(route_lists)
  }

  cluster_file <- file.path(go_terms_dir, "cluster_GO_OR_4_fdr_0.05_down.txt")
  final_gene_sets_file <- file.path(dirname(input_dir), "..", "Outputs", "Final_GO_gene_heatmaps", "final_7_heatmaps_gene_sets.txt")

  route_lists <- list(
    Glycolisis = unique(c(
      route_from_final_heatmaps(final_gene_sets_file, "Heatmap_4_ED_glycolysis"),
      route_from_cluster_terms(cluster_file, c("glycolytic process", "glucose metabolic process", "gluconeogenesis"))
    )),
    Methylcitrate = route_from_cluster_terms(
      cluster_file,
      c("methylcitrate", "propionate metabolic process")
    ),
    Oxidative_phosphorilation = unique(c(
      route_from_final_heatmaps(final_gene_sets_file, "Heatmap_6_ED_OXPHOS"),
      route_from_cluster_terms(cluster_file, c("oxidative phosphorylation", "respiratory electron transport chain"))
    )),
    Pentose_phosphate = route_from_cluster_terms(cluster_file, c("pentose phosphate")),
    TCA = unique(c(
      route_from_final_heatmaps(final_gene_sets_file, "Heatmap_5_ED_TCA"),
      route_from_cluster_terms(cluster_file, c("tricarboxylic acid cycle"))
    ))
  )

  route_lists <- lapply(route_lists, function(x) unique(x[nzchar(x)]))
  route_lists
}

### load LogFC
LogFC_df <- read.table(file.path(input_dir,"txt/LogFC_0.05.txt"))
LogFC_df <- LogFC_df[,17:19]

colnames(LogFC_df) <- c("G_D","G", "G_L")

route_gene_lists <- load_route_gene_lists(routes_dir, go_terms_dir, input_dir)
print(vapply(route_gene_lists, length, integer(1)))

# Create cluster by metabolic route
clusters = lapply(route_gene_lists, function(x) {
  vec <- intersect(x, rownames(LogFC_df))
  LogFC_df_filt <- LogFC_df[vec,, drop = FALSE]
  return(LogFC_df_filt)
})

glycolisis <- clusters[["Glycolisis"]]
methylcitrate <- clusters[["Methylcitrate"]]
oxidative_phosph <- clusters[["Oxidative_phosphorilation"]]
pentose_phosphate <- clusters[["Pentose_phosphate"]]
TCA <- clusters[["TCA"]]

has_route <- function(df) {
  !is.null(df) && nrow(df) > 1 && ncol(df) > 0
}

# Plot
ht_by_route <- function(df,row_title="",column_title="",filename){

  ### Apply hclust function
  hclust_matrix <- as.matrix(df)
  dist_matrix <- dist(hclust_matrix)
  hc_clust <- hclust(dist_matrix,method="ward.D2")
  col_dend <- hclust(dist(t(hclust_matrix)))
  columns_name <- c("G_D","G","G_LCFA")
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

  dir.create(dirname(filename),recursive = T,showWarnings=F)

  pdf(file =filename ,width=size[1]+0.5,height=size[2]+1)
  draw(ht_plt)
  dev.off()

  return(draw(ht_plt))
}

filename <- file.path(output_dir,"001_Figures_paper","15_Figure_4C_Heatmap_glycolisis.pdf")
if (has_route(glycolisis)) glycolisis_plt <- ht_by_route(glycolisis,row_title="",column_title="Interactions",filename)

# filename <- file.path(output_dir,"001_Figures_paper","15_Figure_4C_Heatmap_methylcitrate.pdf")
# methylcitrate_plt <- ht_by_route (methylcitrate,row_title="",column_title="Interactions",filename)

filename <- file.path(output_dir,"001_Figures_paper","15_Figure_4C_Heatmap_oxidative_phosph.pdf")
if (has_route(oxidative_phosph)) oxidative_phosph_plt <- ht_by_route(oxidative_phosph,row_title="",column_title="Interactions",filename)

# filename <- file.path(output_dir,"001_Figures_paper","15_Figure_4C_Heatmap_pentose_phosphate.pdf")
# pentose_phosphate_plt <- ht_by_route (pentose_phosphate,row_title="",column_title="Interactions",filename)

filename <- file.path(output_dir,"001_Figures_paper","15_Figure_4C_Heatmap_TCA.pdf")
if (has_route(TCA)) TCA_plt <- ht_by_route(TCA,row_title="",column_title="Interactions",filename)

#####################
### Bars per gene ###
#####################

ht_by_gene <- function(df,name_of_path){

 for(i in seq_along(rownames(df))){
   print(i)

  hclust_matrix <- as.matrix(df)
  columns_name <- c("G_D","G","G_LCFA")
  colnames(hclust_matrix) <- columns_name

  # max(abs(df))
  # color_breaks <- c(-(max(abs(df))-0.5) , 0, max(abs(df))-0.5 )
  color_breaks <- c(-3 , 0, 3 )

  # my_palette <- c( "yellow",
  #                  "blue",
  #                  colorRampPalette(rev(brewer.pal(8, "Spectral")))(n = length(color_breaks[2:4])-1),
  #                  "red")
  col_fun = colorRamp2(color_breaks, c("blue","white","red"))
  # col_fun = colorRamp2(color_breaks, my_palette)

  ht_plt <- Heatmap(t(hclust_matrix[i,]),
                    na_col = "grey2",
                    col = col_fun,
                    # split = k,
                    name="Log2(OR)",
                    column_order = columns_name,
                    show_column_names = T,
                    column_names_gp = gpar(fontsize = 12),
                    row_names_gp = gpar(fontsize = 6),
                    column_title = "",
                    column_title_side = "bottom",
                    row_dend_reorder=F,
                    column_dend_reorder = F,
                    border_gp = gpar(col = "black", lty = 2),
                    # heatmap_height = unit(6, "cm"),
                    # heatmap_width = unit(8, "cm"),
                    width = unit(0.4, "snpc"),
                    height = unit(0.1, "snpc"),
                    row_title = row_title,
                    row_title_gp = gpar(fontize = 2),
                    # row_dend_side = "right",
                    row_names_side = "left",
                    # row_dend_width = unit(2, "cm"),
                    show_row_names = T ,
                    show_row_dend = F,
                    heatmap_legend_param = list(title = "LogFC",
                                                title_position = "leftcenter-rot",
                                                labels_gp = gpar(font = 3),
                                                title_gp = gpar( fontsize = 8),
                                                at= c(-3, 0, 3)))


  draw(ht_plt)

  size <- calc_ht_size(ht_plt)

  dir.create(file.path(output_dir,"7_Bars_GO_map",name_of_path),recursive = T,showWarnings=F)
  filename <- file.path(output_dir,"7_Bars_GO_map",name_of_path,paste0(rownames(df[i,]),".pdf"))


  pdf(file =filename ,width=size[1]+0.5,height=size[2]+1)
  draw(ht_plt)
  dev.off()

 }
}

# Glycolisis
row_title <- ""
if (has_route(glycolisis)) ht_by_gene(df=glycolisis,name_of_path="Glycolisis")

# Methylcitrate
if (has_route(methylcitrate)) ht_by_gene(df=methylcitrate,name_of_path="Methylcitrate")

# Oxidative_phosphorilation
if (has_route(oxidative_phosph)) ht_by_gene(df=oxidative_phosph,name_of_path="Oxidative_phosphorilation")

# Pentose_phosphate
if (has_route(pentose_phosphate)) ht_by_gene(df=pentose_phosphate,name_of_path="Pentose_phosphate")

# TCA
if (has_route(TCA)) ht_by_gene(df=TCA,name_of_path="TCA")

#################################
### Trends by clusters Fig. ###
#################################

trend_function <- function(df,name_of_path){

# Compute the mean by cluster
result_mean <- df %>%
  summarize(
    Mean_G_D = mean(G_D),
    Mean_G = mean(G),
    Mean_G_LCFA = mean(G_L),
    sd_G_D = sd(G_D),
    sd_G = sd(G),
    sd_G_LCFA = sd(G_L)
  )

### inspect the result
print(result_mean)
write.table(result_mean,file = file.path(input_dir,paste0(name_of_path,"_mean_per_cluster.txt")))

result_sd <- result_mean[4:6]
result_mean <- result_mean[1:3]

data_mean <- pivot_longer(result_mean,cols = c(Mean_G_D,Mean_G,Mean_G_LCFA),names_to = "Conditions",values_to = "Mean")
data_mean$Conditions <- factor(data_mean$Conditions,levels = c("Mean_G_D", "Mean_G", "Mean_G_LCFA") )

data_sd <- pivot_longer(result_sd,cols = c(sd_G_D,sd_G,sd_G_LCFA),names_to = "Conditions",values_to = "sd")
data_sd$Conditions <- factor(data_sd$Conditions,levels = c("Mean_G_D", "Mean_G", "Mean_G_LCFA") )

data <- cbind(data_mean, data_sd[,2])

myarrow=arrow(angle = 15, ends = "last", type = "closed")
col <- c("green", "red", "blue")

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

filename <- file.path(output_dir,"001_Figures_paper","Trends",paste0(name_of_path,"_plot_mean.pdf"))
dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
plot_and_save(data, filename)
return(data)
}

# Glycolisis
route_means <- list()
if (has_route(glycolisis)) route_means[["Glycolysis"]] <- trend_function(df=glycolisis,name_of_path="Glycolisis")

# Methylcitrate
if (has_route(methylcitrate)) route_means[["Methylcitrate"]] <- trend_function(df=methylcitrate,name_of_path="Methylcitrate")

# Oxidative_phosphorilation
if (has_route(oxidative_phosph)) route_means[["Oxidative_phosphorylation"]] <- trend_function(df=oxidative_phosph,name_of_path="Oxidative_phosphorilation")

# Pentose_phosphate
if (has_route(pentose_phosphate)) route_means[["Pentose_phosphate"]] <- trend_function(df=pentose_phosphate,name_of_path="Pentose_phosphate")

# TCA
if (has_route(TCA)) route_means[["TCA"]] <- trend_function(df=TCA,name_of_path="TCA")

dataframe <- bind_rows(route_means, .id = "Routes")

myarrow=arrow(angle = 15, ends = "last", type = "closed")
route_colors <- c(
  "Glycolysis" = "red",
  "Methylcitrate" = "purple",
  "Oxidative_phosphorylation" = "green",
  "Pentose_phosphate" = "orange",
  "TCA" = "blue"
)
route_colors <- route_colors[intersect(names(route_colors), unique(dataframe$Routes))]
x <- dataframe$Conditions

# dataframe$Cond_route <- paste(dataframe$Conditions,dataframe$Routes,sep="_")
### Plot eerything in a same plot

mean_plt_fun <- ggplot(data = dataframe,aes(x = x,y = abs(Mean), color = Routes))+
  geom_line(aes(group = Routes),arrow=myarrow,lwd=2)+
  geom_point(aes(color = Routes),size=4)+
  scale_color_manual(name = "Routes",values = route_colors)+
  # geom_errorbar(aes(ymin = abs(Mean)-sd, ymax=abs(Mean)+sd),
  #               linewidth = 0.6, width = 0.3, position = position_dodge(0.2),alpha=0.5)+
  ylab("Mean")+
  ggtitle("")+
  theme_classic()+
  theme(
    # axis.title.x = element_text(size = 14),  # Tamaño del título del eje X
    # axis.title.y = element_text(size = 14),  # Tamaño del título del eje Y
    axis.text.x = element_text(size = 12),   # Tamaño de los valores del eje X
    axis.text.y = element_text(size = 12),   # Tamaño de los valores del eje Y
    # plot.title = element_text(size = 16)     # Tamaño del título del gráfico
  )

print(mean_plt_fun)

filename <- file.path(output_dir,"001_Figures_paper","16_Figure_4C_Trends_per_route_error_bar.pdf")
dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)

pdf(file =filename,width = 9)
print(mean_plt_fun)
dev.off()
# getwd()
