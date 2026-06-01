#############################
### 0. Load dependencies. ###
#############################
{
  library(tidyverse)
  library(readxl)
  library(writexl)
  library(ComplexHeatmap)
  library(circlize)
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
input_dir  <- "Inputs/002_Processed_data"
output_dir <- "Outputs"

###########################
####### 3. Load data ######
###########################

### Matrix of contrasts
contrast <- read.table(file.path(input_dir,"txt/contrast_matrix.txt"))

### feature data
feat_path = file.path(input_dir,"txt/feature_data_filtered.txt")
feature_data = read.table(feat_path)

### metadata
meta_path <- file.path(input_dir,"txt/metadata_32_samples.txt")
metadata = read.table(meta_path)

### non-coding reads
reads_nc <- read.table(file.path(input_dir,"txt/reads_nc.txt"),check.names = F)

### Contrast with the nomenclature used
df_name_contrast <- read.table(file.path(input_dir,"txt/contrasts_nomenclature.txt"))

### load the list
MyData <- readRDS(file = file.path(input_dir,"RDS/Contrasts_stat_0.05.RDS"))
# mydf <- readRDS("Contrasts_stat.RDS")

### load the logFc,BH and lfcSE dfs
LogFC_df <- read.table(file.path(input_dir,"txt/LogFC_0.05.txt"))
BH_df    <- read.table(file.path(input_dir,"txt/BH_0.05.txt"))
lfcSE_df <- read.table(file.path(input_dir,"txt/lfcSE_0.05.txt"))

### Extract the 8 Iron effects @  EXP and STAT
# Iron_effect_in_EXP_G_D.log2FoldChange_shrunken
# Iron_effect_in_EXP_G.log2FoldChange_shrunken
# Iron_effect_in_EXP_G_L.log2FoldChange_shrunken
# Iron_effect_in_EXP_L.log2FoldChange_shrunken
# Iron_effect_in_STAT_G_D.log2FoldChange_shrunken
# Iron_effect_in_STAT_G.log2FoldChange_shrunken
# Iron_effect_in_STAT_G_L.log2FoldChange_shrunken
# Iron_effect_in_STAT_L.log2FoldChange_shrunken

vec_iron <- c(5:8)
LogFC_df_iron <- LogFC_df[,vec_iron]
BH_df_iron    <- BH_df[,vec_iron]
lfcSE_df_iron <- lfcSE_df[,vec_iron]

### Change the columns names
colnames(LogFC_df_iron) <- paste0(substring(colnames(LogFC_df_iron),1,nchar(colnames(LogFC_df_iron))-24),".LogFC")
# colnames(BH_df_iron) <- paste0(colnames(BH_df_iron),".BH")
colnames(lfcSE_df_iron) <- paste0(substring(colnames(lfcSE_df_iron),1,nchar(colnames(lfcSE_df_iron))-9))

length(which(rownames(LogFC_df_iron)!=rownames(BH_df_iron)))
length(which(rownames(LogFC_df_iron)!=rownames(lfcSE_df_iron)))

###Create Iron effects
iron_effects <- data.frame(LogFC_df_iron,BH_df_iron)

### Save the Iron effects dataframe
write.table(iron_effects,file.path(input_dir,"txt/Iron_effects.txt"))
iron_effects_export <- iron_effects %>%
  rownames_to_column(var = "Gene")
write_xlsx(iron_effects_export,file.path(input_dir,"xlsx/Iron_effects.xlsx"))

rv_non_coding <- grep("RVnc", rownames(iron_effects))
rv_non_coding_custom <- rownames(reads_nc)

rv_nc_genes <- c(rownames(iron_effects)[rv_non_coding],rv_non_coding_custom)
iron_effects <- iron_effects[which(rownames(iron_effects) %in% rv_nc_genes),]
rv_non_coding <- data.frame(iron_effects)

### Save the Iron effects dataframe with non-coding genes
write.table(rv_non_coding,file.path(input_dir,"txt/Iron_effects_RVnc.txt"))
write_xlsx(rv_non_coding,file.path(input_dir,"xlsx/Iron_effects_RVnc.xlsx"))

{
  rv_non_coding$Pattern_G_D_stat   = "Background"
  rv_non_coding$Pattern_G_stat     = "Background"
  rv_non_coding$Pattern_G_L_stat   = "Background"
  rv_non_coding$Pattern_L_stat     = "Background"
  threshold <- 0.05
}

col_logFC <- c(1:4)
for (i in col_logFC){

  #Upregulated per condition
  rv_non_coding[which((rv_non_coding[(i+4)] < threshold) &
                       (rv_non_coding[i] > 0)),(i+8)] <- "UP"
  #Downregulated per condition
  rv_non_coding[which(rv_non_coding[(i+4)] < threshold &
                       rv_non_coding[i] < 0),(i+8)] <- "DOWN"

}

DE_genes_vec <- c()

for (i in c(9:12)){

  print(paste("There are", length(which(rv_non_coding[,i]=="UP" |
                                          rv_non_coding[,i]=="DOWN")), "genes in",colnames(rv_non_coding)[i]))

  DE_genes_vec[i-8] <- length(which(rv_non_coding[,i]=="UP" | rv_non_coding[,i]=="DOWN"))
}

names_samples <- gsub('Pattern_','',colnames(rv_non_coding[9:12]))
# col <- c("lightcoral","red2","lightgoldenrod","yellow3","skyblue","blue4",
#          "lawngreen","darkgreen")
col <- c( "#852121", "#8f4a17","#2c682c","#355f91")


DE_genes_df <- data.frame(DE=DE_genes_vec,Contrast=names_samples,color=col)

### Save the list of DE genes in Lipids
write.table(rv_non_coding,file.path(input_dir,"txt/DE_genes_at_stat_RVnc_th<0.05.txt"))

rv_non_coding$RV <- rownames(rv_non_coding)
rv_non_coding <- rv_non_coding %>% relocate(RV, .before = Iron_effect_in_STAT_G_D.LogFC )
write_xlsx(rv_non_coding,file.path(input_dir,"xlsx/DE_genes_at_stat_RVnc_th<0.05.xlsx"))

df <- DE_genes_df
df$Contrast <- factor(df$Contrast,levels = df$Contrast)
breaks_values <- c(pretty(df$DE))

bar_plt <- ggplot(df, aes(y=DE, x=Contrast, fill= Contrast)) +
  geom_bar(stat="identity",width = 0.8)+
  scale_fill_manual(values = col)+
  scale_y_continuous(expand = c(0,0), limits = c(0,35))+
  # scale_color_manual(values = col)+
  # geom_col_pattern(data=df[5:8,],pattern="stripe",
  #                  pattern_fill = color_vec[5:8],
  #                  pattern_angle = 45,
  #                  pattern_spacing = 0.05)+
  # scale_pattern_manual(values = color_vec[5:8])+
  # coord_flip()+
  theme_classic()+
  geom_hline(yintercept=0)+
  # theme(aspect.ratio = .9)+
  labs(x="Contrasts", y = "Number of Genes")+
  theme(
    axis.text.x   = element_text(size=12),
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
bar_plt

bar_pl_1 <- bar_plt+
  geom_shadowtext(
    data = df,
    aes(x=Contrast, y=DE, label = DE),
    hjust = 0.5,vjust = -0.5,
    position = "stack",
    # nudge_x = -0.3,
    colour = "black",
    bg.colour = "white",
    bg.r = 0.2,
    # family = "Econ Sans Cnd",
    size = 5)

bar_pl_1

pdf(file=file.path(output_dir,"001_Figures_paper","barplot_DE_noncod_genes_stat_effects.pdf"),width=7,height=7)
print(bar_pl_1)
dev.off()

############################
# Create Heatmap iron effects
############################

# Create a matrix for the heatmap

iron_effects_matrix <- as.matrix(iron_effects[,c(1:4)])

rownames(iron_effects_matrix) <- iron_effects$RV
# Define the color palette

color_palette <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
# Create the heatmap

columns_name <- colnames(iron_effects_matrix)
column_title <- "Iron Effects on Non-coding Genes"
clust_name <- "RV Non-coding Genes"

ht_plt <- Heatmap(iron_effects_matrix,
                  na_col = "grey2",
                  col = color_palette,
                  # split = k,
                  name="LogFC",
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
                  cluster_rows = TRUE,
                  # cluster_columns = color_branches(col_dend),
                  row_title = clust_name,
                  row_title_gp = gpar(fontize = 2),
                  row_dend_side = "right",
                  row_names_side = "left",
                  # row_dend_width = unit(2, "cm"),
                  show_row_names = T ,
                  show_row_dend = T,
                  cell_fun = function(j, i, x, y, width, height, fill) {
                    grid.text(round(iron_effects[i, j],digits = 2), x, y, gp = gpar(fontsize = 3))},
                  heatmap_legend_param = list(title = "LogFC",
                                              title_position = "leftcenter-rot",
                                              labels_gp = gpar(font = 3),
                                              title_gp = gpar( fontsize = 8)))

draw(ht_plt)
# Save the heatmap to a PDF file
pdf(file=file.path(output_dir,"001_Figures_paper","Heatmap_iron_effects_noncoding_genes.pdf"),width=7,height=7)
draw(ht_plt)
dev.off()



####################################
##  Heatmap version 2 wo numbers  ##
####################################

# Paletas
col_fun_up   <- colorRamp2(c(0, 2), c("linen", "firebrick4"))
col_fun_down <- colorRamp2(c(-2,0), c("#1a528b","linen"))

legend_up   <- Legend(title = "Upregulated", col_fun = col_fun_up, direction = "vertical")
legend_down <- Legend(title = "Downregulated", col_fun = col_fun_down, direction = "vertical")
combined_legend <- packLegend(legend_up, legend_down, direction = "vertical")

# Matrices
logfc_mat <- as.matrix(rv_non_coding[, grep("LogFC", colnames(rv_non_coding))])
colnames(logfc_mat) <- gsub("Iron_effect_in_(STAT_.*?)\\.LogFC", "\\1", colnames(logfc_mat))

pattern_mat <- as.matrix(rv_non_coding[, grep("^Pattern_", names(rv_non_coding))])
patterns <- as.data.frame(pattern_mat)
pattern_mat <- apply(pattern_mat, c(1,2), function(x) gsub('"', '', x))

# Quitar las comillas dobles de los valores
pattern_mat <- apply(pattern_mat, c(1,2), function(x) gsub('"', '', x))

# Layer function
layer_fun <- function(j, i, x, y, w, h, fill) {

  mapply(function(i_, j_, x_, y_, w_, h_) {
    value  <- logfc_mat[i_, j_]
    status <- pattern_mat[i_, j_]

    # Color de fondo según patrón
    color <- if (is.na(status) || status == "Background") {
      "linen"
    } else if (status == "UP") {
      col_fun_up(abs(value))
    } else if (status == "DOWN") {
      col_fun_down((value))
    } else {
      "white"
    }

    # Dibujar celda coloreada
    grid.rect(x_, y_, w_, h_, gp = gpar(fill = color, col = "black", lwd = 0.2))

    # Añadir texto solo si es DE (UP o DOWN)
    if (!is.na(status) && status %in% c("UP", "DOWN")) {
      grid.text(sprintf("%.2f", value), x_, y_, gp = gpar(fontsize = 6))
    }

  }, i, j, x, y, w, h)
}

pt_final <- Heatmap(logfc_mat,
                    rect_gp = gpar(col = "white", lwd = 0.5),
                    name = 'logFC',
                    cluster_rows = FALSE,
                    cluster_columns = FALSE,
                    column_names_gp = gpar(fontsize = 10, fontface = "italic"),
                    row_names_gp = gpar(fontsize = 8),
                    column_names_rot = 90,
                    show_row_dend = FALSE,
                    row_names_side = "left",
                    show_column_dend = FALSE,
                    show_column_names = TRUE,
                    show_heatmap_legend = FALSE,
                    width = unit(2, "cm"),
                    height = unit(20, "cm"),
                    row_gap = unit(1, "cm"),
                    layer_fun = layer_fun)

draw(pt_final, annotation_legend_list = combined_legend)

pdf(file=file.path(output_dir,"001_Figures_paper","Heatmap_noncoding_genes.pdf"),width=6,height=10)
draw(pt_final, annotation_legend_list = combined_legend)
dev.off()

######################################################################
##  Bar plot of DE genes in non-coding genes UP and DOWN regulated  ##
######################################################################

# 1. Convert to long format
pattern_long <- patterns %>%
  tibble::rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "Condition", values_to = "Status")

# 3. Contar genes up y down
counts <- pattern_long %>%
  filter(Status %in% c("UP", "DOWN")) %>%
  group_by(Condition, Status) %>%
  summarise(N = n(), .groups = "drop") %>%
  # Convertir "DOWN" a negativo
  mutate(N = ifelse(Status == "DOWN", -N, N))
# 4. Create the bar plot
# Etiquetas más limpias
counts$Condition <- gsub("Pattern_", "", counts$Condition)
counts$Condition <- factor(counts$Condition,
                             levels = c("G_D_stat", "G_stat", "G_L_stat", "L_stat"))
# Plot
bar_plt_1 <- ggplot(counts, aes(x = N, y = Condition, fill = Status)) +
  geom_col(width = 0.6) +
  scale_fill_manual(values = c("UP" = "firebrick4", "DOWN" = "#1a528b")) +
  geom_vline(xintercept = 0, color = "black") +
  labs(x = "Number of Genes", y = "Condition", title = "Up- and Down-regulated Genes by Condition") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")

bar_plt_1 <- ggplot(counts, aes(x = N, y = Condition, fill = Status)) +
  geom_col(width = 0.6, position = position_stack(reverse = TRUE)) +
  geom_text(aes(label = abs(N)),
            position = position_stack(vjust = 0.5, reverse = TRUE),
            size = 5, color = "white") +
  scale_fill_manual(
    values = c(
      "UP" = "firebrick4",
      "DOWN" = "#1a528b"
    )
  ) +
  geom_vline(xintercept = 0, color = "black") +
  labs(x = "Number of Genes", y = "Condition", title = "Up- and Down-regulated Genes by Condition") +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    panel.background = element_rect(fill = "linen", color = NA),
    plot.background = element_rect(fill = "linen", color = NA)
  )
bar_plt_1
# Save the bar plot

pdf(file=file.path(output_dir,"001_Figures_paper","Barplot_DE_genes_noncoding.pdf"),width=7,height=5)
print(bar_plt_1)
dev.off()


#####################################
##  Heatmap with just DE nc genes  ##
#####################################

# Create a matrix for the heatmap with just DE nc genes

# Identificar genes con al menos un UP o DOWN en pattern_mat
rows_de <- apply(pattern_mat, 1, function(row) any(row %in% c("UP", "DOWN")))

# Filtrar las matrices
logfc_mat_filtered    <- logfc_mat[rows_de, ]
pattern_mat_filtered  <- pattern_mat[rows_de, ]

layer_fun1 <- function(j, i, x, y, w, h, fill) {

  mapply(function(i_, j_, x_, y_, w_, h_) {
    value  <- logfc_mat_filtered[i_, j_]
    status <- pattern_mat_filtered[i_, j_]

    # Color de fondo según patrón
    color <- if (is.na(status) || status == "Background") {
      "linen"
    } else if (status == "UP") {
      col_fun_up(value)
    } else if (status == "DOWN") {
      col_fun_down((value))
    } else {
      "white"
    }

    # Dibujar celda coloreada
    grid.rect(x_, y_, w_, h_, gp = gpar(fill = color, col = "black", lwd = 0.2))

    # Añadir texto solo si es DE (UP o DOWN)
    if (!is.na(status) && status %in% c("UP", "DOWN")) {
      grid.text(sprintf("%.2f", value), x_, y_, gp = gpar(fontsize = 6))
    }

  }, i, j, x, y, w, h)
}

ht_plt3 <- Heatmap(logfc_mat_filtered,
                    rect_gp = gpar(col = "white", lwd = 0.5),
                    name = 'logFC',
                    cluster_rows = T,
                    cluster_columns = FALSE,
                    column_names_gp = gpar(fontsize = 10, fontface = "italic"),
                    row_names_gp = gpar(fontsize = 8),
                    column_names_rot = 90,
                    show_row_dend = FALSE,
                    row_names_side = "left",
                    show_column_dend = FALSE,
                    show_column_names = TRUE,
                    show_heatmap_legend = FALSE,
                    width = unit(2, "cm"),
                    height = unit(20, "cm"),
                    row_gap = unit(1, "cm"),
                    layer_fun = layer_fun1)

draw(ht_plt3, annotation_legend_list = combined_legend)

pdf(file=file.path(output_dir,"001_Figures_paper","Heatmap_just_DE_noncoding_genes.pdf"),width=6,height=10)
draw(ht_plt3, annotation_legend_list = combined_legend)
dev.off()

###########################################################################################


##########################################
##  Same analysis but for interactions  ##
##########################################


vec_stat <- c(17:20)
LogFC_df_stat <- LogFC_df[,vec_stat]
BH_df_stat <- BH_df[,vec_stat]
lfcSE_df_stat <- lfcSE_df[,vec_stat]

### Change the columns names
colnames(LogFC_df_stat) <- paste0(c("Interaction_G_D.LogFC",
                                    "Interaction_G.LogFC",
                                    "Interaction_G_L.LogFC",
                                    "Interaction_L.LogFC"))
colnames(BH_df_stat)    <- paste0(c("Interaction_G_D.BH",
                                    "Interaction_G.BH",
                                    "Interaction_G_L.BH",
                                    "Interaction_L.BH"))
colnames(lfcSE_df_stat) <- paste0(c("Interaction_G_D.lfcSE",
                                    "Interaction_G.lfcSE",
                                    "Interaction_G_L.lfcSE",
                                    "Interaction_L.lfcSE"))

length(which(rownames(LogFC_df_stat)!=rownames(BH_df_stat)))
length(which(rownames(LogFC_df_stat)!=rownames(lfcSE_df_stat)))

###Create int effects
int_effects <- data.frame(LogFC_df_stat,BH_df_stat)

### Save the int effects dataframe
write.table(int_effects,file.path(input_dir,"txt/Interaction_effects.txt"))
int_effects_export <- int_effects %>%
  rownames_to_column(var = "Gene")
write_xlsx(int_effects_export,file.path(input_dir,"xlsx/Interaction_effects.xlsx"))

rv_non_coding <- grep("RVnc", rownames(int_effects))
rv_non_coding_custom <- rownames(reads_nc)

rv_nc_genes <- c(rownames(int_effects)[rv_non_coding],rv_non_coding_custom)
int_effects <- int_effects[which(rownames(int_effects) %in% rv_nc_genes),]
rv_non_coding <- data.frame(int_effects)

### Save the Iron effects dataframe with non-coding genes
write.table(rv_non_coding,file.path(input_dir,"txt/Interaction_effects_RVnc.txt"))
write_xlsx(rv_non_coding,file.path(input_dir,"xlsx/Interaction_effects_RVnc.xlsx"))

{
  rv_non_coding$Pattern_G_D_int   = "Background"
  rv_non_coding$Pattern_G_int    = "Background"
  rv_non_coding$Pattern_G_L_int   = "Background"
  rv_non_coding$Pattern_L_int     = "Background"
  threshold <- 0.05
}

col_logFC <- c(1:4)
for (i in col_logFC){

  #Upregulated per condition
  rv_non_coding[which((rv_non_coding[(i+4)] < threshold) &
                        (rv_non_coding[i] > 0)),(i+8)] <- "UP"
  #Downregulated per condition
  rv_non_coding[which(rv_non_coding[(i+4)] < threshold &
                        rv_non_coding[i] < 0),(i+8)] <- "DOWN"

}

DE_genes_vec <- c()

for (i in c(9:12)){

  print(paste("There are", length(which(rv_non_coding[,i]=="UP" |
                                          rv_non_coding[,i]=="DOWN")), "genes in",colnames(rv_non_coding)[i]))

  DE_genes_vec[i-8] <- length(which(rv_non_coding[,i]=="UP" | rv_non_coding[,i]=="DOWN"))
}

names_samples <- gsub('Pattern_','',colnames(rv_non_coding[9:12]))
# col <- c("lightcoral","red2","lightgoldenrod","yellow3","skyblue","blue4",
#          "lawngreen","darkgreen")
col <- c( "#852121", "#8f4a17","#2c682c","#355f91")


DE_genes_df <- data.frame(DE=DE_genes_vec,Contrast=names_samples,color=col)

### Save the list of DE genes in Lipids
write.table(rv_non_coding,file.path(input_dir,"txt/DE_genes_at_Interaction_RVnc_th<0.05.txt"))

rv_non_coding$RV <- rownames(rv_non_coding)
rv_non_coding <- rv_non_coding %>% relocate(RV, .before = Interaction_G_D.LogFC )
write_xlsx(rv_non_coding,file.path(input_dir,"xlsx/DE_genes_at_Interaction_RVnc_th<0.05.xlsx"))

df <- DE_genes_df
df$Contrast <- factor(df$Contrast,levels = df$Contrast)
breaks_values <- c(pretty(df$DE))

bar_plt <- ggplot(df, aes(y=DE, x=Contrast, fill= Contrast)) +
  geom_bar(stat="identity",width = 0.8)+
  scale_fill_manual(values = col)+
  scale_y_continuous(expand = c(0,0), limits = c(0,35))+
  # scale_color_manual(values = col)+
  # geom_col_pattern(data=df[5:8,],pattern="stripe",
  #                  pattern_fill = color_vec[5:8],
  #                  pattern_angle = 45,
  #                  pattern_spacing = 0.05)+
  # scale_pattern_manual(values = color_vec[5:8])+
  # coord_flip()+
  theme_classic()+
  geom_hline(yintercept=0)+
  # theme(aspect.ratio = .9)+
  labs(x="Contrasts", y = "Number of Genes")+
  theme(
    axis.text.x   = element_text(size=12),
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
bar_plt

bar_pl_1 <- bar_plt+
  geom_shadowtext(
    data = df,
    aes(x=Contrast, y=DE, label = DE),
    hjust = 0.5,vjust = -0.5,
    position = "stack",
    # nudge_x = -0.3,
    colour = "black",
    bg.colour = "white",
    bg.r = 0.2,
    # family = "Econ Sans Cnd",
    size = 5)

bar_pl_1

pdf(file=file.path(output_dir,"001_Figures_paper","barplot_DE_noncod_genes_Interaction_effects.pdf"),width=7,height=7)
print(bar_pl_1)
dev.off()


############################
# Create Heatmap iron effects
############################

# Create a matrix for the heatmap

int_effects_matrix <- as.matrix(int_effects[,c(1:4)])
# Define the color palette

color_palette <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
# Create the heatmap

columns_name <- colnames(int_effects_matrix)
column_title <- "Interaction Effects on Non-coding Genes"
clust_name <- "RV Non-coding Genes"

ht_plt <- Heatmap(int_effects_matrix,
                  na_col = "grey2",
                  col = color_palette,
                  # split = k,
                  name="LogFC",
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
                  cluster_rows = TRUE,
                  # cluster_columns = color_branches(col_dend),
                  row_title = clust_name,
                  row_title_gp = gpar(fontize = 2),
                  row_dend_side = "right",
                  row_names_side = "left",
                  # row_dend_width = unit(2, "cm"),
                  show_row_names = T ,
                  show_row_dend = T,
                  cell_fun = function(j, i, x, y, width, height, fill) {
                    grid.text(round(iron_effects[i, j],digits = 2), x, y, gp = gpar(fontsize = 3))},
                  heatmap_legend_param = list(title = "LogFC",
                                              title_position = "leftcenter-rot",
                                              labels_gp = gpar(font = 3),
                                              title_gp = gpar( fontsize = 8)))

draw(ht_plt)
# Save the heatmap to a PDF file
pdf(file=file.path(output_dir,"001_Figures_paper","Heatmap_Interaction_effects_noncoding_genes.pdf"),width=7,height=7)
draw(ht_plt)
dev.off()



####################################
##  Heatmap version 2 wo numbers  ##
####################################

# Paletas
col_fun_up   <- colorRamp2(c(0, 2), c("linen", "firebrick4"))
col_fun_down <- colorRamp2(c(-2,0), c("#1a528b","linen"))

legend_up   <- Legend(title = "Upregulated", col_fun = col_fun_up, direction = "vertical")
legend_down <- Legend(title = "Downregulated", col_fun = col_fun_down, direction = "vertical")
combined_legend <- packLegend(legend_up, legend_down, direction = "vertical")

# Matrices
logfc_mat <- as.matrix(rv_non_coding[, grep("LogFC", colnames(rv_non_coding))])
colnames(logfc_mat) <- gsub(".LogFC", "", colnames(logfc_mat))

pattern_mat <- as.matrix(rv_non_coding[, grep("^Pattern_", names(rv_non_coding))])
patterns <- as.data.frame(pattern_mat)
pattern_mat <- apply(pattern_mat, c(1,2), function(x) gsub('"', '', x))

# remove quotes from the values
pattern_mat <- apply(pattern_mat, c(1,2), function(x) gsub('"', '', x))

# Layer function
layer_fun <- function(j, i, x, y, w, h, fill) {

  mapply(function(i_, j_, x_, y_, w_, h_) {
    value  <- logfc_mat[i_, j_]
    status <- pattern_mat[i_, j_]

    # Color de fondo según patrón
    color <- if (is.na(status) || status == "Background") {
      "linen"
    } else if (status == "UP") {
      col_fun_up(abs(value))
    } else if (status == "DOWN") {
      col_fun_down((value))
    } else {
      "white"
    }

    # Dibujar celda coloreada
    grid.rect(x_, y_, w_, h_, gp = gpar(fill = color, col = "black", lwd = 0.2))

    # Añadir texto solo si es DE (UP o DOWN)
    if (!is.na(status) && status %in% c("UP", "DOWN")) {
      grid.text(sprintf("%.2f", value), x_, y_, gp = gpar(fontsize = 6))
    }

  }, i, j, x, y, w, h)
}

pt_final <- Heatmap(logfc_mat,
                    rect_gp = gpar(col = "white", lwd = 0.5),
                    name = 'logFC',
                    cluster_rows = FALSE,
                    cluster_columns = FALSE,
                    column_names_gp = gpar(fontsize = 10, fontface = "italic"),
                    row_names_gp = gpar(fontsize = 8),
                    column_names_rot = 90,
                    show_row_dend = FALSE,
                    row_names_side = "left",
                    show_column_dend = FALSE,
                    show_column_names = TRUE,
                    show_heatmap_legend = FALSE,
                    width = unit(2, "cm"),
                    height = unit(20, "cm"),
                    row_gap = unit(1, "cm"),
                    layer_fun = layer_fun)

draw(pt_final, annotation_legend_list = combined_legend)

pdf(file=file.path(output_dir,"001_Figures_paper","Heatmap_noncoding_genes_Interaction.pdf"),width=6,height=10)
draw(pt_final, annotation_legend_list = combined_legend)
dev.off()

######################################################################
##  Bar plot of DE genes in non-coding genes UP and DOWN regulated  ##
######################################################################

# 1. Convert to long format
pattern_long <- patterns %>%
  tibble::rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "Condition", values_to = "Status")

# 3. Contar genes up y down
counts <- pattern_long %>%
  filter(Status %in% c("UP", "DOWN")) %>%
  group_by(Condition, Status) %>%
  summarise(N = n(), .groups = "drop") %>%
  # Convertir "DOWN" a negativo
  mutate(N = ifelse(Status == "DOWN", -N, N))
# 4. Create the bar plot
# Etiquetas más limpias
counts$Condition <- gsub("Pattern_", "", counts$Condition)
counts$Condition <- factor(counts$Condition,
                           levels = c("G_D_int", "G_int", "G_L_int", "L_int"))
# Plot
bar_plt_1 <- ggplot(counts, aes(x = N, y = Condition, fill = Status)) +
  geom_col(width = 0.6) +
  scale_fill_manual(values = c("UP" = "firebrick4", "DOWN" = "#1a528b")) +
  geom_vline(xintercept = 0, color = "black") +
  labs(x = "Number of Genes", y = "Condition", title = "Up- and Down-regulated Genes by Condition") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")

bar_plt_1 <- ggplot(counts, aes(x = N, y = Condition, fill = Status)) +
  geom_col(width = 0.6, position = position_stack(reverse = TRUE)) +
  geom_text(aes(label = abs(N)),
            position = position_stack(vjust = 0.5, reverse = TRUE),
            size = 5, color = "white") +
  scale_fill_manual(
    values = c(
      "UP" = "firebrick4",
      "DOWN" = "#1a528b"
    )
  ) +
  geom_vline(xintercept = 0, color = "black") +
  labs(x = "Number of Genes", y = "Condition", title = "Up- and Down-regulated Genes by Condition") +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    panel.background = element_rect(fill = "linen", color = NA),
    plot.background = element_rect(fill = "linen", color = NA)
  )
bar_plt_1
# Save the bar plot

pdf(file=file.path(output_dir,"001_Figures_paper","Barplot_DE_genes_noncoding_Interactions.pdf"),width=7,height=5)
print(bar_plt_1)
dev.off()


#####################################
##  Heatmap with just DE nc genes  ##
#####################################

# Create a matrix for the heatmap with just DE nc genes

# Identificar genes con al menos un UP o DOWN en pattern_mat
rows_de <- apply(pattern_mat, 1, function(row) any(row %in% c("UP", "DOWN")))

# Filtrar las matrices
logfc_mat_filtered    <- logfc_mat[rows_de, ]
pattern_mat_filtered  <- pattern_mat[rows_de, ]

layer_fun1 <- function(j, i, x, y, w, h, fill) {

  mapply(function(i_, j_, x_, y_, w_, h_) {
    value  <- logfc_mat_filtered[i_, j_]
    status <- pattern_mat_filtered[i_, j_]

    # Color de fondo según patrón
    color <- if (is.na(status) || status == "Background") {
      "linen"
    } else if (status == "UP") {
      col_fun_up(value)
    } else if (status == "DOWN") {
      col_fun_down((value))
    } else {
      "white"
    }

    # Dibujar celda coloreada
    grid.rect(x_, y_, w_, h_, gp = gpar(fill = color, col = "black", lwd = 0.2))

    # Añadir texto solo si es DE (UP o DOWN)
    if (!is.na(status) && status %in% c("UP", "DOWN")) {
      grid.text(sprintf("%.2f", value), x_, y_, gp = gpar(fontsize = 6))
    }

  }, i, j, x, y, w, h)
}

ht_plt3 <- Heatmap(logfc_mat_filtered,
                   rect_gp = gpar(col = "white", lwd = 0.5),
                   name = 'logFC',
                   cluster_rows = T,
                   cluster_columns = FALSE,
                   column_names_gp = gpar(fontsize = 10, fontface = "italic"),
                   row_names_gp = gpar(fontsize = 8),
                   column_names_rot = 90,
                   show_row_dend = FALSE,
                   row_names_side = "left",
                   show_column_dend = FALSE,
                   show_column_names = TRUE,
                   show_heatmap_legend = FALSE,
                   width = unit(2, "cm"),
                   height = unit(20, "cm"),
                   row_gap = unit(1, "cm"),
                   layer_fun = layer_fun1)

draw(ht_plt3, annotation_legend_list = combined_legend)

pdf(file=file.path(output_dir,"001_Figures_paper","Heatmap_just_DE_noncoding_genes_Interaction.pdf"),width=6,height=10)
draw(ht_plt3, annotation_legend_list = combined_legend)
dev.off()


#################################
##  EXPONENTIAL CASE DE GENES  ##
#################################

vec_exp <- c(1:4)
LogFC_df_exp <- LogFC_df[,vec_exp]
BH_df_exp <- BH_df[,vec_exp]
lfcSE_df_exp <- lfcSE_df[,vec_exp]

### Change the columns names
colnames(LogFC_df_exp) <- paste0(substring(colnames(LogFC_df_exp),1,nchar(colnames(LogFC_df_exp))-24),".LogFC")
# colnames(BH_df_exp) <- paste0(colnames(BH_df_exp),".BH")
colnames(lfcSE_df_exp) <- paste0(substring(colnames(lfcSE_df_exp),1,nchar(colnames(lfcSE_df_exp))-9))

length(which(rownames(LogFC_df_exp)!=rownames(BH_df_exp)))
length(which(rownames(LogFC_df_exp)!=rownames(lfcSE_df_exp)))

###Create exp effects
exp_effects <- data.frame(LogFC_df_exp,BH_df_exp)

### Save the exp effects dataframe
write.table(exp_effects,file.path(input_dir,"txt/Exp_effects.txt"))
exp_effects_export <- exp_effects %>%
  rownames_to_column(var = "Gene")
write_xlsx(exp_effects_export,file.path(input_dir,"xlsx/Exp_effects.xlsx"))

rv_non_coding <- grep("RVnc", rownames(exp_effects))
rv_non_coding_custom <- rownames(reads_nc)

rv_nc_genes <- c(rownames(exp_effects)[rv_non_coding],rv_non_coding_custom)
exp_effects <- exp_effects[which(rownames(exp_effects) %in% rv_nc_genes),]
rv_non_coding <- data.frame(exp_effects)

### Save the exp effects dataframe with non-coding genes
write.table(rv_non_coding,file.path(input_dir,"txt/Exp_effects_RVnc.txt"))
rv_non_coding_export <- rv_non_coding %>%
  rownames_to_column(var = "Gene")
write_xlsx(rv_non_coding_export,file.path(input_dir,"xlsx/Exp_effects_RVnc.xlsx"))

{
  rv_non_coding$Pattern_G_D_stat   = "Background"
  rv_non_coding$Pattern_G_stat     = "Background"
  rv_non_coding$Pattern_G_L_stat   = "Background"
  rv_non_coding$Pattern_L_stat     = "Background"
  threshold <- 0.05
}

col_logFC <- c(1:4)
for (i in col_logFC){

  #Upregulated per condition
  rv_non_coding[which((rv_non_coding[(i+4)] < threshold) &
                        (rv_non_coding[i] > 0)),(i+8)] <- "UP"
  #Downregulated per condition
  rv_non_coding[which(rv_non_coding[(i+4)] < threshold &
                        rv_non_coding[i] < 0),(i+8)] <- "DOWN"

}

DE_genes_vec <- c()

for (i in c(9:12)){

  print(paste("There are", length(which(rv_non_coding[,i]=="UP" |
                                          rv_non_coding[,i]=="DOWN")), "genes in",colnames(rv_non_coding)[i]))

  DE_genes_vec[i-8] <- length(which(rv_non_coding[,i]=="UP" | rv_non_coding[,i]=="DOWN"))
}

names_samples <- gsub('Pattern_','',colnames(rv_non_coding[9:12]))
# col <- c("lightcoral","red2","lightgoldenrod","yellow3","skyblue","blue4",
#          "lawngreen","darkgreen")
col <- c( "#852121", "#8f4a17","#2c682c","#355f91")


DE_genes_df <- data.frame(DE=DE_genes_vec,Contrast=names_samples,color=col)

# ### Save the list of DE genes in Lipids
# write.table(rv_non_coding,file.path(input_dir,"txt/DE_genes_at_stat_RVnc_th<0.05.txt"))
#
# rv_non_coding$RV <- rownames(rv_non_coding)
# rv_non_coding <- rv_non_coding %>% relocate(RV, .before = exp_effect_in_STAT_G_D.LogFC )
# write_xlsx(rv_non_coding,file.path(input_dir,"xlsx/DE_genes_at_stat_RVnc_th<0.05.xlsx"))
#
# df <- DE_genes_df
# df$Contrast <- factor(df$Contrast,levels = df$Contrast)
# breaks_values <- c(pretty(df$DE))
#
# bar_plt <- ggplot(df, aes(y=DE, x=Contrast, fill= Contrast)) +
#   geom_bar(stat="identity",width = 0.8)+
#   scale_fill_manual(values = col)+
#   scale_y_continuous(expand = c(0,0), limits = c(0,35))+
#   # scale_color_manual(values = col)+
#   # geom_col_pattern(data=df[5:8,],pattern="stripe",
#   #                  pattern_fill = color_vec[5:8],
#   #                  pattern_angle = 45,
#   #                  pattern_spacing = 0.05)+
#   # scale_pattern_manual(values = color_vec[5:8])+
#   # coord_flip()+
#   theme_classic()+
#   geom_hline(yintercept=0)+
#   # theme(aspect.ratio = .9)+
#   labs(x="Contrasts", y = "Number of Genes")+
#   theme(
#     axis.text.x   = element_text(size=12),
#     axis.text.y   = element_text(size=12),
#     axis.title.x  = element_text(size=12),
#     axis.title.y  = element_text(size=12),
#     # axis.ticks.y = element_blank(),
#     #panel.background = element_blank(),
#     #panel.grid.major = element_blank(),
#     #panel.grid.minor = element_blank(),
#     #axis.line = element_line(colour = "black"),
#     panel.border = element_rect(colour = "black", fill=NA, linewidth=1,linetype="solid"),
#     legend.title=element_blank(),
#     legend.position="top",
#     #legend.text=element_text(size=14),
#     legend.key.size = unit(1, 'lines'))
# bar_plt
#
# bar_pl_1 <- bar_plt+
#   geom_shadowtext(
#     data = df,
#     aes(x=Contrast, y=DE, label = DE),
#     hjust = 0.5,vjust = -0.5,
#     position = "stack",
#     # nudge_x = -0.3,
#     colour = "black",
#     bg.colour = "white",
#     bg.r = 0.2,
#     # family = "Econ Sans Cnd",
#     size = 5)
#
# bar_pl_1
#
# pdf(file=file.path(output_dir,"001_Figures_paper","barplot_DE_noncod_genes_stat_effects.pdf"),width=7,height=7)
# print(bar_pl_1)
# dev.off()
#
# ###############################################################################
# ##  Bar plot of DE genes in non-coding genes UP and DOWN regulated  EXP case ##
# ###############################################################################
#
# # 1. Convert to long format
# pattern_long <- patterns %>%
#   tibble::rownames_to_column("Gene") %>%
#   pivot_longer(-Gene, names_to = "Condition", values_to = "Status")
#
# # 3. Contar genes up y down
# counts <- pattern_long %>%
#   filter(Status %in% c("UP", "DOWN")) %>%
#   group_by(Condition, Status) %>%
#   summarise(N = n(), .groups = "drop") %>%
#   # Convertir "DOWN" a negativo
#   mutate(N = ifelse(Status == "DOWN", -N, N))
# # 4. Create the bar plot
# # Etiquetas más limpias
# counts$Condition <- gsub("Pattern_", "", counts$Condition)
# counts$Condition <- factor(counts$Condition,
#                            levels = c("G_D_stat", "G_stat", "G_L_stat", "L_stat"))
# # Plot
# bar_plt_1 <- ggplot(counts, aes(x = N, y = Condition, fill = Status)) +
#   geom_col(width = 0.6) +
#   scale_fill_manual(values = c("UP" = "firebrick4", "DOWN" = "#1a528b")) +
#   geom_vline(xintercept = 0, color = "black") +
#   labs(x = "Number of Genes", y = "Condition", title = "Up- and Down-regulated Genes by Condition") +
#   theme_minimal(base_size = 13) +
#   theme(legend.position = "top")
#
# bar_plt_1 <- ggplot(counts, aes(x = N, y = Condition, fill = Status)) +
#   geom_col(width = 0.6, position = position_stack(reverse = TRUE)) +
#   geom_text(aes(label = abs(N)),
#             position = position_stack(vjust = 0.5, reverse = TRUE),
#             size = 5, color = "white") +
#   scale_fill_manual(
#     values = c(
#       "UP" = "firebrick4",
#       "DOWN" = "#1a528b"
#     )
#   ) +
#   geom_vline(xintercept = 0, color = "black") +
#   labs(x = "Number of Genes", y = "Condition", title = "Up- and Down-regulated Genes by Condition") +
#   theme_minimal(base_size = 13) +
#   theme(
#     legend.position = "top",
#     panel.background = element_rect(fill = "linen", color = NA),
#     plot.background = element_rect(fill = "linen", color = NA)
#   )
# bar_plt_1
# # Save the bar plot
#
# pdf(file=file.path(output_dir,"001_Figures_paper","Barplot_DE_genes_noncoding.pdf"),width=7,height=5)
# print(bar_plt_1)
# dev.off()
#
