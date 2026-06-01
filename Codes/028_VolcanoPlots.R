#############################
### 0. Load dependencies. ###
#############################
{
  library(readxl)
  library(writexl)
  library(tidyverse)
  library(shadowtext)
  # library(xlsx)
  library(openxlsx)
  library(ggrepel)
  library(patchwork)

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
input_dir <- "Inputs/002_Processed_data"
output_dir <- "Outputs"

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
BH_df    <- read.table(file.path(input_dir,"txt/BH_0.05.txt"))
lfcSE_df <- read.table(file.path(input_dir,"txt/lfcSE_0.05.txt"))

### Extract the 4 interactions
# 17 Effect_of_Fe_in_response_to_Growth_arrest_G_D.log2FoldChange_shrunken
# 18   Effect_of_Fe_in_response_to_Growth_arrest_G.log2FoldChange_shrunken
# 19 Effect_of_Fe_in_response_to_Growth_arrest_G_L.log2FoldChange_shrunken
# 20   Effect_of_Fe_in_response_to_Growth_arrest_L.log2FoldChange_shrunken

vec_int <- c(17:20)
LogFC_df_int <- LogFC_df[,vec_int]
BH_df_int    <- BH_df[,vec_int]
lfcSE_df_int <- lfcSE_df[,vec_int]

colnames(LogFC_df_int) <- paste0(substring(colnames(LogFC_df_int),1,nchar(colnames(LogFC_df_int))-24),".LogFC")
# colnames(BH_df_int) <- paste0(colnames(int_BH),".BH")
colnames(lfcSE_df_int) <- paste0(substring(colnames(lfcSE_df_int),1,nchar(colnames(lfcSE_df_int))-9))

length(which(rownames(LogFC_df_int)!=rownames(BH_df_int)))
length(which(rownames(LogFC_df_int)!=rownames(lfcSE_df_int)))


# I want to create volcano plots for each of the 4 interactions, using the logFC, BH


for (i in 1:ncol(LogFC_df_int)) {
  # i=1
  logFC_col <- colnames(LogFC_df_int)[i]
  BH_col    <- paste0(substring(logFC_col,1,nchar(logFC_col)-6),".BH")
  lfcSE_col <- paste0(substring(logFC_col,1,nchar(logFC_col)-6),".lfcSE")

  df_plot <- data.frame(feature=rownames(LogFC_df_int),
                        logFC=LogFC_df_int[,logFC_col],
                        BH=BH_df_int[,BH_col],
                        lfcSE=lfcSE_df_int[,lfcSE_col])

  # Create the volcano plot using ggplot2
  p <- ggplot(df_plot, aes(x = logFC, y = -log10(BH))) +
    geom_point() +
    theme_minimal() +
    labs(title = paste("Volcano Plot for", logFC_col),
         x = "Log Fold Change",
         y = "-Log10(BH Adjusted P-value)")
  p
  # Save the plot as a PNG file
# pdf(file.path(output_dir,paste0("Volcano_plot_",logFC_col,".pdf")), width = 8, height = 6)
# print(p)
# dev.off()
}

#####
dir.create(file.path(output_dir, "Volcano_plots"), showWarnings = FALSE, recursive = TRUE)

for (i in 1:ncol(LogFC_df_int)) {

  logFC_col <- colnames(LogFC_df_int)[i]
  BH_col    <- paste0(substring(logFC_col, 1, nchar(logFC_col) - 6), ".BH")
  lfcSE_col <- paste0(substring(logFC_col, 1, nchar(logFC_col) - 6), ".lfcSE")

  df_plot <- data.frame(
    feature = rownames(LogFC_df_int),
    logFC   = LogFC_df_int[, logFC_col],
    BH      = BH_df_int[, BH_col],
    lfcSE   = lfcSE_df_int[, lfcSE_col],
    stringsAsFactors = FALSE
  )

  # Remove NA
  df_plot <- df_plot[!is.na(df_plot$logFC) & !is.na(df_plot$BH), ]

  # Avoid Inf in -log10(BH)
  min_nonzero_bh <- min(df_plot$BH[df_plot$BH > 0], na.rm = TRUE)
  df_plot$BH_plot <- ifelse(df_plot$BH == 0, min_nonzero_bh, df_plot$BH)
  df_plot$minus_log10_BH <- -log10(df_plot$BH_plot)

  # Classification
  df_plot$category <- "Not significant"
  df_plot$category[df_plot$BH < 0.05 & df_plot$logFC > 1]  <- "Upregulated"
  df_plot$category[df_plot$BH < 0.05 & df_plot$logFC < -1] <- "Downregulated"

  df_plot$category <- factor(
    df_plot$category,
    levels = c("Downregulated", "Not significant", "Upregulated")
  )

  # Optional: label top genes
  df_plot$label <- ""
  top_up <- head(df_plot[df_plot$category == "Upregulated", ][order(df_plot[df_plot$category == "Upregulated", ]$BH), "feature"], 8)
  top_down <- head(df_plot[df_plot$category == "Downregulated", ][order(df_plot[df_plot$category == "Downregulated", ]$BH), "feature"], 8)
  genes_to_label <- c(top_up, top_down)
  df_plot$label[df_plot$feature %in% genes_to_label] <- df_plot$feature[df_plot$feature %in% genes_to_label]

  p <- ggplot(df_plot, aes(x = logFC, y = minus_log10_BH)) +
    geom_point(aes(color = category), size = 1.8, alpha = 0.8) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.5, color = "grey40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.5, color = "grey40") +
    geom_text_repel(
      aes(label = label),
      size = 3,
      max.overlaps = 20,
      box.padding = 0.3,
      point.padding = 0.2,
      segment.linewidth = 0.3,
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = c(
        "Downregulated" = "#3C8EC8",
        "Not significant" = "grey75",
        "Upregulated" = "#D1382C"
      )
    ) +
    labs(
      title = gsub("\\.LogFC$", "", logFC_col),
      x = expression(log[2]~fold~change),
      y = expression(-log[10]~adjusted~italic(P))
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 10),
      legend.title = element_blank(),
      legend.position = "top",
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
      axis.line = element_blank()
    )

  pdf(file.path(output_dir, "Volcano_plots", paste0("Volcano_plot_", logFC_col, ".pdf")),
      width = 6.5, height = 5.5)
  print(p)
  dev.off()
}

#### Volcano plots with shared axes and legend in a 2x2 panel
dir.create(file.path(output_dir, "Volcano_plots_shared"), showWarnings = FALSE, recursive = TRUE)

global_x_lim <- max(abs(as.matrix(LogFC_df_int)), na.rm = TRUE)
x_limits <- c(-global_x_lim, global_x_lim)
x_breaks <- pretty(x_limits, n = 5)

global_y_lim <- max(-log10(as.matrix(BH_df_int[BH_df_int > 0])), na.rm = TRUE)
y_limits <- c(0, global_y_lim)
y_breaks <- pretty(y_limits, n = 5)
volcano_list <- list()

for (i in 1:ncol(LogFC_df_int)) {

  logFC_col <- colnames(LogFC_df_int)[i]
  BH_col    <- paste0(substring(logFC_col, 1, nchar(logFC_col) - 6), ".BH")
  lfcSE_col <- paste0(substring(logFC_col, 1, nchar(logFC_col) - 6), ".lfcSE")

  df_plot <- data.frame(
    feature = rownames(LogFC_df_int),
    logFC   = LogFC_df_int[, logFC_col],
    BH      = BH_df_int[, BH_col],
    lfcSE   = lfcSE_df_int[, lfcSE_col],
    stringsAsFactors = FALSE
  )

  # Remove missing values
  df_plot <- df_plot[!is.na(df_plot$logFC) & !is.na(df_plot$BH), ]

  # Avoid Inf in -log10(BH)
  min_nonzero_bh <- min(df_plot$BH[df_plot$BH > 0], na.rm = TRUE)
  df_plot$BH_plot <- ifelse(df_plot$BH == 0, min_nonzero_bh, df_plot$BH)
  df_plot$minus_log10_BH <- -log10(df_plot$BH_plot)

  # Classification
  df_plot$category <- "Not significant"
  df_plot$category[df_plot$BH < 0.05 & df_plot$logFC > 1]  <- "Upregulated"
  df_plot$category[df_plot$BH < 0.05 & df_plot$logFC < -1] <- "Downregulated"

  df_plot$category <- factor(
    df_plot$category,
    levels = c("Downregulated", "Not significant", "Upregulated")
  )

  # Optional labels
  df_plot$label <- ""
  top_up <- head(
    df_plot[df_plot$category == "Upregulated", ][order(df_plot[df_plot$category == "Upregulated", ]$BH), "feature"],
    8
  )
  top_down <- head(
    df_plot[df_plot$category == "Downregulated", ][order(df_plot[df_plot$category == "Downregulated", ]$BH), "feature"],
    8
  )

  genes_to_label <- c(top_up, top_down)
  df_plot$label[df_plot$feature %in% genes_to_label] <- df_plot$feature[df_plot$feature %in% genes_to_label]

  # Cleaner title
  plot_title <- gsub("\\.LogFC$", "", logFC_col)
  plot_title <- gsub("Effect_of_Fe_in_response_to_Growth_arrest_", "", plot_title)

  p <- ggplot(df_plot, aes(x = logFC, y = minus_log10_BH)) +
    geom_point(aes(color = category), size = 1.8, alpha = 0.8) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.5, color = "grey40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.5, color = "grey40") +
    geom_text_repel(
      aes(label = label),
      size = 3,
      max.overlaps = 20,
      box.padding = 0.3,
      point.padding = 0.2,
      segment.linewidth = 0.3,
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = c(
        "Downregulated" = "#3C8EC8",
        "Not significant" = "grey75",
        "Upregulated" = "#D1382C"
      )
    ) +
    scale_x_continuous(
      limits = x_limits,
      breaks = x_breaks
    ) +
    scale_y_continuous(
      limits = y_limits,
      breaks = y_breaks
    ) +
    labs(
      title = plot_title,
      x = expression(log[2]~fold~change),
      y = expression(-log[10]~adjusted~italic(P))
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9),
      legend.title = element_blank(),
      legend.position = "top",
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
      axis.line = element_blank()
    )

  volcano_list[[i]] <- p
}

# Remove repeated axis info
volcano_list[[2]] <- volcano_list[[2]] +
  labs(y = NULL) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

volcano_list[[4]] <- volcano_list[[4]] +
  labs(y = NULL) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

volcano_list[[1]] <- volcano_list[[1]] +
  labs(x = NULL) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

volcano_list[[2]] <- volcano_list[[2]] +
  labs(x = NULL) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# Combine 2x2 with shared legend
volcano_panel <- (volcano_list[[1]] + volcano_list[[2]]) /
  (volcano_list[[3]] + volcano_list[[4]]) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")

# Save
pdf(file.path(output_dir, "Volcano_plots_shared", "Volcano_plots_2x2_panel.pdf"),
    width = 11, height = 8.5)
print(volcano_panel)
dev.off()

volcano_panel

### Save in svg format for better quality

svglite::svglite(
  file.path(output_dir, "Volcano_plots_shared", "Volcano_plots_2x2_panel.svg"),
  width = 11,
  height = 8.5
)
print(volcano_panel)
dev.off()


# Keep y-axis in all plots, remove repeated x-axis except bottom plot
volcano_list_col <- volcano_list

volcano_list_col[[1]] <- volcano_list_col[[1]] +
  labs(x = NULL) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

volcano_list_col[[2]] <- volcano_list_col[[2]] +
  labs(x = NULL) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

volcano_list_col[[3]] <- volcano_list_col[[3]] +
  labs(x = NULL) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# Combine in one column
volcano_panel_1col <- volcano_list_col[[1]] /
  volcano_list_col[[2]] /
  volcano_list_col[[3]] /
  volcano_list_col[[4]] +
  plot_layout(guides = "collect") &

  theme(legend.position = "top")

# Save PDF: tall and narrow
pdf(
  file.path(output_dir, "Volcano_plots_shared", "Volcano_plots_1col_panel.pdf"),
  width = 4,
  height = 14
)
print(volcano_panel_1col)
dev.off()

# Save SVG
svglite::svglite(
  file.path(output_dir, "Volcano_plots_shared", "Volcano_plots_1col_panel.svg"),
  width = 4,
  height = 14
)

print(volcano_panel_1col)
dev.off()
volcano_panel_1col
