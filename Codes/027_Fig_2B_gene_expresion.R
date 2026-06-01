#############################
### 0. Load dependencies. ###
#############################

{
  library(tidyverse)
  library(ggplot2)
  library(DESeq2)
}

############################
### 1. Declare functions ###
############################

plot_one_gene_all_panels <- function(gene, expr_mat, sample_info, outdir = NULL) {

  gene_class_i <- gene_class[gene]
  panel_col    <- class_colors[gene_class_i]

  if (!gene %in% rownames(expr_mat)) {
    stop(paste("Gene not found:", gene))
  }

  df <- data.frame(
    Sample = colnames(expr_mat),
    Expression = as.numeric(expr_mat[gene, colnames(expr_mat)]),
    stringsAsFactors = FALSE
  )

  df <- merge(df, sample_info, by = "Sample")

  carbon_levels <- c("G-D", "G", "G-LCFA", "LCFA")
  fe_levels     <- c("Fe+", "Fe-")

  plot_list <- list()

  for (carb in carbon_levels) {

    df_carb <- df[df$Carbon == carb, ]

    ymin_pair   <- min(df_carb$Expression, na.rm = TRUE)
    ymax_pair   <- max(df_carb$Expression, na.rm = TRUE)
    breaks_pair <- pretty(c(ymin_pair, ymax_pair), n = 4)

    for (fe_status in fe_levels) {

      df_sub <- df_carb[df_carb$Fe == fe_status, ]

      p <- ggplot(df_sub, aes(x = Phase, y = Expression,
                              color = paste(Carbon, Fe, sep = "_"),
                              fill  = paste(Carbon, Fe, sep = "_"))) +
        geom_point(
          shape = 21,
          size = 3.5,
          stroke = 1,
          position = position_jitter(width = 0.08, height = 0)
        ) +
        scale_color_manual(values = carbon_fe_colors, drop = FALSE) +
        scale_fill_manual(values = carbon_fe_fill, drop = FALSE) +
        guides(color = "none", fill = "none") +
        scale_y_continuous(
          limits = c(ymin_pair, ymax_pair),
          breaks = breaks_pair,
          labels = breaks_pair
        ) +
        theme_classic() +
        labs(
          title = paste(gene, "-", carb, "-", fe_status),
          x = NULL,
          y = "VST expression"
        ) +
        theme(
          plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(colour = panel_col, fill = NA, linewidth = 1.5),
          axis.line = element_blank(),
          legend.position = "none"
        )

      fe_clean  <- gsub("\\+", "plus", fe_status)
      fe_clean  <- gsub("-", "minus", fe_clean)
      plot_name <- paste(gene, carb, fe_clean, sep = "_")

      plot_list[[plot_name]] <- p

      if (!is.null(outdir)) {
        pdf(file = file.path(outdir, paste0(plot_name, ".pdf")),
            width = 3.5,
            height = 3.5)
        print(p)
        dev.off()
      }
    }
  }

  return(plot_list)
}

### For the paired plots
plot_one_gene_paired_panels <- function(gene, expr_mat, sample_info, outdir = NULL) {

  gene_class_i <- gene_class[gene]
  panel_col    <- class_colors[gene_class_i]

  if (!gene %in% rownames(expr_mat)) {
    stop(paste("Gene not found:", gene))
  }

  df <- data.frame(
    Sample = colnames(expr_mat),
    Expression = as.numeric(expr_mat[gene, colnames(expr_mat)]),
    stringsAsFactors = FALSE
  )

  # carbon_panel_colors <- c(
  #   "G-D"    = "#5f4290",
  #   "G"      = "#3c8ec8",
  #   "G-LCFA" = "#559d3f",
  #   "LCFA"   = "#d1382c"
  # )

  df <- merge(df, sample_info, by = "Sample")

  carbon_levels <- c("G-D", "G", "G-LCFA", "LCFA")
  pair_list <- list()

  for (carb in carbon_levels) {

    df_carb <- df[df$Carbon == carb, ]

    # panel_col <- carbon_panel_colors[carb]

    ymin_pair   <- min(df_carb$Expression, na.rm = TRUE)
    ymax_pair   <- max(df_carb$Expression, na.rm = TRUE)
    breaks_pair <- pretty(c(ymin_pair, ymax_pair), n = 4)

    df_fe_plus  <- df_carb[df_carb$Fe == "Fe+", ]
    df_fe_minus <- df_carb[df_carb$Fe == "Fe-", ]

    p_left <- ggplot(df_fe_plus, aes(x = Phase, y = Expression,
                                     color = paste(Carbon, Fe, sep = "_"),
                                     fill  = paste(Carbon, Fe, sep = "_"))) +
      geom_point(
        shape = 21,
        size = 3.5,
        stroke = 1,
        position = position_jitter(width = 0.08, height = 0)
      ) +
      scale_color_manual(values = carbon_fe_colors, drop = FALSE) +
      scale_fill_manual(values = carbon_fe_fill, drop = FALSE) +
      guides(color = "none", fill = "none") +
      scale_y_continuous(
        limits = c(ymin_pair, ymax_pair),
        breaks = breaks_pair,
        labels = breaks_pair
      ) +
      theme_classic() +
      labs(
        title = paste(gene, "-", carb, "- Fe+"),
        x = NULL,
        y = "VST expression"
      ) +
      theme(
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = panel_col, fill = NA, linewidth = 1.5),
        axis.line = element_blank(),
        legend.position = "none"
      )

    p_right <- ggplot(df_fe_minus, aes(x = Phase, y = Expression,
                                       color = paste(Carbon, Fe, sep = "_"),
                                       fill  = paste(Carbon, Fe, sep = "_"))) +
      geom_point(
        shape = 21,
        size = 3.5,
        stroke = 1,
        position = position_jitter(width = 0.08, height = 0)
      ) +
      scale_color_manual(values = carbon_fe_colors, drop = FALSE) +
      scale_fill_manual(values = carbon_fe_fill, drop = FALSE) +
      guides(color = "none", fill = "none") +
      scale_y_continuous(
        limits = c(ymin_pair, ymax_pair),
        breaks = breaks_pair,
        labels = breaks_pair
      ) +
      theme_classic() +
      labs(
        title = paste(gene, "-", carb, "- Fe-"),
        x = NULL,
        y = NULL
      ) +
      theme(
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = panel_col, fill = NA, linewidth = 1.5),
        axis.line = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none"
      )

    paired_plot <- p_left + p_right + plot_layout(ncol = 2)

    pair_name <- paste(gene, carb, "paired", sep = "_")
    pair_list[[pair_name]] <- paired_plot

    if (!is.null(outdir)) {
      pdf(file = file.path(outdir, paste0(pair_name, ".pdf")),
          width = 7,
          height = 3.5)
      print(paired_plot)
      dev.off()
    }
  }

  return(pair_list)
}
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

reads    <- read.table(file.path(input_dir, "txt/reads_32_samples.txt"), check.names = FALSE)
metadata <- read.table(file.path(input_dir, "txt/metadata_32_samples.txt"), check.names = FALSE)

#####################################################
### 4. Build DESeq2 object and obtain VST matrix ####
#####################################################

# Reorder metadata to match count matrix columns
metadata_dds <- metadata[match(colnames(reads), rownames(metadata)), , drop = FALSE]

# Safety check
stopifnot(all(rownames(metadata_dds) == colnames(reads)))

# Correct grouping variable
metadata_dds$group <- factor(
  metadata_dds$Culture,
  levels = c("C1", "C2", "C5", "C6", "C7", "C8", "C11", "C12")
)

metadata_dds$Iron <- factor(
  metadata_dds$Iron,
  levels = c("YES", "NO")
)

group <- factor(metadata_dds$short_setup)

### Design matrix
design=model.matrix(~0+group,data=factor(metadata$short_setup))
colnames(design) <- levels(group)

# Build DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = as.matrix(reads),
  colData   = metadata_dds,
  design    = design
)

# VST for visualization
vsd     <- vst(dds, blind = FALSE)
vst_mat <- assay(vsd)

########################################
### 5. Build plotting sample metadata ###
########################################

sample_info <- data.frame(
  Sample = colnames(vst_mat),
  stringsAsFactors = FALSE
)

sample_info$Condition <- sub("_Fe_.*", "", sample_info$Sample)
sample_info$Fe  <- ifelse(grepl("_Fe_YES_", sample_info$Sample), "Fe+", "Fe-")
sample_info$Rep <- sub(".*_(\\d+)$", "\\1", sample_info$Sample)

sample_info$Carbon <- NA
sample_info$Phase  <- NA

sample_info$Carbon[sample_info$Condition %in% c("C1", "C2")]   <- "G-D"
sample_info$Carbon[sample_info$Condition %in% c("C5", "C6")]   <- "G"
sample_info$Carbon[sample_info$Condition %in% c("C7", "C8")]   <- "G-LCFA"
sample_info$Carbon[sample_info$Condition %in% c("C11", "C12")] <- "LCFA"

sample_info$Phase[sample_info$Condition %in% c("C1", "C5", "C7", "C11")] <- "EXP"
sample_info$Phase[sample_info$Condition %in% c("C2", "C6", "C8", "C12")] <- "STAT"

sample_info$Carbon <- factor(sample_info$Carbon, levels = c("G-D", "G", "G-LCFA", "LCFA"))
sample_info$Phase  <- factor(sample_info$Phase, levels = c("EXP", "STAT"))
sample_info$Fe     <- factor(sample_info$Fe, levels = c("Fe+", "Fe-"))

sample_info$Carbon_Fe <- paste(sample_info$Carbon, sample_info$Fe, sep = "_")

##########################################
### 6. Define genes and panel colors #####
##########################################

# Rv3146, Rv3149, Rv3150, Rv3153, Rv3157 <- ED
# Rv1994c, Rv2641 <- EU
# Rv3030 <- DD
# Rv1617 Rv2430c Rv2431c Rv3060c Rv3061c <- DU

class_colors <- c(
  "DU" = "#FFEFD9",
  "DD" = "#A59FCB",
  "EU" = "#FFA373",
  "ED" = "#50486D"
)

gene_class <- c(
  "Rv3146"  = "ED",
  "Rv3149"  = "ED",
  "Rv3150"  = "ED",
  "Rv3153"  = "ED",
  "Rv3157"  = "ED",
  "Rv1994c" = "EU",
  "Rv2641"  = "EU",
  "Rv3030"  = "DD",
  "Rv1617"  = "DU",
  "Rv2430c" = "DU",
  "Rv2431c" = "DU",
  "Rv3060c" = "DU",
  "Rv3061c" = "DU"
)

genes_of_interest <- c(
  "Rv3146", "Rv3149", "Rv3150", "Rv3153", "Rv3157",
  "Rv1994c", "Rv2641",
  "Rv3030",
  "Rv1617", "Rv2430c", "Rv2431c", "Rv3060c", "Rv3061c"
)

##################################
### 7. Define point aesthetics ###
##################################

carbon_fe_colors <- c(
  "G-D_Fe+"    = "#d3c7d9",
  "G-D_Fe-"    = "#5f4290",
  "G_Fe+"      = "#76c6e7",
  "G_Fe-"      = "#3c8ec8",
  "G-LCFA_Fe+" = "#bbdd93",
  "G-LCFA_Fe-" = "#559d3f",
  "LCFA_Fe+"   = "#ee9f9b",
  "LCFA_Fe-"   = "#d1382c"
)

carbon_fe_fill <- c(
  "G-D_Fe+"    = "#d3c7d9",
  "G-D_Fe-"    = "white",
  "G_Fe+"      = "#76c6e7",
  "G_Fe-"      = "white",
  "G-LCFA_Fe+" = "#bbdd93",
  "G-LCFA_Fe-" = "white",
  "LCFA_Fe+"   = "#ee9f9b",
  "LCFA_Fe-"   = "white"
)

############################
### 8. Generate panels #####
############################

dir.create(file.path(output_dir, "Gene_panels_Fig_2B"),
           recursive = TRUE, showWarnings = FALSE)

all_plots <- lapply(genes_of_interest, function(g) {
  plot_one_gene_all_panels(
    gene = g,
    expr_mat = vst_mat,
    sample_info = sample_info,
    outdir = file.path(output_dir, "Gene_panels_Fig_2B")
  )
})

names(all_plots) <- genes_of_interest

##########################
### 9. Quick checks ######
##########################

print(table(metadata_dds$group))
print(any(is.na(metadata_dds$group)))
print(all(colnames(vst_mat) == sample_info$Sample))
print(unique(sample_info$Carbon_Fe))

all_plots["Rv3030"]

### Paired plots for Fe+ vs Fe- within each carbon source
library(patchwork)
dir.create(file.path(output_dir, "Gene_panels_Fig_2B_paired"),
           recursive = TRUE, showWarnings = FALSE)

all_paired_plots <- lapply(genes_of_interest, function(g) {
  plot_one_gene_paired_panels(
    gene = g,
    expr_mat = vst_mat,
    sample_info = sample_info,
    outdir = file.path(output_dir, "Gene_panels_Fig_2B_paired")
  )
})

names(all_paired_plots) <- genes_of_interest
all_paired_plots["Rv3030"]
