######################
### 0. Dependencies ###
######################

library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(tidyverse)
library(shadowtext)
library(grid)

dcols <- function(x) data.frame(colnames(x))


##########################
### 1. Helper Functions ###
##########################

calc_ht_size <- function(ht, unit = "inch") {
  pdf(NULL)
  ht <- draw(ht)
  w  <- convertX(ComplexHeatmap:::width(ht),  unit, valueOnly = TRUE)
  h  <- convertY(ComplexHeatmap:::height(ht), unit, valueOnly = TRUE)
  dev.off()
  c(w, h)
}

make_row_labels <- function(genes, feature_data) {
  sym <- feature_data[genes, "symbol"]
  ifelse(
    !is.na(sym) & nchar(trimws(sym)) > 0,
    trimws(sym),
    genes
  )
}

extract_go_genes <- function(go_list, terms, exact = TRUE) {

  genes <- c()

  for (nm in names(go_list)) {

    df <- go_list[[nm]]

    df_sub <- if (exact) {
      df[df$description %in% terms, , drop = FALSE]
    } else {
      df[grepl(paste(terms, collapse = "|"), df$description, ignore.case = TRUE), , drop = FALSE]
    }

    if (nrow(df_sub) > 0) {
      genes_i <- trimws(unlist(strsplit(df_sub$genes_in_query, ",")))
      genes <- c(genes, genes_i)
    }
  }

  unique(genes)
}

prepare_growth_arrest_matrix <- function(genes, LogFC_df_iron) {

  genes <- unique(intersect(genes, rownames(LogFC_df_iron)))

  if (length(genes) == 0) {
    stop("None of the selected genes are present in LogFC_df_iron.")
  }

  mat <- as.matrix(LogFC_df_iron[genes, , drop = FALSE])

  # Columns 1–4 = Fe+, columns 5–8 = Fe-
  # Reorder as Fe+/Fe- pairs per carbon source
  mat <- mat[, c(1, 5, 2, 6, 3, 7, 4, 8), drop = FALSE]

  colnames(mat) <- c(
    "G-D_Fe+", "G-D_Fe-",
    "G_Fe+", "G_Fe-",
    "G-LCFA_Fe+", "G-LCFA_Fe-",
    "LCFA_Fe+", "LCFA_Fe-"
  )

  mat
}

make_membership_matrix <- function(genes, enrichment_list) {

  genes_per_contrast <- lapply(enrichment_list, function(df) {
    unique(trimws(unlist(strsplit(df$genes_in_query, ","))))
  })

  col_names <- c(
    "G-D_Fe+", "G-D_Fe-",
    "G_Fe+", "G_Fe-",
    "G-LCFA_Fe+", "G-LCFA_Fe-",
    "LCFA_Fe+", "LCFA_Fe-"
  )

  membership_mat <- matrix(
    FALSE,
    nrow = length(genes),
    ncol = length(col_names),
    dimnames = list(genes, col_names)
  )

  if (length(genes_per_contrast) >= 1) {
    membership_mat[genes %in% genes_per_contrast[[1]], c(1, 2)] <- TRUE
  }
  if (length(genes_per_contrast) >= 2) {
    membership_mat[genes %in% genes_per_contrast[[2]], c(3, 4)] <- TRUE
  }
  if (length(genes_per_contrast) >= 3) {
    membership_mat[genes %in% genes_per_contrast[[3]], c(5, 6)] <- TRUE
  }

  membership_mat
}

make_anno_vec <- function(x) ifelse(x, "Presence", "Absence")

make_term_annotation <- function(genes,
                                 term_gene_list,
                                 annotation_fontsize = 8,
                                 anno_size = unit(4, "mm")) {

  anno_df <- as.data.frame(
    lapply(term_gene_list, function(g) genes %in% g),
    row.names = genes
  )

  anno_df[] <- lapply(anno_df, make_anno_vec)

  anno_col_map <- c("Presence" = "black", "Absence" = "white")

  rowAnnotation(
    df = anno_df,
    col = setNames(
      replicate(ncol(anno_df), anno_col_map, simplify = FALSE),
      colnames(anno_df)
    ),
    annotation_name_gp = gpar(fontsize = annotation_fontsize),
    simple_anno_size = anno_size,
    border = TRUE,
    show_legend = FALSE
  )
}

make_presence_legend <- function() {
  Legend(
    labels = c("Absence", "Presence"),
    title  = "",
    graphics = list(
      function(x, y, w, h) grid.rect(x, y, w, h, gp = gpar(fill = "white", col = "black", lwd = 1)),
      function(x, y, w, h) grid.rect(x, y, w, h, gp = gpar(fill = "black", col = "black", lwd = 1))
    )
  )
}

plot_gene_response_heatmap <- function(mat,
                                       title         = "genes",
                                       high_col      = "#cb181d",
                                       mid_col       = "#fcae91",
                                       low_col       = "white",
                                       cluster_rows  = TRUE,
                                       show_numbers  = FALSE,
                                       row_fontsize  = 11,
                                       col_fontsize  = 8,
                                       row_labels    = NULL,
                                       cap           = NULL,
                                       highlight_mat = NULL,
                                       highlight_col = "#FFA373",
                                       highlight_lwd = 3,
                                       right_annotation = NULL) {

  mat_plot <- mat
  mat_plot[!is.finite(mat_plot)] <- NA

  max_abs <- max(abs(mat_plot), na.rm = TRUE)
  if (!is.null(cap)) max_abs <- min(max_abs, cap)
  if (!is.finite(max_abs) || max_abs == 0) max_abs <- 1

  col_fun <- circlize::colorRamp2(
    c(0, max_abs / 2, max_abs),
    c(low_col, mid_col, high_col)
  )

  if (is.null(row_labels)) row_labels <- rownames(mat_plot)

  ht <- Heatmap(
    mat_plot,
    name              = "log2FC",
    col               = col_fun,
    na_col            = "white",
    column_order      = colnames(mat_plot),
    column_split      = factor(
      c("G-D", "G-D", "G", "G", "G-LCFA", "G-LCFA", "LCFA", "LCFA"),
      levels = c("G-D", "G", "G-LCFA", "LCFA")
    ),
    cluster_columns   = FALSE,
    cluster_rows      = cluster_rows,
    show_row_dend     = FALSE,
    row_names_side    = "left",
    column_gap        = unit(2, "mm"),
    row_labels        = row_labels,
    row_names_gp      = gpar(fontsize = row_fontsize),
    column_names_gp   = gpar(fontsize = col_fontsize),
    column_title      = title,
    column_title_side = "top",
    border_gp         = gpar(col = "black", lty = 2),
    show_row_names    = TRUE,
    show_column_names = TRUE,
    cell_fun = function(j, i, x, y, width, height, fill) {

      if (!is.null(highlight_mat) && isTRUE(highlight_mat[i, j])) {
        grid.rect(
          x, y, width, height,
          gp = gpar(col = highlight_col, fill = NA, lwd = highlight_lwd)
        )
      }

      if (show_numbers && !is.na(mat_plot[i, j])) {
        grid.text(round(mat_plot[i, j], 2), x, y, gp = gpar(fontsize = 6))
      }
    },
    heatmap_legend_param = list(
      title          = "log2FC",
      title_position = "leftcenter-rot",
      labels_gp      = gpar(fontsize = 6),
      title_gp       = gpar(fontsize = 8)
    )
  )

  if (!is.null(right_annotation)) {
    ht <- ht + right_annotation
  }

  ht
}

make_final_go_heatmap <- function(go_list,
                                  term_groups,
                                  LogFC_df_iron,
                                  feature_data,
                                  title,
                                  high_col,
                                  mid_col,
                                  highlight_col,
                                  output_file,
                                  use_abs = FALSE,
                                  add_annotation = TRUE,
                                  cap = NULL,
                                  cluster_rows = TRUE,
                                  row_fontsize = 11,
                                  col_fontsize = 8,
                                  highlight_lwd = 3,
                                  exact = TRUE) {

  term_gene_list <- lapply(term_groups, function(terms) {
    extract_go_genes(go_list, terms = terms, exact = exact)
  })

  genes <- unique(unlist(term_gene_list))

  if (length(genes) == 0) {
    warning(paste("No genes found for", title))
    return(NULL)
  }

  mat <- prepare_growth_arrest_matrix(genes, LogFC_df_iron)

  if (use_abs) {
    mat <- abs(mat)
  }

  row_ha <- NULL
  if (add_annotation) {
    row_ha <- make_term_annotation(
      genes = rownames(mat),
      term_gene_list = term_gene_list
    )
  }

  ht <- plot_gene_response_heatmap(
    mat              = mat,
    title            = title,
    high_col         = high_col,
    mid_col          = mid_col,
    row_labels       = make_row_labels(rownames(mat), feature_data),
    highlight_mat    = make_membership_matrix(rownames(mat), go_list),
    highlight_col    = highlight_col,
    highlight_lwd    = highlight_lwd,
    right_annotation = row_ha,
    row_fontsize     = row_fontsize,
    col_fontsize     = col_fontsize,
    cap              = cap,
    cluster_rows     = cluster_rows
  )

  size <- calc_ht_size(ht)

  pdf(output_file, width = size[1] + 1, height = size[2] + 1)
  if (add_annotation) {
    draw(ht, annotation_legend_list = list(make_presence_legend()))
  } else {
    draw(ht)
  }
  dev.off()

  list(
    ht = ht,
    genes = rownames(mat),
    term_gene_list = term_gene_list
  )
}


##############################
### 2. Paths & Directories  ###
##############################

input_folder <- "Inputs/002_Processed_data/GO_terms"
input_dir    <- "Inputs/002_Processed_data"
output_dir   <- "Outputs"

final_heatmap_dir <- file.path(output_dir, "Final_GO_gene_heatmaps")
dir.create(final_heatmap_dir, recursive = TRUE, showWarnings = FALSE)


####################
### 3. Load Data  ###
####################

lista_EU <- readRDS(file.path(input_folder, "lista_up_GO_enrichment_OR_4_fdr_0.05.rds"))
lista_ED <- readRDS(file.path(input_folder, "lista_down_GO_enrichment_OR_4_fdr_0.05.rds"))

feature_data <- read.table(file.path(input_dir, "txt", "feature_data_filtered.txt"))

myData   <- readRDS(file.path(input_dir, "RDS/Contrasts_stat_0.05.RDS"))
LogFC_df <- read.table(file.path(input_dir, "txt", "LogFC_0.05.txt"))
BH_df    <- read.table(file.path(input_dir, "txt", "BH_0.05.txt"))
lfcSE_df <- read.table(file.path(input_dir, "txt", "lfcSE_0.05.txt"))

vec_iron      <- 9:16
LogFC_df_iron <- LogFC_df[, vec_iron]
BH_df_iron    <- BH_df[, vec_iron]
lfcSE_df_iron <- lfcSE_df[, vec_iron]


##############################
### 4. Heatmap 1: EU Iron  ###
##############################

hm1 <- make_final_go_heatmap(
  go_list = lista_EU,
  term_groups = list(
    siderophore_metabolism = c("siderophore metabolic process"),
    siderophore_biosynthesis = c("siderophore biosynthetic process"),
    iron_starvation = c("cellular response to iron ion starvation")
  ),
  LogFC_df_iron = LogFC_df_iron,
  feature_data = feature_data,
  title = "EU siderophore / iron deprivation genes",
  high_col = "#cb181d",
  mid_col = "#fcae91",
  highlight_col = "#FFA373",
  output_file = file.path(final_heatmap_dir, "Heatmap_1_EU_siderophore_iron.pdf"),
  use_abs = FALSE,
  add_annotation = TRUE,
  highlight_lwd = 3.5,
  row_fontsize = 11
)


################################
### 5. Heatmap 2: EU Metals  ###
################################

metal_terms <- list(
  metal_homeostasis = c("metal ion homeostasis", "cellular metal ion homeostasis"),
  metal_response = c("response to metal ion", "cellular response to metal ion"),
  cadmium = c("response to cadmium ion"),
  copper = c("response to copper ion"),
  metal_transport = c("transition metal ion transport")
)

metal_term_gene_list <- lapply(metal_terms, function(terms) {
  extract_go_genes(lista_EU, terms = terms, exact = TRUE)
})

metal_genes <- unique(unlist(metal_term_gene_list))

# Remove genes already plotted in Heatmap 1
if (!is.null(hm1)) {
  metal_genes <- setdiff(metal_genes, hm1$genes)
}

# Optional: curated genes from Carmen figure
curated_metal_genes <- c(
  "ctpC", "ctpG", "ctpV", "Rv2936", "lpqS", "mymT",
  "csoR", "ricR", "cmtR", "cadI", "arsC",
  "smtB", "zur", "kmtR", "clgR",
  "bfrA", "furA", "hspX",
  "clpP1", "clpP2",
  "acr2", "Rv2640c", "Rv2642",
  "Rv0190", "Rv0846c", "Rv0968", "Rv0970",
  "Rv1993", "Rv1995",
  "pacL1"
)

# Only keep curated genes if they are already locus tags in LogFC_df_iron
curated_metal_genes <- intersect(curated_metal_genes, rownames(LogFC_df_iron))

metal_genes <- unique(c(metal_genes, curated_metal_genes))
metal_genes <- intersect(metal_genes, rownames(LogFC_df_iron))

mat_metal <- prepare_growth_arrest_matrix(metal_genes, LogFC_df_iron)

row_ha_metal <- make_term_annotation(
  genes = rownames(mat_metal),
  term_gene_list = metal_term_gene_list
)

hm2_ht <- plot_gene_response_heatmap(
  mat              = mat_metal,
  title            = "EU metal response genes",
  high_col         = "#cb181d",
  mid_col          = "#fcae91",
  row_labels       = make_row_labels(rownames(mat_metal), feature_data),
  highlight_mat    = make_membership_matrix(rownames(mat_metal), lista_EU),
  highlight_col    = "#FFA373",
  highlight_lwd    = 3.5,
  right_annotation = row_ha_metal,
  row_fontsize     = 11,
  col_fontsize     = 8,
  cluster_rows     = TRUE
)

size <- calc_ht_size(hm2_ht)

pdf(
  file.path(final_heatmap_dir, "Heatmap_2_EU_metals.pdf"),
  width = size[1] + 1,
  height = size[2] + 1
)
draw(hm2_ht, annotation_legend_list = list(make_presence_legend()))
dev.off()

hm2 <- list(
  ht = hm2_ht,
  genes = rownames(mat_metal),
  term_gene_list = metal_term_gene_list
)


################################
### 6. Heatmap 3: EU Lipids  ###
################################

hm3 <- make_final_go_heatmap(
  go_list = lista_EU,
  term_groups = list(
    cholesterol = c("cholesterol metabolic process", "cholesterol biosynthetic process"),
    acylglycerol = c("acylglycerol metabolic process", "acylglycerol biosynthetic process"),
    triglyceride = c("triglyceride metabolic process", "triglyceride biosynthetic process"),
    neutral_lipid = c("neutral lipid metabolic process", "neutral lipid biosynthetic process"),
    response_to_lipid = c("response to lipid")
  ),
  LogFC_df_iron = LogFC_df_iron,
  feature_data = feature_data,
  title = "EU lipid response genes",
  high_col = "#cb181d",
  mid_col = "#fcae91",
  highlight_col = "#FFA373",
  output_file = file.path(final_heatmap_dir, "Heatmap_3_EU_lipids.pdf"),
  use_abs = FALSE,
  add_annotation = TRUE,
  cap = 6,
  highlight_lwd = 3.5,
  row_fontsize = 11
)


###################################
### 7. Heatmap 4: ED Glycolysis ###
###################################

hm4 <- make_final_go_heatmap(
  go_list = lista_ED,
  term_groups = list(
    glycolysis = c("glycolytic process")
  ),
  LogFC_df_iron = LogFC_df_iron,
  feature_data = feature_data,
  title = "ED glycolysis genes",
  high_col = "#08519c",
  mid_col = "#9ecae1",
  highlight_col = "#50486D",
  output_file = file.path(final_heatmap_dir, "Heatmap_4_ED_glycolysis.pdf"),
  use_abs = TRUE,
  add_annotation = FALSE,
  highlight_lwd = 3.5,
  row_fontsize = 11
)


##############################
### 8. Heatmap 5: ED TCA   ###
##############################

hm5 <- make_final_go_heatmap(
  go_list = lista_ED,
  term_groups = list(
    tca_cycle = c("tricarboxylic acid cycle")
  ),
  LogFC_df_iron = LogFC_df_iron,
  feature_data = feature_data,
  title = "ED TCA cycle genes",
  high_col = "#08519c",
  mid_col = "#9ecae1",
  highlight_col = "#50486D",
  output_file = file.path(final_heatmap_dir, "Heatmap_5_ED_TCA.pdf"),
  use_abs = TRUE,
  add_annotation = FALSE,
  highlight_lwd = 3.5,
  row_fontsize = 11
)


#################################
### 9. Heatmap 6: ED OXPHOS   ###
#################################

hm6 <- make_final_go_heatmap(
  go_list = lista_ED,
  term_groups = list(
    oxphos = c("respiratory electron transport chain")
  ),
  LogFC_df_iron = LogFC_df_iron,
  feature_data = feature_data,
  title = "ED respiratory electron transport genes",
  high_col = "#08519c",
  mid_col = "#9ecae1",
  highlight_col = "#50486D",
  output_file = file.path(final_heatmap_dir, "Heatmap_6_ED_OXPHOS.pdf"),
  use_abs = TRUE,
  add_annotation = FALSE,
  highlight_lwd = 3.5,
  row_fontsize = 11
)


##########################################
### 10. Heatmap 7: ED Amino Acids      ###
##########################################

hm7 <- make_final_go_heatmap(
  go_list = lista_ED,
  term_groups = list(
    amino_acid_metabolism = c(
      "amino acid metabolic process",
      "amino acid biosynthetic process",
      "amino acid catabolic process"
    ),
    aspartate_family = c(
      "aspartate family amino acid biosynthetic process"
    ),
    threonine = c(
      "threonine metabolic process",
      "threonine biosynthetic process"
    )
  ),
  LogFC_df_iron = LogFC_df_iron,
  feature_data = feature_data,
  title = "ED amino acid metabolism genes",
  high_col = "#08519c",
  mid_col = "#9ecae1",
  highlight_col = "#50486D",
  output_file = file.path(final_heatmap_dir, "Heatmap_7_ED_amino_acids.pdf"),
  use_abs = TRUE,
  add_annotation = TRUE,
  highlight_lwd = 3.5,
  row_fontsize = 11,
  exact = TRUE
)


##############################
### 11. Save Final Objects  ###
##############################

final_heatmaps <- list(
  Heatmap_1_EU_siderophore_iron = hm1,
  Heatmap_2_EU_metals = hm2,
  Heatmap_3_EU_lipids = hm3,
  Heatmap_4_ED_glycolysis = hm4,
  Heatmap_5_ED_TCA = hm5,
  Heatmap_6_ED_OXPHOS = hm6,
  Heatmap_7_ED_amino_acids = hm7
)

saveRDS(
  final_heatmaps,
  file.path(final_heatmap_dir, "final_7_heatmaps_objects.rds")
)

final_gene_sets <- lapply(final_heatmaps, function(x) {
  if (is.null(x)) return(character(0))
  x$genes
})

saveRDS(
  final_gene_sets,
  file.path(final_heatmap_dir, "final_7_heatmaps_gene_sets.rds")
)

final_gene_sets_df <- stack(final_gene_sets)
colnames(final_gene_sets_df) <- c("gene", "heatmap")

write.table(
  final_gene_sets_df,
  file.path(final_heatmap_dir, "final_7_heatmaps_gene_sets.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


##############################
### 12. Sanity Checks       ###
##############################

cat("\nNumber of genes per final heatmap:\n")
print(lengths(final_gene_sets))

cat("\nHeatmap 3 lipid term overlaps:\n")
if (!is.null(hm3)) {
  print(sapply(hm3$term_gene_list, length))
  print(length(intersect(hm3$term_gene_list$acylglycerol,
                         hm3$term_gene_list$triglyceride)))
  print(length(intersect(hm3$term_gene_list$acylglycerol,
                         hm3$term_gene_list$neutral_lipid)))
  print(length(intersect(hm3$term_gene_list$triglyceride,
                         hm3$term_gene_list$neutral_lipid)))
}

cat("\nHeatmap 7 amino acid term overlaps:\n")
if (!is.null(hm7)) {
  print(sapply(hm7$term_gene_list, length))
  print(length(intersect(hm7$term_gene_list$amino_acid_metabolism,
                         hm7$term_gene_list$aspartate_family)))
  print(length(intersect(hm7$term_gene_list$amino_acid_metabolism,
                         hm7$term_gene_list$threonine)))
  print(length(intersect(hm7$term_gene_list$aspartate_family,
                         hm7$term_gene_list$threonine)))
}

cat("\nAvailable ED terms containing amino / threonine / respiratory / glycolytic / tricarboxylic:\n")
ed_terms <- unique(unlist(lapply(lista_ED, function(x) x$description)))
print(grep("amino|threonine|respiratory|glycolytic|tricarboxylic", ed_terms, value = TRUE, ignore.case = TRUE))
