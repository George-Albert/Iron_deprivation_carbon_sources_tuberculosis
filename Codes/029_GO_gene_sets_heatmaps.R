
######################
### 0. Dependencies ###
######################

library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(tidyverse)
library(shadowtext)

# Utility: quick column name inspection
dcols <- function(x) data.frame(colnames(x))


##########################
### 1. Helper Functions ###
##########################

# Returns rendered width and height of a ComplexHeatmap object in inches.
calc_ht_size <- function(ht, unit = "inch") {
  pdf(NULL)
  ht <- draw(ht)
  w  <- convertX(ComplexHeatmap:::width(ht),  unit, valueOnly = TRUE)
  h  <- convertY(ComplexHeatmap:::height(ht), unit, valueOnly = TRUE)
  dev.off()
  c(w, h)
}

# Builds "symbol (locus_tag)" display labels for heatmap rows.
# Falls back to bare locus_tag when no gene symbol is available.
make_row_labels <- function(genes, feature_data) {
  sym <- feature_data[genes, "symbol"]
  ifelse(
    !is.na(sym) & nchar(trimws(sym)) > 0,
    paste0(trimws(sym), " (", genes, ")"),
    genes
  )
}

# Extracts unique gene IDs from a GO enrichment list for given term descriptions.
# exact = TRUE: exact string match; FALSE: case-insensitive regex.
extract_go_genes <- function(go_list, terms, exact = TRUE) {
  genes <- c()
  for (nm in names(go_list)) {
    df     <- go_list[[nm]]
    df_sub <- if (exact) {
      df[df$description %in% terms, , drop = FALSE]
    } else {
      df[grepl(paste(terms, collapse = "|"), df$description, ignore.case = TRUE), , drop = FALSE]
    }
    if (nrow(df_sub) > 0)
      genes <- c(genes, trimws(unlist(strsplit(df_sub$genes_in_query, ","))))
  }
  unique(genes)
}

# Subsets LogFC_df_iron for the given genes and reorders columns to pair
# Fe+ / Fe- conditions side by side for each carbon source: G-D, G, G-LCFA, LCFA.
prepare_growth_arrest_matrix <- function(genes, LogFC_df_iron) {
  genes <- intersect(genes, rownames(LogFC_df_iron))
  if (length(genes) == 0) stop("None of the selected genes are present in LogFC_df_iron.")

  mat <- as.matrix(LogFC_df_iron[genes, , drop = FALSE])
  # Columns 1–4 = Fe+, 5–8 = Fe-; reorder to interleave Fe+ / Fe- per carbon source
  mat <- mat[, c(1, 5, 2, 6, 3, 7, 4, 8), drop = FALSE]
  colnames(mat) <- c(
    "G-D_Fe+", "G-D_Fe-", "G_Fe+",    "G_Fe-",
    "G-LCFA_Fe+", "G-LCFA_Fe-", "LCFA_Fe+", "LCFA_Fe-"
  )
  mat
}

# Builds a ComplexHeatmap for a gene x condition log2FC matrix.
# Row names on the left, no row dendrogram, columns grouped by carbon source.
# row_labels: optional character vector; defaults to rownames(mat).
plot_gene_response_heatmap <- function(mat,
                                       title        = "EU genes",
                                       high_col     = "#cb181d",
                                       low_col      = "white",
                                       cluster_rows = TRUE,
                                       show_numbers = FALSE,
                                       row_fontsize = 6,
                                       col_fontsize = 7,
                                       row_labels   = NULL,
                                       cap          = NULL,
                                       eu_mat       = NULL) {
  mat_plot <- mat
  mat_plot[!is.finite(mat_plot)] <- NA
  max_abs <- max(abs(mat_plot), na.rm = TRUE)

  # If cap is provided, clamp the color scale at that value.
  # Cells above the cap are shown in high_col but the scale is not stretched.
  if (!is.null(cap)) max_abs <- min(max_abs, cap)

  col_fun <- circlize::colorRamp2(
    c(0, max_abs / 2, max_abs),
    c(low_col, "#fcae91", high_col)
  )

  if (is.null(row_labels)) row_labels <- rownames(mat_plot)

  ComplexHeatmap::Heatmap(
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
      # Orange border on cells where the gene is EU in that contrast
      if (!is.null(eu_mat) && isTRUE(eu_mat[i, j]))
        grid.rect(x, y, width, height, gp = gpar(col = "#FFA373", fill = NA, lwd = 1.5))
      # Optionally overlay numeric log2FC values
      if (show_numbers && !is.na(mat_plot[i, j]))
        grid.text(round(mat_plot[i, j], 2), x, y, gp = gpar(fontsize = 6))
    },
    heatmap_legend_param = list(
      title          = "log2FC",
      title_position = "leftcenter-rot",
      labels_gp      = gpar(fontsize = 6),
      title_gp       = gpar(fontsize = 8)
    )
  )
}

# Builds a boolean matrix marking EU cells (TRUE = gene is EU in that contrast).
# Columns are paired Fe+/Fe- per carbon source, matching prepare_growth_arrest_matrix().
# lista_EU[[1]] = G-D, [[2]] = G, [[3]] = G-LCFA; LCFA has no EU list → always FALSE.
# Depends on global: lista_EU.
make_eu_matrix <- function(genes) {
  eu_per_contrast <- lapply(lista_EU, function(df) {
    unique(trimws(unlist(strsplit(df$genes_in_query, ","))))
  })

  col_names <- c("G-D_Fe+", "G-D_Fe-", "G_Fe+",    "G_Fe-",
                 "G-LCFA_Fe+", "G-LCFA_Fe-", "LCFA_Fe+", "LCFA_Fe-")

  eu_mat <- matrix(FALSE, nrow = length(genes), ncol = 8,
                   dimnames = list(genes, col_names))

  eu_mat[genes %in% eu_per_contrast[[1]], c(1, 2)] <- TRUE  # G-D
  eu_mat[genes %in% eu_per_contrast[[2]], c(3, 4)] <- TRUE  # G
  eu_mat[genes %in% eu_per_contrast[[3]], c(5, 6)] <- TRUE  # G-LCFA

  eu_mat
}

# Resolves gene symbols or locus tags to rownames of feature_data.
# Returns only successfully resolved, unique locus tags.
# Depends on globals: feature_data, sym_to_locus (built in Section 3).
resolve_ids <- function(ids) {
  loci <- sapply(ids, function(id) {
    if      (id %in% rownames(feature_data))   id
    else if (id %in% names(sym_to_locus))      unname(sym_to_locus[id])
    else                                       NA_character_
  }, USE.NAMES = FALSE)
  unique(loci[!is.na(loci)])
}

# Converts a logical vector to "Presence" / "Absence" character vector.
make_anno_vec <- function(x) ifelse(x, "Presence", "Absence")


##############################
### 2. Paths & Directories  ###
##############################

input_folder    <- "Inputs/002_Processed_data/GO_terms"
input_dir       <- "Inputs/002_Processed_data"
output_dir      <- "Outputs"
out_heatmap_dir <- file.path(output_dir, "EU_GO_gene_heatmaps")
dir.create(out_heatmap_dir, recursive = TRUE, showWarnings = FALSE)


####################
### 3. Load Data  ###
####################

# GO enrichment results: enhanced (EU) and damped (ED) gene lists per contrast
lista_EU <- readRDS(file.path(input_folder, "lista_up_GO_enrichment_OR_4_fdr_0.05.rds"))
lista_ED <- readRDS(file.path(input_folder, "lista_down_GO_enrichment_OR_4_fdr_0.05.rds"))

# Gene annotation table (locus tags, gene symbols, functional info)
feature_data <- read.table(file.path(input_dir, "txt", "feature_data_filtered.txt"))

# Symbol-to-locus lookup derived from feature_data (used by resolve_ids)
fd_sym       <- feature_data[nchar(trimws(coalesce(feature_data$symbol, ""))) > 0, ]
sym_to_locus <- setNames(rownames(fd_sym), trimws(fd_sym$symbol))

# Differential expression results for all contrasts (FDR < 0.05)
myData   <- readRDS(file.path(input_dir, "RDS/Contrasts_stat_0.05.RDS"))
LogFC_df <- read.table(file.path(input_dir, "txt/LogFC_0.05.txt"))
BH_df    <- read.table(file.path(input_dir, "txt/BH_0.05.txt"))
lfcSE_df <- read.table(file.path(input_dir, "txt/lfcSE_0.05.txt"))

# Columns 9–16 correspond to the 8 iron-deprivation contrasts
vec_iron      <- 9:16
LogFC_df_iron <- LogFC_df[, vec_iron]
BH_df_iron    <- BH_df[,    vec_iron]
lfcSE_df_iron <- lfcSE_df[, vec_iron]


##########################
### 4. Define Gene Sets ###
##########################

### 4a. EU GO term sets -------------------------------------------------------
# Each entry maps a set name to one or more GO term descriptions.
# Genes are extracted from lista_EU in Section 5.
eu_sets <- list(

  set1_siderophore = c(
    "siderophore metabolic process",
    "siderophore biosynthetic process"
  ),

  set2_iron_starvation = c(
    "cellular response to iron ion starvation"
  ),

  set3_siderophore_iron_starvation = c(
    "siderophore metabolic process",
    "siderophore biosynthetic process",
    "cellular response to iron ion starvation"
  ),

  set4_metal_homeostasis = c(
    "metal ion homeostasis",
    "cellular metal ion homeostasis"
  ),

  set5_metal_response = c(
    "response to metal ion",
    "cellular response to metal ion"
  ),

  set6_cadmium_response = c(
    "response to cadmium ion"
  ),

  set7_copper_response = c(
    "response to copper ion"
  ),

  set8_transition_metal_transport = c(
    "transition metal ion transport"
  ),

  # set9 pending definition
  # eu_sets$set9_lipid_top_terms <- get_terms_by_position(lista_EU, c(1, 2, 4, 9))

  set10_cholesterol = c(
    "cholesterol metabolic process",
    "cholesterol biosynthetic process"
  ),

  set11_acylglycerol = c(
    "acylglycerol metabolic process",
    "acylglycerol biosynthetic process"
  ),

  set12_triglyceride = c(
    "triglyceride metabolic process",
    "triglyceride biosynthetic process"
  ),

  set13_neutral_lipid = c(
    "neutral lipid metabolic process",
    "neutral lipid biosynthetic process"
  ),

  set14_response_to_lipids = c(
    "response to lipid"
  )
)

### 4b. Curated metal gene lists ----------------------------------------------
# Gene IDs can be symbols or locus tags; resolved to locus tags via resolve_ids.
# Categories: Cu (copper), Zn (zinc), Cd (cadmium), As (arsenic),
#             Metal (general metal-response genes not specific to one metal).
metal_gene_list <- list(
  Cu    = c("Rv0190", "Rv2963", "Rv0846c", "lpqS", "mymT",
            "csoR", "Rv0968", "ctpV", "Rv0970", "ctpG"),
  Zn    = c("smtB", "zur"),
  Cd    = c("cadI", "Rv2642", "arsC", "cmtR"),
  As    = c("arsC"),
  Metal = c("bfrA", "furA", "hspX", "kmtR")
)

metal_gene_list_loci <- lapply(metal_gene_list, resolve_ids)


#####################################################
### 5. Generate Individual EU Set Heatmaps        ###
#####################################################
# For each GO set: extract genes from lista_EU, build the log2FC matrix,
# plot and export the heatmap. Gene sets and heatmaps are stored in lists
# for downstream use.

eu_heatmaps  <- list()
eu_gene_sets <- list()

for (set_name in names(eu_sets)) {

  genes_i <- extract_go_genes(lista_EU, terms = eu_sets[[set_name]], exact = TRUE)
  eu_gene_sets[[set_name]] <- genes_i

  if (length(genes_i) == 0) {
    warning(paste("No genes found for", set_name))
    next
  }

  mat_i <- prepare_growth_arrest_matrix(genes_i, LogFC_df_iron)
  ht_i  <- plot_gene_response_heatmap(
    mat          = mat_i,
    title        = set_name,
    high_col     = "#cb181d",
    cluster_rows = TRUE,
    row_labels   = make_row_labels(rownames(mat_i), feature_data),
    eu_mat       = make_eu_matrix(rownames(mat_i))
  )
  eu_heatmaps[[set_name]] <- ht_i

  size <- calc_ht_size(ht_i)
  pdf(file.path(out_heatmap_dir, paste0(set_name, ".pdf")),
      width = size[1] + 0.7, height = size[2] + 0.7)
  draw(ht_i)
  dev.off()
}


##################################################################
### 6. Regenerate Lipid Sets with Capped Color Scale         ###
##################################################################
# Sets 11–13 share a single extreme outlier (~10.9) that compresses
# the rest of the color scale. Capping at 6 (≈ 95th percentile)
# preserves the outlier in max color while improving differentiation.

capped_sets <- c("set11_acylglycerol", "set12_triglyceride", "set13_neutral_lipid")

for (s in capped_sets) {
  mat_i <- prepare_growth_arrest_matrix(eu_gene_sets[[s]], LogFC_df_iron)
  ht_i  <- plot_gene_response_heatmap(
    mat          = mat_i,
    title        = s,
    high_col     = "#cb181d",
    cluster_rows = TRUE,
    row_labels   = make_row_labels(rownames(mat_i), feature_data),
    cap          = 6,
    eu_mat       = make_eu_matrix(rownames(mat_i))
  )
  eu_heatmaps[[s]] <- ht_i
  size <- calc_ht_size(ht_i)
  pdf(file.path(out_heatmap_dir, paste0(s, ".pdf")),
      width = size[1] + 0.7, height = size[2] + 0.7)
  draw(ht_i)
  dev.off()
}


##############################################################
### 8. Augment Metal-Related Sets with Curated Genes       ###
##############################################################
# The main loop builds eu_gene_sets from GO enrichment results only.
# Here we add curated genes from metal_gene_list_loci to each biologically
# relevant set and regenerate those heatmap PDFs.

metal_set_mapping <- list(
  set4_metal_homeostasis          = unique(unlist(metal_gene_list_loci)),
  set5_metal_response             = unique(unlist(metal_gene_list_loci)),
  set6_cadmium_response           = unique(c(metal_gene_list_loci$Cd, metal_gene_list_loci$As)),
  set7_copper_response            = metal_gene_list_loci$Cu,
  set8_transition_metal_transport = unique(c(metal_gene_list_loci$Cu,
                                             metal_gene_list_loci$Zn,
                                             metal_gene_list_loci$Cd))
)

for (s in names(metal_set_mapping)) {
  # Only add genes with LogFC data not already in the set
  new_genes <- intersect(
    setdiff(unique(metal_set_mapping[[s]]), eu_gene_sets[[s]]),
    rownames(LogFC_df_iron)
  )
  if (length(new_genes) == 0) next

  eu_gene_sets[[s]] <- unique(c(eu_gene_sets[[s]], new_genes))

  mat_i <- prepare_growth_arrest_matrix(eu_gene_sets[[s]], LogFC_df_iron)
  ht_i  <- plot_gene_response_heatmap(
    mat          = mat_i,
    title        = s,
    high_col     = "#cb181d",
    cluster_rows = TRUE,
    row_labels   = make_row_labels(rownames(mat_i), feature_data),
    eu_mat       = make_eu_matrix(rownames(mat_i))
  )
  eu_heatmaps[[s]] <- ht_i

  size <- calc_ht_size(ht_i)
  pdf(file.path(out_heatmap_dir, paste0(s, ".pdf")),
      width = size[1] + 0.7, height = size[2] + 0.7)
  draw(ht_i)
  dev.off()
}


##################################################################
### 9. Combined Metal Response Heatmap with Metal Annotation   ###
##################################################################
# Builds a single heatmap from the union of GO-derived metal genes
# (sets 5, 6, 7) and all curated metal_gene_list genes.
# A right-side binary annotation indicates Cu/Zn/Cd/As/Metal membership.

metal_genes     <- unique(c(eu_gene_sets$set5_metal_response,
                             eu_gene_sets$set6_cadmium_response,
                             eu_gene_sets$set7_copper_response))
all_metal_genes <- unique(c(metal_genes, unlist(metal_gene_list_loci)))
mat_metal       <- prepare_growth_arrest_matrix(all_metal_genes, LogFC_df_iron)

# Binary annotation: Cu/Zn/Cd/As from curated metal_gene_list_loci;
# Metal column from GO terms "response to metal ion" / "cellular response to metal ion"
# directly in lista_EU — a gene can be TRUE in more than one column.
genes_response_to_metal <- extract_go_genes(
  lista_EU,
  terms = c("response to metal ion", "cellular response to metal ion"),
  exact = TRUE
)

metal_annotation_df <- as.data.frame(
  lapply(metal_gene_list_loci[c("Cu","Zn","Cd","As")],
         function(loci) rownames(mat_metal) %in% loci),
  row.names = rownames(mat_metal)
)
metal_annotation_df$Metal <- rownames(mat_metal) %in% genes_response_to_metal

# Row annotation with per-cell black borders; auto-legends suppressed
anno_col_map <- c("Presence" = "black", "Absence" = "white")

row_ha_metal <- rowAnnotation(
  Cu    = anno_simple(make_anno_vec(metal_annotation_df$Cu),    col = anno_col_map, gp = gpar(col = "black")),
  Zn    = anno_simple(make_anno_vec(metal_annotation_df$Zn),    col = anno_col_map, gp = gpar(col = "black")),
  Cd    = anno_simple(make_anno_vec(metal_annotation_df$Cd),    col = anno_col_map, gp = gpar(col = "black")),
  As    = anno_simple(make_anno_vec(metal_annotation_df$As),    col = anno_col_map, gp = gpar(col = "black")),
  Metal = anno_simple(make_anno_vec(metal_annotation_df$Metal), col = anno_col_map, gp = gpar(col = "black")),
  annotation_name_gp = gpar(fontsize = 6),
  border      = TRUE,
  show_legend = rep(FALSE, 5)
)

# Custom legend: explicit grid.rect() to ensure black border is visible on white tile
lgd_metal <- Legend(
  labels   = c("Absence", "Presence"),
  title    = "",
  graphics = list(
    function(x, y, w, h) grid.rect(x, y, w, h, gp = gpar(fill = "white", col = "black", lwd = 1)),
    function(x, y, w, h) grid.rect(x, y, w, h, gp = gpar(fill = "black", col = "black", lwd = 1))
  )
)

# Annotation placed to the RIGHT of the heatmap using +
ht_metal <- plot_gene_response_heatmap(
  mat          = mat_metal,
  title        = "Metal response EU genes",
  high_col     = "#cb181d",
  cluster_rows = TRUE,
  row_labels   = make_row_labels(rownames(mat_metal), feature_data),
  eu_mat       = make_eu_matrix(rownames(mat_metal))
) + row_ha_metal

pdf(file.path(out_heatmap_dir, "combined_metal_response_cadmium_copper.pdf"),
    width = 7, height = 6)
draw(ht_metal, annotation_legend_list = list(lgd_metal))
dev.off()


##############################
### 10. Export Gene Sets    ###
##############################

# Save final eu_gene_sets (after curation patches) as RDS and flat text
saveRDS(eu_gene_sets,
        file.path(out_heatmap_dir, "EU_GO_gene_sets_used_for_heatmaps.rds"))

eu_gene_sets_df <- stack(eu_gene_sets)
colnames(eu_gene_sets_df) <- c("gene", "set")
write.table(eu_gene_sets_df,
            file.path(out_heatmap_dir, "EU_GO_gene_sets_used_for_heatmaps.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
