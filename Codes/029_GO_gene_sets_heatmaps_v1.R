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
    paste0(trimws(sym), " (", genes, ")"),
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
  # Reorder as carbon-paired Fe+/Fe-
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

  # lista[[1]] = G-D
  # lista[[2]] = G
  # lista[[3]] = G-LCFA
  # LCFA has no enrichment list here
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

plot_gene_response_heatmap <- function(mat,
                                       title         = "genes",
                                       high_col      = "#cb181d",
                                       mid_col       = "#fcae91",
                                       low_col       = "white",
                                       cluster_rows  = TRUE,
                                       show_numbers  = FALSE,
                                       row_fontsize  = 6,
                                       col_fontsize  = 7,
                                       row_labels    = NULL,
                                       cap           = NULL,
                                       highlight_mat = NULL,
                                       highlight_col = "#FFA373") {

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

  Heatmap(
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
          gp = gpar(col = highlight_col, fill = NA, lwd = 1.5)
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
}

generate_go_heatmaps <- function(go_list,
                                 go_sets,
                                 LogFC_df_iron,
                                 feature_data,
                                 output_dir,
                                 prefix,
                                 high_col,
                                 mid_col,
                                 highlight_col,
                                 exact_default = TRUE,
                                 exact_overrides = list(),
                                 cap_sets = NULL,
                                 cap_value = NULL) {

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  heatmaps  <- list()
  gene_sets <- list()

  for (set_name in names(go_sets)) {

    exact_i <- exact_default
    if (set_name %in% names(exact_overrides)) {
      exact_i <- exact_overrides[[set_name]]
    }

    genes_i <- extract_go_genes(
      go_list = go_list,
      terms   = go_sets[[set_name]],
      exact   = exact_i
    )

    genes_i <- unique(genes_i)
    gene_sets[[set_name]] <- genes_i

    if (length(genes_i) == 0) {
      warning(paste("No genes found for", prefix, set_name))
      next
    }

    mat_i <- prepare_growth_arrest_matrix(
      genes        = genes_i,
      LogFC_df_iron = LogFC_df_iron
    )

    cap_i <- NULL
    if (!is.null(cap_sets) && set_name %in% cap_sets) {
      cap_i <- cap_value
    }

    ht_i <- plot_gene_response_heatmap(
      mat           = mat_i,
      title         = paste(prefix, set_name),
      high_col      = high_col,
      mid_col       = mid_col,
      low_col       = "white",
      cluster_rows  = TRUE,
      row_labels    = make_row_labels(rownames(mat_i), feature_data),
      cap           = cap_i,
      highlight_mat = make_membership_matrix(rownames(mat_i), go_list),
      highlight_col = highlight_col
    )

    heatmaps[[set_name]] <- ht_i

    size <- calc_ht_size(ht_i)

    pdf(
      file = file.path(output_dir, paste0(prefix, "_", set_name, ".pdf")),
      width = size[1] + 0.7,
      height = size[2] + 0.7
    )
    draw(ht_i)
    dev.off()
  }

  saveRDS(
    gene_sets,
    file.path(output_dir, paste0(prefix, "_GO_gene_sets_used_for_heatmaps.rds"))
  )

  gene_sets_df <- stack(gene_sets)
  colnames(gene_sets_df) <- c("gene", "set")

  write.table(
    gene_sets_df,
    file.path(output_dir, paste0(prefix, "_GO_gene_sets_used_for_heatmaps.txt")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  list(
    heatmaps  = heatmaps,
    gene_sets = gene_sets
  )
}

resolve_ids <- function(ids) {
  loci <- sapply(ids, function(id) {
    if (id %in% rownames(feature_data)) {
      id
    } else if (id %in% names(sym_to_locus)) {
      unname(sym_to_locus[id])
    } else {
      NA_character_
    }
  }, USE.NAMES = FALSE)

  unique(loci[!is.na(loci)])
}

make_anno_vec <- function(x) ifelse(x, "Presence", "Absence")


##############################
### 2. Paths & Directories  ###
##############################

input_folder <- "Inputs/002_Processed_data/GO_terms"
input_dir    <- "Inputs/002_Processed_data"
output_dir   <- "Outputs"

out_eu_heatmap_dir <- file.path(output_dir, "EU_GO_gene_heatmaps")
out_ed_heatmap_dir <- file.path(output_dir, "ED_GO_gene_heatmaps")

dir.create(out_eu_heatmap_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_ed_heatmap_dir, recursive = TRUE, showWarnings = FALSE)


####################
### 3. Load Data  ###
####################

lista_EU <- readRDS(file.path(input_folder, "lista_up_GO_enrichment_OR_4_fdr_0.05.rds"))
lista_ED <- readRDS(file.path(input_folder, "lista_down_GO_enrichment_OR_4_fdr_0.05.rds"))

feature_data <- read.table(file.path(input_dir, "txt", "feature_data_filtered.txt"))

fd_sym <- feature_data[nchar(trimws(coalesce(feature_data$symbol, ""))) > 0, ]
sym_to_locus <- setNames(rownames(fd_sym), trimws(fd_sym$symbol))

myData   <- readRDS(file.path(input_dir, "RDS/Contrasts_stat_0.05.RDS"))
LogFC_df <- read.table(file.path(input_dir, "txt/LogFC_0.05.txt"))
BH_df    <- read.table(file.path(input_dir, "txt/BH_0.05.txt"))
lfcSE_df <- read.table(file.path(input_dir, "txt/lfcSE_0.05.txt"))

vec_iron      <- 9:16
LogFC_df_iron <- LogFC_df[, vec_iron]
BH_df_iron    <- BH_df[, vec_iron]
lfcSE_df_iron <- lfcSE_df[, vec_iron]


##########################
### 4. Define EU Sets  ###
##########################

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


##########################
### 5. Define ED Sets  ###
##########################

ed_sets <- list(

  set1_glycolysis = c(
    "glycolytic process"
  ),

  set2_tca_cycle = c(
    "tricarboxylic acid cycle"
  ),

  set3_respiratory_etc = c(
    "respiratory electron transport chain"
  ),

  set4_amino_acid_metabolism = c(
    "amino acid metabolic process",
    "amino acid biosynthetic process",
    "amino acid catabolic process"
  ),

  set5_aspartate_family = c(
    "aspartate family amino acid biosynthetic process"
  ),

  set6_threonine = c(
    "threonine metabolic process",
    "threonine biosynthetic process"
  )
)

ed_exact_overrides <- list(
  set4_amino_acid_metabolism = FALSE
)


############################################
### 6. Generate Individual EU Heatmaps   ###
############################################

res_eu <- generate_go_heatmaps(
  go_list       = lista_EU,
  go_sets       = eu_sets,
  LogFC_df_iron = LogFC_df_iron,
  feature_data  = feature_data,
  output_dir    = out_eu_heatmap_dir,
  prefix        = "EU",
  high_col      = "#cb181d",
  mid_col       = "#fcae91",
  highlight_col = "#FFA373",
  exact_default = TRUE,
  cap_sets      = c("set11_acylglycerol", "set12_triglyceride", "set13_neutral_lipid"),
  cap_value     = 6
)

eu_heatmaps  <- res_eu$heatmaps
eu_gene_sets <- res_eu$gene_sets


############################################
### 7. Generate Individual ED Heatmaps   ###
############################################

res_ed <- generate_go_heatmaps(
  go_list         = lista_ED,
  go_sets         = ed_sets,
  LogFC_df_iron   = abs(LogFC_df_iron),
  feature_data    = feature_data,
  output_dir      = out_ed_heatmap_dir,
  prefix          = "ED",
  high_col        = "#08519c",
  mid_col         = "#9ecae1",
  highlight_col   = "#50486D",
  exact_default   = TRUE,
  exact_overrides = ed_exact_overrides
)

ed_heatmaps  <- res_ed$heatmaps
ed_gene_sets <- res_ed$gene_sets


##############################################################
### 8. Curated Metal Gene Lists for EU Metal Heatmaps      ###
##############################################################

metal_gene_list <- list(
  Cu    = c("Rv0190", "Rv2963", "Rv0846c", "lpqS", "mymT",
            "csoR", "Rv0968", "ctpV", "Rv0970", "ctpG"),
  Zn    = c("smtB", "zur"),
  Cd    = c("cadI", "Rv2642", "arsC", "cmtR"),
  As    = c("arsC"),
  Metal = c("bfrA", "furA", "hspX", "kmtR")
)

metal_gene_list_loci <- lapply(metal_gene_list, resolve_ids)

metal_set_mapping <- list(
  set4_metal_homeostasis          = unique(unlist(metal_gene_list_loci)),
  set5_metal_response             = unique(unlist(metal_gene_list_loci)),
  set6_cadmium_response           = unique(c(metal_gene_list_loci$Cd, metal_gene_list_loci$As)),
  set7_copper_response            = metal_gene_list_loci$Cu,
  set8_transition_metal_transport = unique(c(
    metal_gene_list_loci$Cu,
    metal_gene_list_loci$Zn,
    metal_gene_list_loci$Cd
  ))
)

for (s in names(metal_set_mapping)) {

  if (!s %in% names(eu_gene_sets)) next

  new_genes <- intersect(
    setdiff(unique(metal_set_mapping[[s]]), eu_gene_sets[[s]]),
    rownames(LogFC_df_iron)
  )

  if (length(new_genes) == 0) next

  eu_gene_sets[[s]] <- unique(c(eu_gene_sets[[s]], new_genes))

  mat_i <- prepare_growth_arrest_matrix(eu_gene_sets[[s]], LogFC_df_iron)

  ht_i <- plot_gene_response_heatmap(
    mat           = mat_i,
    title         = paste("EU", s),
    high_col      = "#cb181d",
    mid_col       = "#fcae91",
    cluster_rows  = TRUE,
    row_labels    = make_row_labels(rownames(mat_i), feature_data),
    highlight_mat = make_membership_matrix(rownames(mat_i), lista_EU),
    highlight_col = "#FFA373"
  )

  eu_heatmaps[[s]] <- ht_i

  size <- calc_ht_size(ht_i)

  pdf(
    file.path(out_eu_heatmap_dir, paste0("EU_", s, ".pdf")),
    width = size[1] + 0.7,
    height = size[2] + 0.7
  )
  draw(ht_i)
  dev.off()
}


##################################################################
### 9. Combined Metal Response Heatmap with Annotation         ###
##################################################################

metal_genes <- unique(c(
  eu_gene_sets$set5_metal_response,
  eu_gene_sets$set6_cadmium_response,
  eu_gene_sets$set7_copper_response
))

all_metal_genes <- unique(c(metal_genes, unlist(metal_gene_list_loci)))
all_metal_genes <- intersect(all_metal_genes, rownames(LogFC_df_iron))

mat_metal <- prepare_growth_arrest_matrix(all_metal_genes, LogFC_df_iron)

genes_response_to_metal <- extract_go_genes(
  lista_EU,
  terms = c("response to metal ion", "cellular response to metal ion"),
  exact = TRUE
)

metal_annotation_df <- as.data.frame(
  lapply(metal_gene_list_loci[c("Cu", "Zn", "Cd", "As")], function(loci) {
    rownames(mat_metal) %in% loci
  }),
  row.names = rownames(mat_metal)
)

metal_annotation_df$Metal <- rownames(mat_metal) %in% genes_response_to_metal

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

lgd_metal <- Legend(
  labels = c("Absence", "Presence"),
  title  = "",
  graphics = list(
    function(x, y, w, h) grid.rect(x, y, w, h, gp = gpar(fill = "white", col = "black", lwd = 1)),
    function(x, y, w, h) grid.rect(x, y, w, h, gp = gpar(fill = "black", col = "black", lwd = 1))
  )
)

ht_metal <- plot_gene_response_heatmap(
  mat           = mat_metal,
  title         = "EU metal response genes",
  high_col      = "#cb181d",
  mid_col       = "#fcae91",
  cluster_rows  = TRUE,
  row_labels    = make_row_labels(rownames(mat_metal), feature_data),
  highlight_mat = make_membership_matrix(rownames(mat_metal), lista_EU),
  highlight_col = "#FFA373"
) + row_ha_metal

pdf(
  file.path(out_eu_heatmap_dir, "EU_combined_metal_response_with_annotation.pdf"),
  width = 7,
  height = 6
)
draw(ht_metal, annotation_legend_list = list(lgd_metal))
dev.off()


##############################
### 10. Export Final Sets  ###
##############################

saveRDS(
  eu_gene_sets,
  file.path(out_eu_heatmap_dir, "EU_GO_gene_sets_used_for_heatmaps_final.rds")
)

eu_gene_sets_df <- stack(eu_gene_sets)
colnames(eu_gene_sets_df) <- c("gene", "set")

write.table(
  eu_gene_sets_df,
  file.path(out_eu_heatmap_dir, "EU_GO_gene_sets_used_for_heatmaps_final.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  ed_gene_sets,
  file.path(out_ed_heatmap_dir, "ED_GO_gene_sets_used_for_heatmaps_final.rds")
)

ed_gene_sets_df <- stack(ed_gene_sets)
colnames(ed_gene_sets_df) <- c("gene", "set")

write.table(
  ed_gene_sets_df,
  file.path(out_ed_heatmap_dir, "ED_GO_gene_sets_used_for_heatmaps_final.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


##############################
### 11. Optional Checks     ###
##############################

cat("\nEU sets generated:\n")
print(lengths(eu_gene_sets))

cat("\nED sets generated:\n")
print(lengths(ed_gene_sets))

cat("\nAvailable ED terms containing amino / threonine / respiratory:\n")
ed_terms <- unique(unlist(lapply(lista_ED, function(x) x$description)))
print(grep("amino|threonine|respiratory|glycolytic|tricarboxylic", ed_terms, value = TRUE, ignore.case = TRUE))
