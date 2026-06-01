
{
  library(tidyverse)
  library(ggplot2)
  library(reshape2)
  library(dplyr)
  library(gridExtra)
  library(ComplexHeatmap)
  library(dendextend)
  library(umap)
  library(rgl)
  library(RColorBrewer)
  library(circlize)
  library(eulerr)
  library(ggrepel)
  library(cowplot)
  library(grid)
  library(ggthemes)
  library(openxlsx)
  library(svglite)


}

############################
### 1. Declare functions ###
############################
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

build_top5_df <- function(ora_results,
                          gene_set = "EU",
                          top_n = 5,
                          sort_by = c("fdr","log2OR"),
                          value_col = "log2OR") {

  sort_by <- match.arg(sort_by)
  conditions <- names(ora_results)

  # 1) top regulones por condición
  top_regs <- lapply(conditions, function(cond) {
    if (!gene_set %in% names(ora_results[[cond]])) return(character(0))

    res <- ora_results[[cond]][[gene_set]]$all
    if (is.null(res) || nrow(res) == 0) return(character(0))

    res <- res[!is.na(res[[value_col]]) & is.finite(res[[value_col]]), ]
    if (nrow(res) == 0) return(character(0))

    if (sort_by == "fdr") {
      res <- res[order(res$fdr), ]
    } else {
      res <- res[order(abs(res[[value_col]]), decreasing = TRUE), ]
    }

    head(res$Regulator, top_n)
  })

  # 2) regulones únicos (≤ 20)
  regulons <- unique(unlist(top_regs))

  # 3) construir df vacío
  df <- data.frame(Regulon = regulons, stringsAsFactors = FALSE)

  # 4) rellenar columnas por condición
  for (cond in conditions) {
    if (!gene_set %in% names(ora_results[[cond]])) {
      df[[cond]] <- NA_real_
      next
    }

    res <- ora_results[[cond]][[gene_set]]$all
    v <- res[[value_col]]
    names(v) <- res$Regulator

    df[[cond]] <- v[df$Regulon]
  }

  return(df)
}
# Extrae una matriz (regulones x condiciones) para una métrica dada de $all
extract_metric_matrix <- function(ora_results,
                                  gene_set = "EU",
                                  value_col = "log2OR",
                                  key_col = c("Regulator","Description")) {

  key_col <- match.arg(key_col)

  conditions <- names(ora_results)

  # Lista de data.frames por condición (solo los que existan)
  df_list <- lapply(conditions, function(cond) {
    if (is.null(ora_results[[cond]]) || !gene_set %in% names(ora_results[[cond]])) return(NULL)

    res <- ora_results[[cond]][[gene_set]]$all
    if (is.null(res) || nrow(res) == 0) return(NULL)

    # asegúrate de que exista la col deseada
    if (!value_col %in% colnames(res)) return(NULL)

    out <- res %>%
      transmute(
        Regulon = .data[[key_col]],
        value   = .data[[value_col]]
      )
    out$condition <- cond
    out
  })
  df_list <- Filter(Negate(is.null), df_list)
  if (length(df_list) == 0) return(NULL)

  long <- bind_rows(df_list)

  # wide: filas regulón, columnas condición
  wide <- long %>%
    tidyr::pivot_wider(names_from = condition, values_from = value)

  # a matrix
  mat <- as.matrix(wide[, -1, drop = FALSE])
  rownames(mat) <- wide$Regulon
  mat
}

# Construye un set de regulones = top_n por condición (según fdr o según |log2OR|)
get_top_regulons_union <- function(ora_results,
                                   gene_set = "EU",
                                   top_n = 5,
                                   sort_by = c("fdr","abs_log2OR")) {

  sort_by <- match.arg(sort_by)
  conditions <- names(ora_results)

  regs <- unlist(lapply(conditions, function(cond) {
    if (is.null(ora_results[[cond]]) || !gene_set %in% names(ora_results[[cond]])) return(character(0))
    res <- ora_results[[cond]][[gene_set]]$all
    if (is.null(res) || nrow(res) == 0) return(character(0))

    res <- res %>% filter(!is.na(log2OR) & is.finite(log2OR))

    if (nrow(res) == 0) return(character(0))

    if (sort_by == "fdr") {
      res <- res[order(res$fdr), ]
    } else {
      res <- res[order(abs(res$log2OR), decreasing = TRUE), ]
    }

    head(res$Regulator, top_n)
  }))

  unique(regs)
}

###################################################
### 2. Set working directory and create folders ###
###################################################
main_wd <- getwd()
input_dir <- "Inputs/002_Processed_data"
output_dir <- "Outputs"

#################
##  Load Data  ##
#################
# myData <- readRDS(file.path(input_dir,"RDS/Contrasts_stat_0.05.RDS"))

LogFC_df <- read.table(file.path(input_dir,"txt/LogFC_0.05.txt"))
BH_df    <- read.table(file.path(input_dir,"txt/BH_0.05.txt"))
lfcSE_df <- read.table(file.path(input_dir,"txt/lfcSE_0.05.txt"))
feature_data <- read.table(file.path(input_dir,"txt/feature_data_filtered.txt"))
# # Define the pattern to remove the parenthesis
# names(myData)
# pattern <- paste(c("[(]", "[)]"), collapse = "|")

### Extract the 4 interactions
# 17 Effect_of_Fe_in_response_to_Growth_arrest_G_D.log2FoldChange_shrunken
# 18   Effect_of_Fe_in_response_to_Growth_arrest_G.log2FoldChange_shrunken
# 19 Effect_of_Fe_in_response_to_Growth_arrest_G_L.log2FoldChange_shrunken
# 20   Effect_of_Fe_in_response_to_Growth_arrest_L.log2FoldChange_shrunken

vec_int <- c(17:20)

LogFC_df_int <- LogFC_df[,vec_int]
BH_df_int <- BH_df[,vec_int]
lfcSE_df_int <- lfcSE_df[,vec_int]

colnames(LogFC_df_int) <- paste0(substring(colnames(LogFC_df_int),1,nchar(colnames(LogFC_df_int))-24),".LogFC")
# colnames(BH_df_int) <- paste0(colnames(int_BH),".BH")
colnames(lfcSE_df_int) <- paste0(substring(colnames(lfcSE_df_int),1,nchar(colnames(lfcSE_df_int))-9))

length(which(rownames(LogFC_df_int)!=rownames(BH_df_int)))
length(which(rownames(LogFC_df_int)!=rownames(lfcSE_df_int)))

### Extract the 8 responses to Growth arrest
vec_iron <- c(9:16)
LogFC_df_iron <- LogFC_df[,vec_iron]
BH_df_iron <- BH_df[,vec_iron]
lfcSE_df_iron <- lfcSE_df[,vec_iron]

colnames(LogFC_df_iron) <- paste0(substring(colnames(LogFC_df_iron),1,nchar(colnames(LogFC_df_iron))-24),".LogFC")
# colnames(BH_df_int) <- paste0(colnames(int_BH),".BH")
colnames(lfcSE_df_iron) <- paste0(substring(colnames(lfcSE_df_iron),1,nchar(colnames(lfcSE_df_iron))-9))

length(which(rownames(LogFC_df_iron)!=rownames(BH_df_iron)))
length(which(rownames(LogFC_df_iron)!=rownames(lfcSE_df_iron)))

###########################
#### 3. Check Coherence ###
###########################

threshold=0.05
threshold_int=0.05

name <- c("C1_C2_Glycerol-Dextrose","C5_C6_Glycerol","C7_C8_Glycerol-LCFA","C11_C12_LCFA")

vec <- c(9:12)

Genes_list_enhanced <- list()

for (i in vec) {

  indice <- i-8
  print(paste0("Calculating ",name[indice], " layer"))

  summary=data.frame(
    logFC_with=LogFC_df[,i],
    logFC_without=LogFC_df[,i+4],
    logFC_int=LogFC_df[,i+8],
    BH_with=BH_df[,i],
    BH_without=BH_df[,i+4],
    BH_int=BH_df[,i+8])
  rownames(summary)=rownames(LogFC_df)

  summary$gene <- rownames(summary)

  {
    summary$label <- "Background"
    summary$label[which((summary$BH_with<threshold & summary$BH_without<threshold) &
                          summary$logFC_with*summary$logFC_without>0)]="Coherent response"

    summary$label[which((summary$BH_with<threshold & summary$BH_without<threshold) &
                          summary$logFC_with*summary$logFC_without<0)]="Flipped response"

    summary=summary[order(summary$label),]

    test=cor.test(summary$logFC_with,summary$logFC_without)
    test$estimate
    # 0.8945709
    test$p.value
    # 0

    length(which(summary$label=="Coherent response"))
    # 2348
    length(which(summary$label=="Coherent response"))/(nrow(summary)-length(which(summary$label=="Background")))
    # 0.9779259
    length(which(summary$BH_int<0.05))
    # 1160
    length(which(summary$BH_without<threshold))
    # 2927

    #Intercept again
    length(which(summary$BH_with<threshold & summary$BH_without<threshold))
    # 2401
  }
  #############################################
  ### Absolute effect sizes  ###
  #############################################
  {
    threshold_int=0.05
    threshold=0.01

    {
      summary$label_int="Background"
      summary$label_int[which( summary$logFC_with>0 & summary$logFC_without>0)]="UP"
      summary$label_int[which( summary$logFC_with>0 & summary$logFC_without>0 &
                                 summary$BH_int<threshold_int & summary$logFC_int>0)]="UP_More_without_Fe"
      summary$label_int[which( summary$logFC_with>0 & summary$logFC_without>0 &
                                 summary$BH_int<threshold_int & summary$logFC_int<0)]="UP_More_with_Fe"

      summary$label_int[which( summary$logFC_with<0 & summary$logFC_without<0)]="DOWN"
      summary$label_int[which( summary$logFC_with<0 & summary$logFC_without<0 &
                                 summary$BH_int<threshold_int & summary$logFC_int>0)]="DOWN_More_with_Fe"
      summary$label_int[which( summary$logFC_with<0 & summary$logFC_without<0 &
                                 summary$BH_int<threshold_int & summary$logFC_int<0)]="DOWN_More_without_Fe"

      summary$label_int[which(summary$logFC_with<0 & summary$logFC_without>0 & summary$BH_int<threshold_int)]="FLIPPED_int_down_up_weak"
      summary$label_int[which(summary$logFC_with>0 & summary$logFC_without<0 & summary$BH_int<threshold_int)]="FLIPPED_int_up_down_weak"

      summary$label_int[which(summary$BH_with<threshold & summary$BH_without<threshold & summary$logFC_with<0 & summary$logFC_without>0 & summary$BH_int<threshold_int)]="FLIPPED_int_down_up_strong"
      summary$label_int[which(summary$BH_with<threshold & summary$BH_without<threshold & summary$logFC_with>0 & summary$logFC_without<0 & summary$BH_int<threshold_int)]="FLIPPED_int_up_down_strong"
    }

    summary_df <- data.frame(summary(factor(summary$label_int)))

    up_chunk=summary[which(summary$label_int %in% c("UP_More_with_Fe","UP_More_without_Fe")),]
    up_chunk$int_rebuilt=abs(up_chunk$logFC_without)-abs(up_chunk$logFC_with)
    up_chunk$label_int="Upregulated_genes"

    down_chunk=summary[which(summary$label_int %in% c("DOWN_More_with_Fe","DOWN_More_without_Fe")),]
    down_chunk$int_rebuilt=abs(down_chunk$logFC_without)-abs(down_chunk$logFC_with)
    down_chunk$label_int="Downregulated_genes"

    # tab=rbind(up_chunk,down_chunk)

    ### Save the genes enhanced_upregulated, dampened_upreg, enhanced_downregulated and dampened_downreg
    # enhanced_upregulated <- rownames(summary[which(summary$label_int=="UP_More_without_Fe"),])
    # damped_upreg <- rownames(summary[which(summary$label_int=="UP_More_with_Fe"),])
    #
    # enhanced_downregulated <- rownames(summary[which(summary$label_int=="DOWN_More_without_Fe"),])
    # damped_downreg <- rownames(summary[which(summary$label_int=="DOWN_More_with_Fe"),])

    enhanced_upregulated <- summary[
      summary$label_int == "UP_More_without_Fe",
      c("gene","logFC_int","BH_int" )]

    damped_upreg <- summary[
      summary$label_int == "UP_More_with_Fe",
      c("gene", "logFC_int","BH_int" )]

    enhanced_downregulated <- summary[
      summary$label_int == "DOWN_More_without_Fe",
      c("gene", "logFC_int","BH_int" )]

    damped_downreg <- summary[
      summary$label_int == "DOWN_More_with_Fe",
      c("gene","logFC_int","BH_int" )]

    dir.create(file.path(output_dir,"Enhanced_and_damped_genes_regulon",name[indice]),recursive = T,showWarnings = F)
    file_name <- c("enhanced_upregulated.txt","damped_upreg.txt","enhanced_downregulated.txt","damped_downreg.txt")
    dir <- file.path(output_dir,"Enhanced_and_damped_genes_regulon",name[indice])

    Genes_list_enhanced[[name[indice]]] <- list(
      EU = enhanced_upregulated,      # Enhanced Up
      ED = enhanced_downregulated,    # Enhanced Down
      DU = damped_upreg,              # Damped Up
      DD = damped_downreg             # Damped Down
    )

    write.table(enhanced_upregulated,file.path(dir, file_name[1]),sep = "\t",row.names = FALSE,quote = FALSE)

    write.table(damped_upreg,file.path(dir, file_name[2]),sep = "\t",row.names = FALSE,quote = FALSE)

    write.table(enhanced_downregulated,file.path(dir, file_name[3]),sep = "\t",row.names = FALSE,quote = FALSE)

    write.table(damped_downreg,file.path(dir, file_name[4]),sep = "\t",row.names = FALSE,quote = FALSE)
  }
}

###########################
##  Enrichment Analysis  ##
###########################

# Load regulons
regulons    <- read.table(file.path(main_wd,"Inputs","001_Raw_data","Minch_regulons.txt"),header = TRUE,sep = "\t",stringsAsFactors = FALSE)
regulons_JS <- openxlsx::read.xlsx(file.path(main_wd,"Inputs","001_Raw_data","Signed_TRN_MTB.xlsx"),sheet = "Network")
regulons_JS <- regulons_JS[,c("Regulator_Rv_index","Target_Rv_index")]

############################
### Merge regulon genes ####
############################

# 1. Primer regulon (vector)
regulon_genes <- unique(regulons[[1]])

# 2. Genes de regulons_JS (ajusta nombre de columna si es necesario)
regulon_genes_js <- unique(regulons_JS$Regulator_Rv_index)
regulon_genes_js <- regulon_genes_js[!regulon_genes_js %in% regulon_genes]

regulons_JS_filtered <- regulons_JS[regulons_JS$Regulator_Rv_index %in% regulon_genes_js, ]
length(unique(regulons_JS_filtered$Regulator_Rv_index))

# # 4. Unión: regulon + genes de JS que no estén en el regulon
# final_genes <- c(regulon_genes,regulon_genes_js)
#
# # 5. Eliminar posibles duplicados (por seguridad)
# final_genes <- unique(final_genes)

regulons_JS_filtered <- regulons_JS_filtered %>%
  rename(Regulator = Regulator_Rv_index, gene = Target_Rv_index)

regulons_final <- rbind(
  regulons[, c("Regulator", "gene")],
  regulons_JS_filtered)

# Map locus_tags to gene symbols
gene_map        <- feature_data$symbol
names(gene_map) <- feature_data$locus_tag

# fullfill missing mappings with locus_tags (if symbol is NA or empty, keep locus_tag)
gene_map[is.na(gene_map) | gene_map == ""] <- names(gene_map)[is.na(gene_map) | gene_map == ""]

regulons_final$gene_name <- gene_map[(regulons_final$Regulator)]
regulons_final$gene_name[is.na(regulons_final$gene_name)] <- regulons_final$Regulator[is.na(regulons_final$gene_name)]

# Create a list of regulons
regulon_list <- split(regulons_final$gene, regulons_final$Regulator)

# regulons$gene_name <- gene_map[regulons$Regulator]
# regulons$gene_name[is.na(regulons$gene_name)] <- regulons$Regulator[is.na(regulons$gene_name)]

############################################################
##  Performe an ORA for each condition and each gene set  ##
##  (enhanced up, enhanced down)                          ##
############################################################
set.seed(12345)

REGULON_ORA <- function(query,
                        TERM2GENE,
                        universe,
                        test_mode = c("enrichment","depletion","two.sided"),
                        pvalue_correction_method = "BH",
                        min_term_size = 1,
                        max_term_size = Inf,
                        merge_duplicates = TRUE,
                        add_haldane = TRUE) {

  test_mode <- match.arg(test_mode)

  # -------------------------------
  # QC + estandarizar IDs
  # -------------------------------
  stopifnot(all(c("gene_name","gene") %in% colnames(TERM2GENE)))

  # query    <- toupper(unique(query))
  # universe <- toupper(unique(universe))
  #
  query    <- unique(query)
  universe <- unique(universe)

  TERM2GENE <- TERM2GENE[, c("gene_name","gene")]
  TERM2GENE$gene_name <- (TERM2GENE$gene_name)
  TERM2GENE$gene      <- (TERM2GENE$gene)
  TERM2GENE <- TERM2GENE[!duplicated(TERM2GENE), ]

  # Mantener solo genes del universo
  TERM2GENE <- TERM2GENE[TERM2GENE$gene %in% universe, ]

  # Query dentro del universo
  query <- query[query %in% universe]
  if (length(query) == 0) stop("After intersecting with universe, query is empty.")

  # -------------------------------
  # regulon_list + tamaños
  # -------------------------------
  regulon_list <- split(TERM2GENE$gene, TERM2GENE$gene_name)
  regulon_list <- lapply(regulon_list, unique)

  term_sizes <- vapply(regulon_list, length, integer(1))
  keep_terms <- names(term_sizes)[term_sizes >= min_term_size & term_sizes <= max_term_size]

  regulon_list <- regulon_list[keep_terms]
  term_sizes   <- term_sizes[keep_terms]

  if (length(regulon_list) == 0) stop("No regulons left after size filtering.")

  # -------------------------------
  # Fisher por regulón
  # -------------------------------
  mode <- switch(test_mode,
                 "enrichment" = "greater",
                 "depletion"  = "less",
                 "two.sided"  = "two.sided")

  N <- length(universe)
  n <- length(query)

  regs <- names(regulon_list)

  pvals   <- numeric(length(regs))
  ORs     <- numeric(length(regs))
  ci_low  <- numeric(length(regs))
  ci_high <- numeric(length(regs))
  k_vec   <- integer(length(regs))   # overlap
  K_vec   <- integer(length(regs))   # size in universe
  members <- character(length(regs))

  for (i in seq_along(regs)) {
    reg <- regs[i]
    reg_genes <- regulon_list[[reg]]

    k <- sum(query %in% reg_genes)         # in both
    K <- length(reg_genes)                 # in regulon within universe
    k_vec[i] <- k
    K_vec[i] <- K

    # 2x2:
    #            In_reg  Not_in_reg
    # In_query      k      n-k
    # Not_in_q    K-k   N-K-n+k
    a <- k
    b <- n - k
    c <- K - k
    d <- N - K - n + k

    tab <- matrix(c(a, c, b, d), nrow = 2, byrow = TRUE)

    ft <- fisher.test(tab, alternative = mode)

    pvals[i]   <- ft$p.value
    ORs[i]     <- unname(ft$estimate)
    ci_low[i]  <- ft$conf.int[1]
    ci_high[i] <- ft$conf.int[2]

    members[i] <- paste(intersect(query, reg_genes), collapse = ", ")
  }

  # Haldane/Anscombe correction SOLO para log2OR si quieres evitar 0/Inf raros
  if (add_haldane) {
    # OR corregido con +0.5 en cada celda
    ORs_h <- numeric(length(regs))
    for (i in seq_along(regs)) {
      reg_genes <- regulon_list[[regs[i]]]
      k <- sum(query %in% reg_genes)
      K <- length(reg_genes)
      a <- k; b <- n - k; c <- K - k; d <- N - K - n + k
      ORs_h[i] <- ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
    }
    log2OR <- log2(ORs_h)
  } else {
    log2OR <- log2(ORs)
  }

  fdr <- p.adjust(pvals, method = pvalue_correction_method)

  results <- data.frame(
    Regulator = regs,
    Size_in_universe = K_vec,
    Number_in_query  = k_vec,
    Percentage       = (k_vec / K_vec) * 100,
    p   = pvals,
    fdr = fdr,
    OR  = ORs,
    Min_value_CI = ci_low,
    Max_value_CI = ci_high,
    Members = members,
    stringsAsFactors = FALSE
  )

  # Duplicados = regulones con el mismo set de genes
  if (merge_duplicates) {
    signature <- vapply(regulon_list[results$Regulator],
                        function(x) paste(sort(x), collapse = "|"),
                        character(1))
    results$uniqueness_index <- as.integer(factor(signature))

    # quedarte con el mejor por firma (menor fdr)
    results <- results[order(results$fdr), ]
    results <- results[!duplicated(results$uniqueness_index), ]
  } else {
    results$uniqueness_index <- seq_len(nrow(results))
  }

  results$log2OR <- log2OR[match(results$Regulator, regs)]
  results <- results[order(results$fdr, -abs(results$log2OR)), ]

  list(
    universe = universe,
    query = query,
    TERM2GENE = TERM2GENE,
    regulon_list = regulon_list,
    result = results
  )
}
# Define universe as all genes tested in the differential expression analysis
universe <- rownames(LogFC_df)

TERM2GENE <- regulons_final[, c("gene_name","gene")] ### these have to be the names of the columns in order to use the ORA  function

ora_results <- list()

for (condition in names(Genes_list_enhanced)) {

  ora_results[[condition]] <- list()

  for (gene_set in names(Genes_list_enhanced[[condition]])) {

    genes_of_interest <- Genes_list_enhanced[[condition]][[gene_set]]$gene

    out <- REGULON_ORA(
      query      = genes_of_interest,
      TERM2GENE  = TERM2GENE,
      universe   = universe,
      test_mode  = "enrichment",
      pvalue_correction_method = "BH",
      min_term_size = 1,
      max_term_size = Inf,
      merge_duplicates = TRUE,
      add_haldane = TRUE
    )

    res <- out$result
    res_sig <- subset(res, fdr < 0.1)

    ora_results[[condition]][[gene_set]] <- list(
      all = res,
      sig = res_sig
    )

    dir.create(file.path(output_dir, "ORA_regulons", condition),
               recursive = TRUE, showWarnings = FALSE)

    write.table(
      res,
      file.path(output_dir, "ORA_regulons", condition, paste0("ORA_", gene_set, ".txt")),
      sep = "\t", row.names = FALSE, quote = FALSE
    )
  }
}

saveRDS(ora_results, file.path(output_dir, "ORA_regulons", "ORA_object.RDS"))

###########################################
##  Create the heatmaps for ORA results  ##
###########################################

mat_EU <- build_top5_df(ora_results, gene_set = "EU")
mat_ED <- build_top5_df(ora_results, gene_set = "ED")

cond_rename <- c(
  "C1_C2_Glycerol-Dextrose" = "G-D",
  "C5_C6_Glycerol"          = "G",
  "C7_C8_Glycerol-LCFA"     = "G-LCFA",
  "C11_C12_LCFA"            = "LCFA"
)

mat_EU_signed <- as.matrix(mat_EU[, -1])
rownames(mat_EU_signed) <- mat_EU$Regulon

mat_ED_signed <- as.matrix(mat_ED[, -1])
rownames(mat_ED_signed) <- mat_ED$Regulon

colnames(mat_EU_signed) <- cond_rename[colnames(mat_EU_signed)]
colnames(mat_ED_signed) <- cond_rename[colnames(mat_ED_signed)]

plot_ora_heatmap_gg <- function(ora_results,
                                gene_set = "EU",
                                top_n = 5,
                                sort_by = c("fdr","abs_log2OR"),
                                mx = 6,                 # máximo fijo para escala de color
                                fdr_cap = 1e-6,         # para que -log10 no explote
                                na_fdr = 1,             # NA fdr -> 1 (punto pequeño)
                                na_for_dist = 0,        # NA SOLO para clustering
                                high_col = "red3",      # EU: rojo; ED: azul
                                mid_col  = "grey95",
                                low_col  = "white",
                                base_size = 9) {

  sort_by <- match.arg(sort_by)

  # 1) Define regulones a mostrar (top_n por condición, unidos)
  regs_keep <- get_top_regulons_union(ora_results, gene_set, top_n, sort_by)
  if (length(regs_keep) == 0) stop("No regulons found for this gene_set/top_n.")

  # 2) Matrices: log2OR + fdr
  mat_or  <- extract_metric_matrix(ora_results, gene_set, value_col = "log2OR", key_col = "Regulator")
  mat_fdr <- extract_metric_matrix(ora_results, gene_set, value_col = "fdr",    key_col = "Regulator")

  # Alinear y subset a regulones de interés
  mat_or  <- mat_or[intersect(regs_keep, rownames(mat_or)), , drop = FALSE]
  mat_fdr <- mat_fdr[rownames(mat_or), colnames(mat_or), drop = FALSE]  # mismo orden

  # 3) Clustering filas (sin dendrograma, solo orden)
  mat_dist <- mat_or
  mat_dist[is.na(mat_dist) | !is.finite(mat_dist)] <- na_for_dist
  hc <- hclust(dist(mat_dist), method = "ward.D2")
  row_order <- rownames(mat_or)[hc$order]

  # 4) Long format para ggplot
  df_or <- as.data.frame(mat_or) %>%
    mutate(Regulon = rownames(.)) %>%
    pivot_longer(-Regulon, names_to = "Condition", values_to = "log2OR")

  df_fdr <- as.data.frame(mat_fdr) %>%
    mutate(Regulon = rownames(.)) %>%
    pivot_longer(-Regulon, names_to = "Condition", values_to = "fdr")

  df <- df_or %>%
    left_join(df_fdr, by = c("Regulon","Condition"))

  # Factor order = clustering
  df$Regulon   <- factor(df$Regulon, levels = row_order, ordered = TRUE)
  df$Condition <- factor(df$Condition, levels = colnames(mat_or))

  # fdr: NA -> 1 y cap mínimo
  df$fdr_plot <- df$fdr
  df$fdr_plot[is.na(df$fdr_plot) | !is.finite(df$fdr_plot)] <- na_fdr
  df$fdr_plot <- pmax(df$fdr_plot, fdr_cap)

  # 5) Plot (tile + punto como en tu ejemplo)
  p <- ggplot(df, aes(x = Condition, y = Regulon)) +
    geom_tile(color = "black", fill = "grey99", linewidth = 0.6) +
    geom_point(aes(size = -log10(fdr_plot), color = log2OR), shape = 15) +
    coord_equal() +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    # scale_size(range = c(1, 5), name = expression(-log[10](FDR))) +
    scale_size(
      range = c(1, 5),
      name  = "FDR",
      breaks = c(1, 2, 3, 4, 5),        # en escala -log10
      labels = function(x) {
        formatC(10^(-x), format = "E", digits = 1)
      }
    ) +
    scale_color_gradient2(
      low = low_col, mid = mid_col, high = high_col,
      midpoint = 0,
      limits = c(-mx, mx),
      oob = scales::squish,
      name = "log2(OR)"
    ) +
    theme_tufte(base_size = base_size, base_family = "Helvetica") +
    theme(
      axis.ticks = element_blank(),
      axis.text.x = element_text(size = base_size * 0.85, angle = 270, hjust = 0, colour = "grey40"),
      axis.text.y = element_text(size = base_size * 0.8),
      legend.position = "left",
      plot.margin = unit(c(t = 0.2, r = 0.2, b = 0.2, l = 0.2), "cm")
    ) +
    labs(x = NULL, y = NULL)

  list(plot = p, mat_or = mat_or, mat_fdr = mat_fdr, hc = hc)
}

save_ora_heatmap_pdf <- function(p, file, width = 8, height = 10) {
  pdf(file, width = width, height = height)
  print(p)
  dev.off()
}

ora_dir <- file.path(output_dir, "ORA_regulons")

mx_fixed <- 5  # pon un valor fijo para EU/ED (o calcula uno global y luego lo fijas)

# EU (rojo)
eu <- plot_ora_heatmap_gg(
  ora_results,
  gene_set = "EU",
  top_n = 5,
  sort_by = "fdr",
  mx = mx_fixed,
  high_col = "red3",
  low_col = "#0C6291"
)

save_ora_heatmap_pdf(
  eu$plot,
  file.path(ora_dir, "Heatmap_EU_gg_no_dend.pdf"),
  width = 8, height = 10
)

# ED (azul)
ed <- plot_ora_heatmap_gg(
  ora_results,
  gene_set = "ED",
  top_n = 5,
  sort_by = "fdr",
  mx = mx_fixed,
  high_col = "red3",
  low_col = "#0C6291"
)

save_ora_heatmap_pdf(
  ed$plot,
  file.path(ora_dir, "Heatmap_ED_gg_no_dend.pdf"),
  width = 8, height = 10
)

#########################
##### Other options:##### this is the one i used for the final version of the paper, with clustering and dendrograms.
#########################
plot_ora_heatmap_ch <- function(ora_results,
                                gene_set = "EU",
                                top_n = 5,
                                sort_by = c("fdr", "abs_log2OR"),
                                high_col = "red3",     # EU: rojo; ED: azul
                                mid_col  = "#fcae91",  # color intermedio
                                low_col  = "white",
                                k = 6,
                                na_shift = 0.5,        # cuánto bajar los NA en clustering
                                cell_fontsize = 3,
                                row_fontsize = 6,
                                col_fontsize = 6,
                                legend_fontsize = 6,
                                legend_title_fontsize = 8,
                                heatmap_width = unit(2, "cm"),
                                column_title = "Conditions",
                                row_title = "Regulons") {

  sort_by <- match.arg(sort_by)

  # 1) Regulones a mostrar
  regs_keep <- get_top_regulons_union(ora_results, gene_set, top_n, sort_by)
  if (length(regs_keep) == 0) stop("No regulons found for this gene_set/top_n.")

  # 2) Matrices: log2OR + fdr
  mat_or  <- extract_metric_matrix(ora_results, gene_set, value_col = "log2OR", key_col = "Regulator")
  mat_fdr <- extract_metric_matrix(ora_results, gene_set, value_col = "fdr",    key_col = "Regulator")

  # Alinear y subset
  regs_keep <- intersect(regs_keep, rownames(mat_or))
  if (length(regs_keep) == 0) stop("None of the selected regulons are present in mat_or.")

  mat_or  <- mat_or[regs_keep, , drop = FALSE]
  mat_fdr <- mat_fdr[rownames(mat_or), colnames(mat_or), drop = FALSE]

  # 3) Preparar matriz de plotting con el mismo tratamiento que en los heatmaps anteriores
  plot_matrix <- mat_or

  pos_inf <- is.infinite(plot_matrix) & plot_matrix > 0
  neg_inf <- is.infinite(plot_matrix) & plot_matrix < 0
  nan_val <- is.nan(plot_matrix)

  max_cap <- ceiling(max(plot_matrix[is.finite(plot_matrix)], na.rm = TRUE))

  plot_matrix[pos_inf] <- max_cap
  plot_matrix[neg_inf] <- NA
  plot_matrix[nan_val] <- NA

  # 4) Matriz para clustering
  cluster_matrix <- plot_matrix
  min_finite <- min(cluster_matrix[is.finite(cluster_matrix)], na.rm = TRUE)
  cluster_matrix[is.na(cluster_matrix)] <- min_finite - na_shift

  dist_matrix <- dist(cluster_matrix)
  hc <- hclust(dist_matrix, method = "ward.D2")

  # 5) Escala de color como en los heatmaps anteriores
  break_vec <- c(0, (max_cap / 2) - 2, max_cap - 2)

  # Seguridad por si max_cap es pequeño
  break_vec <- sort(unique(pmax(break_vec, 0)))
  if (length(break_vec) == 1) break_vec <- c(0, max_cap / 2, max_cap)
  if (length(break_vec) == 2) break_vec <- c(break_vec[1], mean(break_vec), break_vec[2])

  col_fun <- circlize::colorRamp2(
    break_vec,
    c(low_col, mid_col, high_col)
  )

  print(paste0(" the names of the columns in plot_matrix before change are: ", paste(colnames(plot_matrix), collapse = ", ")))
  colnames(plot_matrix) <- c("G-D", "G", "G-LCFA", "LCFA")
  print(paste0(" the names of the columns in plot_matrix after change are: ", paste(colnames(plot_matrix), collapse = ", ")))

  # 6) Heatmap
  ht <- ComplexHeatmap::Heatmap(
    plot_matrix,
    na_col = "white",
    col = col_fun,
    split = k,
    name = "Log2(OR)",
    column_order = colnames(plot_matrix),
    show_column_names = TRUE,
    column_names_gp = grid::gpar(fontsize = col_fontsize),
    row_names_gp = grid::gpar(fontsize = row_fontsize),
    column_title = column_title,
    column_title_side = "bottom",
    border_gp = grid::gpar(col = "black", lty = 2),
    width = heatmap_width,
    cluster_rows = dendextend::color_branches(hc, k = k),
    row_title = row_title,
    row_title_gp = grid::gpar(fontsize = row_fontsize),
    row_dend_side = "right",
    row_names_side = "left",
    show_row_names = TRUE,
    show_row_dend = TRUE,
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (!is.na(plot_matrix[i, j])) {
        grid::grid.text(
          round(plot_matrix[i, j], 2),
          x, y,
          gp = grid::gpar(fontsize = cell_fontsize)
        )
      }
    },
    heatmap_legend_param = list(
      title = "Log2(OR)",
      at = break_vec,
      labels = round(break_vec, 2),
      title_position = "leftcenter-rot",
      labels_gp = grid::gpar(fontsize = legend_fontsize),
      title_gp = grid::gpar(fontsize = legend_title_fontsize)
    )
  )

  # 7) Salida
  list(
    ht = ht,
    plot_matrix = plot_matrix,
    mat_fdr = mat_fdr,
    cluster_matrix = cluster_matrix,
    hc = hc,
    break_vec = break_vec,
    max_cap = max_cap
  )
}

res_eu <- plot_ora_heatmap_ch(
  ora_results = ora_results,
  gene_set = "EU",
  top_n = 5,
  sort_by = "fdr",
  high_col = "#cb181d",
  mid_col = "#fcae91",
  low_col = "white",
  k = 6,
  row_title = "Regulons"
)

ComplexHeatmap::draw(res_eu$ht)

res_ed <- plot_ora_heatmap_ch(
  ora_results = ora_results,
  gene_set = "ED",
  top_n = 5,
  sort_by = "fdr",
  high_col = "#08519c",
  mid_col = "#9ecae1",
  low_col = "white",
  k = 4,
  row_title = "Regulons"
)

ComplexHeatmap::draw(res_ed$ht)

size <- calc_ht_size(res_eu$ht)

### Save as svg for better quality in the paper



pdf(file.path(ora_dir, "ORA_heatmap_EU.pdf"), width = size[1] + 0.5, height = size[2] + 1)
ComplexHeatmap::draw(res_eu$ht)
dev.off()

# Install once if needed
# install.packages("svglite")

svg_file <- file.path(ora_dir, "ORA_heatmap_EU.svg")

svglite::svglite(
  file = svg_file,
  width = size[1] + 0.5,
  height = size[2] + 1
)

ComplexHeatmap::draw(res_eu$ht)

dev.off()

size <- calc_ht_size(res_ed$ht)
pdf(file.path(ora_dir, "ORA_heatmap_ED.pdf"), width = size[1] + 0.5, height = size[2] + 1)
ComplexHeatmap::draw(res_ed$ht)
dev.off()

svg_file_ed <- file.path(ora_dir, "ORA_heatmap_ED.svg")

svglite::svglite(
  file = svg_file_ed,
  width = size[1] + 0.5,
  height = size[2] + 1
)

ComplexHeatmap::draw(res_ed$ht)

dev.off()
