#############################
### 0. Load dependencies. ###
#############################

library(tidyverse)
library(ggplot2)
library(dplyr)
library(patchwork)

############################
### 1. Set directories #####
############################

main_wd <- getwd()
setwd(main_wd)

input_dir  <- "Inputs/002_Processed_data"
output_dir <- "Outputs"

iron_response_to_stat <- read.table(
  file.path(input_dir, "txt/iron_response_to_stat.txt"),
  check.names = FALSE
)

th <- 0.05

########################################
### 2. Define short contrast names #####
########################################

short_names <- c(
  "G-D", "G", "G-LCFA", "LCFA"
)

names(short_names) <- colnames(iron_response_to_stat)[1:4]

# For the "without Fe" block, reuse the same short labels
short_names_without <- c(
  "G-D", "G", "G-LCFA", "LCFA"
)

names(short_names_without) <- colnames(iron_response_to_stat)[5:8]

################################
### 3. Pair combinations #######
################################

combination_names_with    <- combn(colnames(iron_response_to_stat)[1:4], m = 2)
combination_names_without <- combn(colnames(iron_response_to_stat)[5:8], m = 2)

########################################
### 4. Helper plotting functions #######
########################################

make_corr_plot <- function(df_filtered, x_var, y_var, x_lab, y_lab, point_fill = "grey80") {

  test <- cor.test(df_filtered[[x_var]], df_filtered[[y_var]], method = "pearson")

  xy_range <- range(c(df_filtered[[x_var]], df_filtered[[y_var]]), na.rm = TRUE)
  lim_max <- max(abs(xy_range))
  axis_limits <- c(-lim_max, lim_max)

  p <- ggplot(df_filtered, aes(x = .data[[x_var]], y = .data[[y_var]])) +
    geom_point(
      shape = 21,
      size = 1.6,
      stroke = 0.25,
      colour = "black",
      fill = point_fill,
      alpha = 0.75
    ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.4, colour = "grey50") +
    geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = FALSE,
      linewidth = 0.7,
      colour = "black"
    ) +
    geom_vline(xintercept = 0, linewidth = 0.35, colour = "grey40") +
    geom_hline(yintercept = 0, linewidth = 0.35, colour = "grey40") +
    coord_cartesian(xlim = axis_limits, ylim = axis_limits) +
    labs(
      title = paste0(x_lab, " vs ", y_lab),
      subtitle = paste0("R = ", round(test$estimate, 2), "   |   n = ", nrow(df_filtered)),
      x = x_lab,
      y = y_lab
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
      plot.subtitle = element_text(size = 9, hjust = 0.5),
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 8),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7),
      axis.line = element_blank()
    )

  list(plot = p, cor = test$estimate, pval = test$p.value)
}
strip_axes <- function(p, remove_x = FALSE, remove_y = FALSE) {

  if (remove_x) {
    p <- p +
      labs(x = NULL) +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
  }

  if (remove_y) {
    p <- p +
      labs(y = NULL) +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      )
  }

  p
}

#########################################
### 5. Build one plot collection ########
#########################################

build_plot_set <- function(combination_names, short_map, point_fill) {

  plot_list <- list()
  cor_matrix <- matrix(NA, nrow = ncol(combination_names), ncol = 1)
  pvalues_matrix <- matrix(NA, nrow = ncol(combination_names), ncol = 1)
  names_cor_matrix <- character(ncol(combination_names))

  for (col in seq_len(ncol(combination_names))) {

    var1 <- combination_names[1, col]
    var2 <- combination_names[2, col]

    df <- data.frame(
      LogFC1 = iron_response_to_stat[, var1],
      LogFC2 = iron_response_to_stat[, var2],
      BH1    = iron_response_to_stat[, paste0(var1, ".BH")],
      BH2    = iron_response_to_stat[, paste0(var2, ".BH")],
      row.names = rownames(iron_response_to_stat)
    )

    df_filtered <- df[df$BH1 < th & df$BH2 < th, , drop = FALSE]

    short1 <- short_map[[var1]]
    short2 <- short_map[[var2]]

    # Save filtered data
    df_to_corr_dir <- paste0(short1, "_vs_", short2, "_Fig_S1B.txt")
    dir.create(file.path(input_dir, "txt/Correlogram_data_Fig_S1B"), showWarnings = FALSE, recursive = TRUE)
    write.table(
      df_filtered,
      file.path(input_dir, "txt/Correlogram_data_Fig_S1B", df_to_corr_dir)
    )

    # Correlation
    test <- cor.test(df_filtered$LogFC1, df_filtered$LogFC2, method = "pearson")
    names_cor_matrix[col] <- paste0(short1, "_vs_", short2)
    cor_matrix[col, 1] <- test$estimate
    pvalues_matrix[col, 1] <- -log10(test$p.value)

    # Plot
    plot_obj <- make_corr_plot(
      df_filtered = df_filtered,
      x_var = "LogFC1",
      y_var = "LogFC2",
      x_lab = short1,
      y_lab = short2,
      point_fill = point_fill
    )

    p <- plot_obj$plot
    plot_list[[paste0(short1, "_vs_", short2)]] <- p

    fig_repo <- "001_Figures_paper"
    scatter_dir <- "4_Figure_S1B_Scatter_plot_correlogram"
    dir.create(file.path(output_dir, fig_repo, scatter_dir), recursive = TRUE, showWarnings = FALSE)

    pdf(file.path(output_dir, fig_repo, scatter_dir, paste0(short1, "_vs_", short2, ".pdf")),
        width = 4.2, height = 4.0)
    print(p)
    dev.off()
  }

  rownames(cor_matrix) <- names_cor_matrix
  rownames(pvalues_matrix) <- names_cor_matrix

  list(
    plots = plot_list,
    cor_matrix = cor_matrix,
    pvalues_matrix = pvalues_matrix
  )
}

res_with <- build_plot_set(
  combination_names = combination_names_with,
  short_map = short_names,
  point_fill = "grey70"
)

res_without <- build_plot_set(
  combination_names = combination_names_without,
  short_map = short_names_without,
  point_fill = "white"
)

make_staircase_panel <- function(plot_list) {

  p1 <- strip_axes(plot_list[["G-D_vs_G"]],         remove_x = TRUE,  remove_y = FALSE)
  p2 <- strip_axes(plot_list[["G-D_vs_G-LCFA"]],    remove_x = TRUE,  remove_y = TRUE)
  p3 <- strip_axes(plot_list[["G-D_vs_LCFA"]],      remove_x = TRUE,  remove_y = TRUE)
  p4 <- strip_axes(plot_list[["G_vs_G-LCFA"]],      remove_x = TRUE,  remove_y = FALSE)
  p5 <- strip_axes(plot_list[["G_vs_LCFA"]],        remove_x = TRUE,  remove_y = TRUE)
  p6 <- strip_axes(plot_list[["G-LCFA_vs_LCFA"]],   remove_x = FALSE, remove_y = FALSE)

  layout <- "
  ABC
  DEF"
  # layout <- "
  # ABC
  # DE#
  # F##
  # "
  wrap_plots(
    A = p1,
    B = p2,
    C = p3,
    D = p4,
    E = p5,
    F = p6,
    design = layout
  )
}
panel_with <- make_staircase_panel(res_with$plots)
panel_without <- make_staircase_panel(res_without$plots)

fig_repo <- "001_Figures_paper"
scatter_dir <- "4_Figure_S1B_Scatter_plot_correlogram"

pdf(file.path(output_dir, fig_repo, scatter_dir, "Correlogram_with_Fe_staircase.pdf"),
    width = 11, height = 8.5)
print(panel_with)
dev.off()

pdf(file.path(output_dir, fig_repo, scatter_dir, "Correlogram_without_Fe_staircase.pdf"),
    width = 11, height = 8.5)
print(panel_without)
dev.off()

final_panel <- panel_with / panel_without

pdf(file.path(output_dir, fig_repo, scatter_dir, "Correlogram_all_staircase.pdf"),
    width = 12, height = 16)
print(final_panel)
dev.off()

### Save in svg format for paper
svglite::svglite(file.path(output_dir, fig_repo, scatter_dir, "Correlogram_all_staircase.svg"),
    width = 12, height = 16)
