######################
### 0.Dependencies ###
######################
{
  library(tidyverse)
  library(shadowtext)
  library(cowplot)
}

###########################
### 1.Declare Functions ###
###########################

dcols=function(x){data.frame(colnames(x))}

#################################
### 2. Set working directory  ###
#################################

main_wd <- getwd()
setwd(main_wd)
input_dir <- "Inputs/002_Processed_data"
output_dir <- "Outputs"

#####################
### 3. Load data  ###
#####################

### feature data filtered
feat_path = file.path(input_dir,"txt/feature_data_filtered.txt")
feature_data = read.table(feat_path)

# DE genes per condition
DE_genes_per_cond <- read.table(file = file.path(input_dir,"txt/DE_genes_per_cond_th_0.05.txt"))
# ncRNA reads data(to obtain all the gene names)
# reads_nc          <- read.table(file.path(input_dir,"txt/reads_nc.txt"),check.names = F)
# Iron effects data frame
iron_effects      <- read.table(file.path(input_dir,"txt/Iron_effects.txt"))

# LogFC, BH, lfcSE data frames
LogFC_df <- read.table(file.path(input_dir,"txt/LogFC_0.05.txt"))
BH_df    <- read.table(file.path(input_dir,"txt/BH_0.05.txt"))
lfcSE_df <- read.table(file.path(input_dir,"txt/lfcSE_0.05.txt"))

# First: summarize biotypes:
length(which(feature_data$Type=="CDS"))
#3908
length(which(feature_data$Type=="Stable_RNA"))
#104

iron_vec <- c(1:8)
growth_vec <- c(9:16)

# # Subset the DE genes per condition data frame for growth conditions
# df <- DE_genes_per_cond[growth_vec,]
# df$Conditions <- rownames(df)
# df$Conditions <- factor(df$Conditions,levels = df$Conditions)

# Subset data for growth conditions
LogFC_df_growth <- LogFC_df[,growth_vec]
BH_df_growth    <- BH_df[,growth_vec]
lfcSE_df_growth <- lfcSE_df[,growth_vec]

############################################

### Change the columns names
colnames(LogFC_df_growth) <- paste0(substring(colnames(LogFC_df_growth),1,nchar(colnames(LogFC_df_growth))-24),".LogFC")
colnames(lfcSE_df_growth) <- paste0(substring(colnames(lfcSE_df_growth),1,nchar(colnames(lfcSE_df_growth))-9),".lfcSE")

# check that the rownames are the same
length(which(rownames(LogFC_df_growth)!=rownames(BH_df_growth)))
length(which(rownames(LogFC_df_growth)!=rownames(lfcSE_df_growth)))

# Combine LogFC and BH into a single data frame
growth_effects <- data.frame(LogFC_df_growth,BH_df_growth)
write.table(growth_effects,file.path(input_dir,"txt/Growth_effects.txt"),sep="\t",quote=F)

{
  iron_effects$Pattern_G_D_exp    = "Background"
  iron_effects$Pattern_G_exp      = "Background"
  iron_effects$Pattern_G_L_exp    = "Background"
  iron_effects$Pattern_L_exp      = "Background"
  iron_effects$Pattern_G_D_stat   = "Background"
  iron_effects$Pattern_G_stat     = "Background"
  iron_effects$Pattern_G_L_stat   = "Background"
  iron_effects$Pattern_L_stat     = "Background"
  threshold <- 0.05
}

{ growth_effects$Pattern_G_D_with_Fe  = "Background"
  growth_effects$Pattern_G_with_Fe    = "Background"
  growth_effects$Pattern_G_L_with_Fe  = "Background"
  growth_effects$Pattern_L_with_Fe    = "Background"
  growth_effects$Pattern_G_D_no_Fe    = "Background"
  growth_effects$Pattern_G_no_Fe      = "Background"
  growth_effects$Pattern_G_L_no_Fe    = "Background"
  growth_effects$Pattern_L_no_Fe      = "Background"
}

### Separate by direction (UP or DOWN)
col_logFC <- c(1:8)

for (i in col_logFC){

  #Upregulated per condition
  iron_effects[which((iron_effects[(i+8)] < threshold) &
                       (iron_effects[i] > 0)),(i+16)] <- "UP"
  #Downregulated per condition
  iron_effects[which(iron_effects[(i+8)] < threshold &
                       iron_effects[i] < 0),(i+16)] <- "DOWN"

  #Upregulated per condition in growth
  growth_effects[which((growth_effects[(i+8)] < threshold) &
                       (growth_effects[i] > 0)),(i+16)] <- "UP"
  #Downregulated per condition in growth
  growth_effects[which(growth_effects[(i+8)] < threshold &
                       growth_effects[i] < 0),(i+16)] <- "DOWN"

}

# Create df of coding genes and non-coding genes separately for iron effects and growth effects

length(which(feature_data$locus_tag != rownames(iron_effects)))
#0
length(which(feature_data$locus_tag != rownames(growth_effects)))
#0

iron_effects$Type   <- feature_data$Type
growth_effects$Type <- feature_data$Type

iron_effects_non_coding <- iron_effects[which(iron_effects$Type == "Stable_RNA"),]
iron_effects_coding     <- iron_effects[which(iron_effects$Type == "CDS"),]

growth_effects_non_coding <- growth_effects[which(growth_effects$Type == "Stable_RNA"),]
growth_effects_coding     <- growth_effects[which(growth_effects$Type == "CDS"),]


# Save data frames
DE_tables_dir <- file.path(input_dir,"001_DE_tables")
dir.create(DE_tables_dir,showWarnings = F)

write.table(iron_effects_non_coding,
            file.path(DE_tables_dir,"Iron_effects_non_coding.txt"),
            sep="\t",quote=F)
write.table(iron_effects_coding,
            file.path(DE_tables_dir,"Iron_effects_coding.txt"),
            sep="\t",quote=F)
write.table(growth_effects_non_coding,
            file.path(DE_tables_dir,"Growth_effects_non_coding.txt"),
            sep="\t",quote=F)
write.table(growth_effects_coding,
            file.path(DE_tables_dir,"Growth_effects_coding.txt"),
            sep="\t",quote=F)

#################
##  Bar plots  ##
#################

# Create a list of data frames for further for looping if needed
effects_list <- list(iron_effects_non_coding,
                     iron_effects_coding,
                     growth_effects_non_coding,
                     growth_effects_coding,
                     iron_effects,
                     growth_effects)

# Clean names of the list
names(effects_list) <- c(
  "Iron (non-coding)",
  "Iron (coding)",
  "Growth (non-coding)",
  "Growth (coding)",
  "Iron deprivation effects",
  "Growth arrest effects"
)


# Function to plot up and down regulated genes bar plot
plot_updown_bar <- function(df, title_name) {

  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)

  pattern_cols <- grep("^Pattern_", colnames(df), value = TRUE)
  all_conditions <- gsub("^Pattern_", "", pattern_cols)

  pattern_long <- df %>%
    rownames_to_column("Gene") %>%
    pivot_longer(all_of(pattern_cols), names_to = "Condition", values_to = "Status") %>%
    mutate(
      Condition = gsub("^Pattern_", "", Condition),
      Condition = factor(Condition, levels = all_conditions)
    )

  counts <- pattern_long %>%
    filter(Status %in% c("UP", "DOWN")) %>%
    group_by(Condition, Status) %>%
    summarise(N = n(), .groups = "drop") %>%
    complete(
      Condition = all_conditions,
      Status = c("UP", "DOWN"),
      fill = list(N = 0)
    ) %>%
    mutate(
      Condition = factor(Condition, levels = all_conditions),
      label = abs(ifelse(Status == "DOWN", -N, N)),
      N = ifelse(Status == "DOWN", -N, N),
      Condition_label = gsub("^G_D_.*", "G-D", as.character(Condition)),
      Condition_label = gsub("^G_L_.*", "G-LCFA", Condition_label),
      Condition_label = gsub("^G_.*", "G", Condition_label),
      Condition_label = gsub("^L_.*", "LCFA", Condition_label)
    )

  label_map <- counts %>%
    distinct(Condition, Condition_label) %>%
    deframe()

  max_n <- max(abs(counts$N), na.rm = TRUE)
  x_lim <- max_n * 1.20

  counts <- counts %>%
    mutate(
      inside_bar = abs(N) >= 0.12 * max_n,
      x_label = ifelse(
        inside_bar,
        N * 0.5,
        ifelse(N > 0, N + 0.04 * max_n, N - 0.04 * max_n)
      ),
      hjust_lab = ifelse(
        inside_bar,
        0.5,
        ifelse(N > 0, 0, 1)
      ),
      text_col = ifelse(inside_bar, "white", "black")
    )

  ggplot(counts, aes(x = N, y = Condition, fill = Status)) +
    geom_col(width = 0.72, color = "black", linewidth = 0.25) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
    geom_text(
      aes(x = x_label, label = label, hjust = hjust_lab, color = text_col),
      size = 4
    ) +
    scale_y_discrete(labels = label_map) +
    scale_fill_manual(
      values = c("UP" = "firebrick4", "DOWN" = "#1a528b")
    ) +
    scale_color_identity() +
    scale_x_continuous(
      limits = c(-x_lim, x_lim),
      breaks = pretty(c(-max_n, max_n), n = 5),
      labels = abs(pretty(c(-max_n, max_n), n = 5)),
      expand = c(0, 0)
    ) +
    labs(
      x = "Number of genes",
      y = NULL,
      title = title_name
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
      axis.text = element_text(size = 12, colour = "black"),
      axis.title.x = element_text(size = 13),
      axis.title.y = element_blank(),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 11),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
      axis.line = element_blank(),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(10, 20, 10, 20)
    )
}

suplementary_dir <- file.path(output_dir, "001_Figures_paper","Suplementary_Figures")
dir.create(suplementary_dir, showWarnings = FALSE)

for (i in c(1:6)) {

  p <- plot_updown_bar(effects_list[[i]], names(effects_list[i]))

  # Save the bar plot
  pdf(file=file.path(suplementary_dir,paste0("Barplot_DE_genes_effects_",names(effects_list[i]),".pdf")),
      width = 7, height = 5)
  print(p)
  dev.off()
}


#####################
##  Density plots  ##
#####################

# Check that the growth effects objects are in the same order as feature_data
all(rownames(growth_effects)==feature_data$locus_tag)
length(which(rownames(growth_effects)!=feature_data$locus_tag))

th <- 0.05
density_plotter=function(stat,feat=feature_data,threshold=th,culture){
  stat$Type=feat$Type
  stat=stat[which(stat$BH<threshold),]

  pc=stat$LogFC[which(stat$Type=="CDS")]
  npc=stat$LogFC[which(stat$Type=="Stable_RNA")]
  test=wilcox.test(npc,pc,alt="greater")
  print(paste0("Wilcox rank test p=",test$p.value," (alt: npc>pc)"))


  pl=ggplot(stat)+geom_density(aes(x=LogFC,fill=Type),alpha=0.8)+
    xlim(-10,10)+
    xlab(paste0("logFC Response to Growth arrest at ",culture))+
    ylim(0,0.6)+
    theme(legend.position="none")+
    geom_vline(xintercept=0)+
    theme_classic()+
    theme(
      axis.text.y   = element_text(size=14),
      axis.text.x   = element_text(size=14),
      axis.title.y  = element_text(size=14),
      axis.title.x  = element_text(size=14),
      #panel.background = element_blank(),
      #panel.grid.major = element_blank(),
      #panel.grid.minor = element_blank(),
      #axis.line = element_line(colour = "black"),
      panel.border = element_rect(colour = "black", fill=NA, linewidth = 1,linetype="solid"),
      #legend.title=element_blank(),
      legend.position="top",
      #legend.text=element_text(size=14),
      legend.key.size = unit(1, 'lines'))

  return(pl)
}

# Create a vector to loop through the stats
vec_stat <- c(1:4)
carbon_source_name <- c("Glycerol-Dextrose","Glycerol","Glycerol-Lipid","LCFA")

for (i in vec_stat) {

  # We define the indices for each condition
  idx_Fe   <- i
  idx_noFe <- i + 4

  idx_BH_Fe   <- i + 8
  idx_BH_noFe <- i + 12

  effect_name <- carbon_source_name[i]
  print(paste("Processing:", effect_name))

  # Create dataframes for each condition
  stat_noFe <- data.frame(
    LogFC = growth_effects[[idx_noFe]],
    BH    = growth_effects[[idx_BH_noFe]]
  )
  stat_Fe <- data.frame(
    LogFC = growth_effects[[idx_Fe]],
    BH    = growth_effects[[idx_BH_Fe]]
  )

  pl_noFe <- density_plotter(stat_noFe, culture = "Fe/-")
  pl_Fe   <- density_plotter(stat_Fe, culture = "Fe/+")

  combined_pl <- plot_grid(pl_noFe, pl_Fe, nrow = 2, align = "hv")

  pdf(file = file.path(suplementary_dir, paste0("Density_", effect_name, "_Fe_effects.pdf")),
      width = 6, height = 6)
  print(combined_pl)
  dev.off()
}
