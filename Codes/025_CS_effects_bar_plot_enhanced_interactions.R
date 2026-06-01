
######################
### 0.Dependencies ###
######################
{
  library(tidyverse)
  library(shadowtext)
  library(cowplot)
  library(writexl)
}

###################################################
### 2. Set working directory and create folders ###
###################################################
main_wd <- getwd()
setwd(main_wd)
input_dir  <- "Inputs/002_Processed_data"
output_dir <- "Outputs"

safe_path_name <- function(x) {
  x <- gsub("->", "to", x, fixed = TRUE)
  x <- gsub("[<>:\"/\\\\|?*]", "_", x)
  x <- gsub("\\s+", "_", x)
  x
}

###########################
####### 3. Load data ######
###########################
metadata         <- read.table(file.path(input_dir,"txt/metadata_32_samples.txt"),check.names = F)
feature_data     <- read.table(file.path(input_dir,"txt/feature_data_genes_type_added.txt"),check.names = F)
df_name_contrast <- read.table(file.path(input_dir,"txt/contrasts_nomenclature.txt"))

### Lets filter the metadata to retrieve only C7-C8 and C11-C12
metadata_filtered <- metadata %>%
  filter((Culture %in% c("C5","C6","C7","C8")) | (Culture %in% c("C11","C12")))

# Get the batch info
batch <- metadata_filtered$Full_path
batch <- gsub("Mtb_Latency/","",batch)
batch <- gsub("/.*","",batch)

unique(batch)

# "HN00150037_C5_C6_Junio_2021"       "HN00161701_C7_C8_Diciembre_2021"   "HN00170106_C14_C15_C8r_Mayo_2022"
# "HN00189620_C11_C15_Abril_2023"     "HN00166209_C11_C12_C2r_Marzo_2022"

# Get the contrasts names of interest and its file names
file_names <- paste0(safe_path_name(df_name_contrast$title),".txt")
# file_names <- file_names[c(9:42,50:51)]
current_folder <- file.path(output_dir,"Data/txt/th_0.05_th_size_0")
### Get full names including folder path
list.of.files = list.files(current_folder, full.names = TRUE)
### Keep only the basename (file names) matching dataframe column
clean_list <- list.of.files[basename(list.of.files) %in% file_names]
stopifnot(length(clean_list) == length(file_names))
data_name <- basename(clean_list)
data_name <- substr(data_name,1,nchar(data_name)-4)

### Read the data
myData <- lapply(clean_list, read.table)
myData <- setNames(object = myData, data_name)

feature_data_filtered <- feature_data[rownames(myData[[1]]),]
length(which(rownames(feature_data_filtered)!= rownames(myData[[1]])))

# phase_effects <- myData[15:22]
# cs_effect     <- myData[c(1:10)]

names(myData)
order_desired <- safe_path_name(df_name_contrast$title)
myData <- myData[order_desired]
stopifnot(!any(vapply(myData, is.null, logical(1))))

# Define and index and a name to save data
configs <- list(
  c(9, 12, 38),     # Config 1: Glycerol_Dextrose comparison, with Fe+
  c(10, 12, 50),     # Config 2 (Glycerol comparison, with Fe+)
  c(11, 12, 37),       # Config 3 (Glycerol_LCFA comparison, with Fe+)

  c(13, 12, 42),     # Config 1: Glycerol_Dextrose comparison, without Fe-
  c(14, 12, 51),     # Config 2 (Glycerol comparison, without Fe-)
  c(15, 12, 41)       # Config 3 (Glycerol_LCFA comparison, without Fe-)
)

# Nombres para tus comparaciones
names_vec <- c(
  "Glycerol_Dextrose_vs_LCFA_with_Fe",
  "Glycerol_vs_LCFA_with_Fe",
  "Glycerol_LCFA_vs_LCFA_with_Fe",

  "Glycerol_Dextrose_vs_LCFA_without_Fe",
  "Glycerol_vs_LCFA_without_Fe",
  "Glycerol_LCFA_vs_LCFA_without_Fe"
)

int_summary   <- list()
fraction_list <- list()

for (indice in seq_along(configs)) {

  idx <- configs[[indice]]
  with_i    <- idx[1]
  without_i <- idx[2]
  int_i     <- idx[3]

  name_i <- names_vec[indice]

  ### ---- CREATE SUMMARY DF ---- ###

  logFC_with    <- data.frame(logFC_with    = myData[[with_i]][,"log2FoldChange_shrunken"])
  logFC_without <- data.frame(logFC_without = myData[[without_i]][,"log2FoldChange_shrunken"])
  logFC_int     <- data.frame(logFC_int     = myData[[int_i]][,"log2FoldChange_shrunken"])

  BH_with    <- data.frame(BH_with    = myData[[with_i]][,"BH"])
  BH_without <- data.frame(BH_without = myData[[without_i]][,"BH"])
  BH_int     <- data.frame(BH_int     = myData[[int_i]][,"BH"])

  summary <- cbind(logFC_with, logFC_without, logFC_int,
                   BH_with, BH_without, BH_int)
  rownames(summary) <- rownames(myData[[with_i]])

  ### ---- THRESHOLDS ---- ###
  threshold     <- 0.05
  threshold_int <- 0.05

  ### ---- LABELING GENES ---- ###
  summary$label_int <- "Background"

  # DOWN
  summary$label_int[summary$logFC_with<0 & summary$logFC_without<0] <- "DOWN"

  # Non-coherent
  summary$label_int[
    summary$logFC_with>0 & summary$logFC_without<0 & summary$BH_int<threshold_int |
      summary$logFC_with<0 & summary$logFC_without>0 & summary$BH_int<threshold_int
  ] <- "Non_Coherent"

  # Enhanced DOWN
  summary$label_int[
    summary$logFC_with<0 & summary$logFC_without<0 &
      summary$BH_int<threshold_int & summary$logFC_int<0
  ] <- "Enhanced_DOWN"

  # Damped DOWN
  summary$label_int[
    summary$logFC_with<0 & summary$logFC_without<0 &
      summary$BH_int<threshold_int & summary$logFC_int>0
  ] <- "Damped_DOWN"

  # UP
  summary$label_int[
    summary$logFC_with>0 & summary$logFC_without>0
  ] <- "UP"

  # Enhanced UP
  summary$label_int[
    summary$logFC_with>0 & summary$logFC_without>0 &
      summary$BH_int<threshold_int & summary$logFC_int>0
  ] <- "Enhanced_UP"

  # Damped UP
  summary$label_int[
    summary$logFC_with>0 & summary$logFC_without>0 &
      summary$BH_int<threshold_int & summary$logFC_int<0
  ] <- "Damped_UP"

  ### ---- SUMMARY TABLE ---- ###
  summary_df <- data.frame(
    summary = summary(factor(summary$label_int,
                             levels = c("DOWN","UP",
                                        "Damped_DOWN","Damped_UP",
                                        "Enhanced_DOWN","Enhanced_UP",
                                        "Non_Coherent","Background")))
  )

  # Add coherent class
  summary_df["Coherent",] <- sum(summary_df[3:6,])

  ### ---- FRACTIONS ---- ###
  fraction <- numeric(6)

  fraction[1] <- 100*nrow(summary[summary$label_int=="Non_Coherent",]) / sum(summary_df[3:7,])
  fraction[2] <- 100*summary_df["Enhanced_DOWN",] / sum(summary_df[3:7,])
  fraction[3] <- 100*summary_df["Enhanced_UP",]   / sum(summary_df[3:7,])
  fraction[4] <- 100*summary_df["Damped_DOWN",]   / sum(summary_df[3:7,])
  fraction[5] <- 100*summary_df["Damped_UP",]     / sum(summary_df[3:7,])
  fraction[6] <- sum(fraction[2:5])

  fraction_df <- data.frame(
    Description = c("non_coh_vs_coh","Enhanced_DOWN","Enhanced_UP",
                    "Damped_DOWN","Damped_UP","total_coherent_fraction"),
    Fraction = fraction
  )

  ### ---- SAVE OUTPUTS ---- ###
  write.table(
    fraction_df,
    file.path(input_dir, paste0("txt/fraction_genes_carbon_source_effects_", name_i, ".txt")),
    quote = FALSE, sep = "\t"
  )

  # save summary_df or summary if needed
  write.table(summary_df,
              file.path(input_dir, paste0("txt/summary_table_carbon_source_effects_", name_i, ".txt")))

  cat("\n Finished configuration:", name_i, "\n")

  int_summary[[indice]]   <- summary_df
  fraction_list[[indice]] <- fraction_df
}

# Combine all summaries into a final data frame
final_summary <- data.frame(
  int_summary[[1]],
  int_summary[[2]],
  int_summary[[3]],
  int_summary[[4]],
  int_summary[[5]],
  int_summary[[6]]
)

# Give proper column names
colnames(final_summary) <- names_vec

# Save the final summary
write.table(
  final_summary,
  file.path(input_dir, "txt/final_summary_carbon_source_effects.txt"))

# Assign names to fraction_list
names(fraction_list) <- names_vec

# Save the fraction list

write.table(fraction_list,
            file.path(input_dir, "txt/fraction_list_carbon_source_effects.txt"))

names(fraction_list) <- make.names(substr(names(fraction_list), 1, 31), unique = TRUE)
write_xlsx(
  fraction_list,
  path = file.path(input_dir, "xlsx/fraction_list_carbon_source_effects.xlsx")
)

#################################################################
##                       stacked barplot                       ##
#################################################################

categorie <- factor(rownames(final_summary[3:7,]),levels = c("Non_Coherent","Damped_DOWN","Damped_UP",
                                                          "Enhanced_DOWN","Enhanced_UP"))
categorie <- rep(categorie,6)

values       <- c(final_summary[3:7,1],
                  final_summary[3:7,2],
                  final_summary[3:7,3],
                  final_summary[3:7,4],
                  final_summary[3:7,5],
                  final_summary[3:7,6])
interactions <- names_vec
interactions <- rep(colnames(final_summary),each=5)
interactions <- factor(interactions,levels = unique(interactions))

df  <- data.frame(interactions,categorie,values)
breaks_values <- c(pretty(0:max(final_summary)))

# breaks_values <- pretty(0:3000)
color <- c("grey","#A59FCB","#FFEFD9","#50486D","#FFA373")

bar_plt <- ggplot(df, aes(fill=categorie, y=values, x=interactions)) +
  geom_bar(position="stack", stat="identity",width = 0.7, color = "black")+
  coord_flip()+
  scale_y_continuous(expand=c(0,0),breaks = breaks_values,labels = abs(breaks_values),limits=c(0,3000))+
  scale_fill_manual(values = color)+
  theme_classic()+
  geom_hline(yintercept=0)+
  # theme(aspect.ratio = .9)+
  labs(x="Interactions", y = "Number of Genes")+
  theme(
    axis.text.x   = element_text(size=12),
    axis.text.y   = element_text(size=12),
    axis.title.x  = element_text(size=12),
    axis.title.y  = element_text(size=12),
    axis.ticks.y = element_blank(),
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
    aes(x=interactions, y=values, label = values),
    hjust = 2,
    position = "stack",
    # nudge_x = -0.3,
    colour = "black",
    bg.colour = "white",
    bg.r = 0.2,
    # family = "Econ Sans Cnd",
    size = 3)

print(bar_pl_1)

pdf(file = file.path(output_dir,"001_Figures_paper","CS_effects_Stacked_Bar_plot_interactions_enhanced.pdf"),width=8,height=5)
print(bar_pl_1)
dev.off()
