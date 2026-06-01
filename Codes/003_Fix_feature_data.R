#############################
### 0. Load dependencies. ###
#############################

library(dplyr)
library(tibble)
library(writexl)

###################################################
### 2. Set working directory and create folders ###
###################################################

main_wd <- getwd()
setwd(main_wd)
input_dir <- "Inputs/002_Processed_data"
output_dir <- "Outputs"

###########################
####### 3. Load data ######
###########################
reads        <- read.table(file.path(input_dir,"txt/reads_32_samples.txt"),check.names = F)
reads_nc     <- read.table(file.path(input_dir,"txt/reads_nc.txt"),check.names = F)
metadata     <- read.table(file.path(input_dir,"txt/metadata_32_samples.txt"),check.names = F)
feature_data <- read.table(file.path(input_dir,"txt/feature_data_genes.txt"),check.names = F)

### Release 5 (2024-07-11) of the M. tuberculosis H37Rv genome annotation.
mtb_meta     <- read.delim(file.path(input_dir,"txt/Mycobacterium_tuberculosis_H37Rv_txt_v5.txt"),check.names = F)

### Create feature rows for custom ncRNA loci absent from the base annotation.
new_rows <- tibble(locus_tag = rownames(reads_nc), class="ncRNA_custom", symbol=NA, product=NA, genome=NA,
                   description=NA, gene_type="ncRNA", gene_name=NA,
                   gene_id=NA, transcript_id=NA, GO_terms=NA,
                   KEGG_pathways=NA, COG_category=NA)

# Add the new columns to the new_rows tibble
for (col in setdiff(names(feature_data), c("locus_tag","class","gene_type"))) {
  new_rows[[col]] <- NA
}

# Bind the new rows to the existing feature_data
feature_data <- bind_rows(feature_data, new_rows)

rownames(feature_data) <- feature_data$locus_tag
feature_data <- feature_data[rownames(reads),]

stopifnot(identical(rownames(feature_data), rownames(reads)))
stopifnot(identical(colnames(reads), rownames(metadata)))

unique(feature_data$class)
# "ncRNA_custom"   "protein_coding" "tRNA" "pseudogene"     "ncRNA"
# "rRNA"           "misc_RNA"
#  "other"


### Create directory to save gene names by class.
gene_names_by_class_dir <- file.path(output_dir, "gene_names_by_class")
dir.create(gene_names_by_class_dir, showWarnings = FALSE)

### Save gene names by class.
for (gene_class in unique(feature_data$class)) {
  gene_names <- feature_data %>%
    filter(class == gene_class) %>%
    pull(locus_tag)

  write.table(gene_names, file = file.path(gene_names_by_class_dir, paste0("gene_names_", gene_class, ".txt")),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  write_xlsx(as.data.frame(gene_names), path = file.path(gene_names_by_class_dir, paste0("gene_names_", gene_class, ".xlsx")))
}


### Create pseudogene data frame.
pseudogene_data <- feature_data %>%
  filter(class == "pseudogene")

### Search for pseudogene information in the MTB metadata.
mtb_pseudogene_data <- mtb_meta[mtb_meta$Locus %in% pseudogene_data$locus_tag, ]
mtb_pseudogene_data <- mtb_pseudogene_data[order(mtb_pseudogene_data$Locus), ]

### Check for pseudogenes not found in the MTB metadata.
missing <- pseudogene_data[!pseudogene_data$locus_tag %in% mtb_meta$Locus , ]
# Rv1887a this is the one missing pseudogene in the MTB metadata

### Remove missing pseudogenes from the pseudogene data.
pseudogene_data <- pseudogene_data[-which(pseudogene_data$locus_tag %in% missing$locus_tag), ]

### Ensure the order matches.
stopifnot(identical(mtb_pseudogene_data$Locus, pseudogene_data$locus_tag))

### Update pseudogene data with MTB metadata.
pseudogene_data$New_Feature   <- mtb_pseudogene_data$Feature
pseudogene_data$New_Name      <- mtb_pseudogene_data$Name
pseudogene_data$New_Function  <- mtb_pseudogene_data$Function
pseudogene_data$Is_Pseudogene <- mtb_pseudogene_data$Is_Pseudogene
pseudogene_data$New_Functional_Category <- mtb_pseudogene_data$Functional_Category

### Save the updated pseudogene data.
write.table(pseudogene_data, file.path(input_dir, "txt/pseudogene_data_updated.txt"), sep = "\t")
write_xlsx (pseudogene_data, path = file.path(input_dir, "xlsx/pseudogene_data_updated.xlsx"))

### Check the whole feature data against the MTB metadata.

unique(mtb_meta$Feature)
# "ncRNA"      "tRNA"       "CDS"        "rRNA"       "promoter"   "misc_RNA"
#  "-10_signal" "-35_signal"

### Check how many features are in each class.
class_summary <- feature_data %>%
  group_by(class) %>%
  summarise(count = n())

# class             count
# <chr>             <int>
# 1 misc_RNA          2
# 2 ncRNA             20
# 3 ncRNA_custom     183
# 4 other              2
# 5 protein_coding  3906
# 6 pseudogene        30
# 7 rRNA               3
# 8 tRNA              45

mtb_meta <- mtb_meta[,-35] # Remove the empty repeated column
class_summary_mtb <- mtb_meta %>%
  group_by(Feature) %>%
  summarise(count = n())

# Feature        count
# <chr>           <int>
# 1 -10_signal      4
# 2 -35_signal      4
# 3 CDS            4031
# 4 misc_RNA        2
# 5 ncRNA          92
# 6 promoter        6
# 7 rRNA            3
# 8 tRNA           45

### Add the Type column used downstream to separate coding and stable RNA features.
feature_data <- feature_data %>%
  mutate(Type = case_when(
    class == "protein_coding" ~ "CDS",
    class == "pseudogene"     ~ "CDS",
    class == "other"          ~ "CDS",
    class == "ncRNA_custom" ~ "Stable_RNA",
    class == "ncRNA"        ~ "Stable_RNA",
    class == "misc_RNA"     ~ "Stable_RNA",
    class == "tRNA"         ~ "Stable_RNA",
    class == "rRNA"         ~ "Stable_RNA")) %>%
  relocate(Type, .after = class)

unique(feature_data$Type)

length(which(feature_data$Type == "CDS"))          # 3938
length(which(feature_data$Type == "Stable_RNA"))   # 253

write.table(feature_data, file.path(input_dir,"txt/feature_data_genes_type_added.txt"), sep = "\t")
