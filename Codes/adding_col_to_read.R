
x <- metadata$Full_path
batch <- gsub("^Mtb_Latency/|/RawFASTQ.*$", "", x)
unique_batches <- unique(batch)

# [1] "HN00153515_C1_C2_Julio_2021"       "HN00166209_C11_C12_C2r_Marzo_2022" "HN00150037_C5_C6_Junio_2021"
# [4] "HN00161701_C7_C8_Diciembre_2021"   "HN00170106_C14_C15_C8r_Mayo_2022"  "HN00189620_C11_C15_Abril_2023"

feature_data <- read.table(file.path(input_dir,"txt/feature_data_genes_type_added.txt"))


length(which(rownames(feature_data) != rownames(reads_32)))

reads_32$Biotype <- feature_data$Type
reads_32$low_expr_filter <- "NO"
reads_32$low_expr_filter[rownames(reads_32) %in% rownames(feature_data)] <- "YES"

write.xlsx(as.data.frame(reads_32), file = file.path(input_dir, "xlsx/reads_32_with_biotype.xlsx"))
