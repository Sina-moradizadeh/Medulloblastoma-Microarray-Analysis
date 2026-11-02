# Set Work Directory
setwd("E:/Brain Cancer/GSE_124814")
dir()

#Load Packages
library(tidyverse)
library(readxl)


# Load Data
des_table <- read_xlsx("GSE124814_sample_descriptions.xlsx",skip = 1)
view(des_table)

my_data <- read.table("Raw_Data.txt" , header = T)

identical(
  colnames(my_data),
  des_table$`Sample name`
)

# Meta Data
meta_data <- des_table |>
  select(1, 4) |>
  rename(Sample_name = 'Sample name', Sample_type = 'source name') |>
  column_to_rownames(var = "Sample_name")

table(meta_data$Sample_type)
view(meta_data)

write.table(meta_data, "MetaData.txt", quote = F, sep = "\t")

data_N <- meta_data |>
  filter(Sample_type == "Normal")
write.table(data_N, "Sample_N.txt", quote = F, sep = "\t")

data_M <- meta_data |>
  filter(Sample_type == "Medulloblastoma")
write.table(data_M, "Sample_M.txt", quote = F, sep = "\t")

# PCA

pdf("PCA_03.pdf", width = 8, height = 6)
# 1
gene_vars <- apply(my_data, 1, var, na.rm = TRUE)
top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:2000]

pca_result <- prcomp(t(my_data[top_genes, ]), center = TRUE, scale. = TRUE)

# 2
pca_scores <- as.data.frame(pca_result$x)
pca_scores$Group <- meta_data$Sample_type

# 3
explained_var <- summary(pca_result)$importance[2, ] * 100

# 4
ggplot(pca_scores, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 2, alpha = 0.7) +
  stat_ellipse(aes(group = Group), level = 0.95, linetype = "dashed") +
  scale_color_manual(values = c("Normal" = "#2ECC71", "Medulloblastoma" = "#E74C3C")) +
  theme_minimal() +
  labs(
    title = "PCA: Normal vs Medulloblastoma",
    x = paste0("PC1 (", round(explained_var[1], 2), "%)"),
    y = paste0("PC2 (", round(explained_var[2], 2), "%)"),
    color = "Group"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position = "top"
  )

dev.off()

