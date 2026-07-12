# Set Work Directory
setwd("E:/Brain Cancer/GSE_124814")
dir()

# Load Packages
library(tidyverse)
library(pheatmap)

# Load Data
deg_fil <- read.delim("Filtered_DEG.txt", header = T)
View(deg_fil)

meta_data <- read.delim("MetaData.txt", header = T)
view(meta_data)

my_data <- read.table("Raw_Data.txt" , header = T)

deg_fil_2 <- deg_fil |>
  mutate(
    Significance = case_when(
      adj.P.Val < 0.05 & logFC > 1 ~ "Up-regulated",
      adj.P.Val < 0.05 & logFC < -1 ~ "Down-regulated",
      TRUE ~ "Not significant"
    )
  )

view(deg_fil_2)
write.csv(deg_fil_2, "DEG_Fil_2.csv", quote = F)

# Volcano Plot

pdf("05_Volcano_Plot.pdf", width = 8, height = 6)
deg_fil <- deg_fil |>
  mutate(
    Significance = case_when(
      adj.P.Val < 0.05 & logFC > 1 ~ "Up-regulated",
      adj.P.Val < 0.05 & logFC < -1 ~ "Down-regulated",
      TRUE ~ "Not significant"
    )
  )


ggplot(deg_fil, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
  geom_point(size = 1.5, alpha = 0.7) +
  scale_color_manual(
    values = c("Up-regulated" = "#E74C3C", 
               "Down-regulated" = "#2E86C1", 
               "Not significant" = "gray80")
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.5) +
  theme_minimal() +
  labs(
    title = "Volcano Plot",
    subtitle = "Differential Expression Analysis",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = "Regulation"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.position = "top",
    panel.grid.minor = element_blank()
  )
dev.off()


# Heatmap

pdf("05_pheatmap.pdf", height = 8, width = 6)


deg_fil <- deg_fil |>
  mutate(
    Regulation = case_when(
      logFC > 0 ~ "Up-regulated",
      logFC < 0 ~ "Down-regulated",
      TRUE ~ "Not significant"
    )
  )

deg_top <- deg_fil |>
  arrange(desc(abs(logFC))) |>
  head(100)

top_genes <- rownames(deg_top)
heatmap_data <- my_data[top_genes, ]


annotation_row <- data.frame(
  Regulation = deg_top$Regulation
)
rownames(annotation_row) <- rownames(deg_top)


annotation_col <- data.frame(
  Group = meta_data$Sample_type
)
rownames(annotation_col) <- colnames(heatmap_data)


annotation_colors <- list(
  Regulation = c("Up-regulated" = "#E74C3C", "Down-regulated" = "#2E86C1"),
  Group = c("Normal" = "#2ECC71", "Medulloblastoma" = "#E74C3C")
)


pheatmap(heatmap_data,
         scale = "row",
         clustering_method = "ward.D2",
         show_rownames = FALSE,
         show_colnames = FALSE,
         annotation_col = annotation_col,
         annotation_row = annotation_row,
         annotation_colors = annotation_colors,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Heatmap of Top 100 DEGs",
         fontsize = 8,
         border_color = NA)

dev.off()


















