# Enrichment Analysis

# Load Data set
install.packages("msigdbr")
library(msigdbr)
library(tidyverse)


# Load Data
setwd("E:/Brain Cancer/GSE_124814")
dir()

# Filter Data
data <- read.table("DEG.txt", sep = "\t")
nrow(data)
View(data)

data_subset <- subset(data, data$adj.P.Val < 0.05)
nrow(data_subset)
View(data_subset)

data_subset_RowNames <- data_subset |>
  rownames_to_column(var = "GeneSymbol")
view(data_subset_RowNames)

write.table(data_subset_RowNames, "adj.P.val_filter.txt", quote = F, sep = "\t")

# C3: mi-RNA targets, Transcription Factor targets

msigdb_data <- msigdbr(species = "Homo sapiens", category = "C3")
view(msigdb_data)


# Cluster / Enrich

library(clusterProfiler)

gene_symbol_Sub <- msigdb_data |>
  select(gs_name, gene_symbol)
view(gene_symbol_Sub)
nrow(gene_symbol_Sub)

gene_symbol_list <- data_subset_RowNames$GeneSymbol
length(gene_symbol_list)


enricher_run <- enricher(gene = gene_symbol_list, TERM2GENE = gene_symbol_Sub)
result <- enricher_run@result
view(result)
nrow(result)

write.table(result, "mi_RNA_Enrichmen_result.txt", quote = F, sep = "\t")



# C5: Gene Ontology(GO)

msigdbr_data_GO <- msigdbr(species = "Homo sapiens", category = "C5")
View(msigdbr_data_GO)

# Cluster / Enrich

library(clusterProfiler)

gene_symbol_Sub_C5 <- msigdbr_data_GO |>
  select(gs_name, gene_symbol)
View(gene_symbol_Sub_C5)
nrow(gene_symbol_Sub_C5)


gene_symbol_list <- data_subset_RowNames$GeneSymbol
length(gene_symbol_list)


enricher_run_C5 <- enricher(gene = gene_symbol_list, TERM2GENE = gene_symbol_Sub_C5)
result_C5 <- enricher_run_C5@result
View(result_C5)
nrow(result_C5)

write.table(result_C5, "GO_Enrichmen_result.txt", quote = F, sep = "\t")



# C2: KEGG

msigdbr_data_KEGG <- msigdbr(species = "Homo sapiens",
                             category = "C2")
View(msigdbr_data_KEGG)

# Cluster / Enrich

library(clusterProfiler)

gene_symbol_Sub_C2 <- msigdbr_data_KEGG |>
  select(gs_name, gene_symbol)
View(gene_symbol_Sub_C2)
nrow(gene_symbol_Sub_C2)


gene_symbol_list <- data_subset_RowNames$GeneSymbol
length(gene_symbol_list)


enricher_run_C2 <- enricher(gene = gene_symbol_list,
                            TERM2GENE = gene_symbol_Sub_C2)
result_C2 <- enricher_run_C2@result
View(result_C2)
nrow(result_C2)

write.table(result_C2, "KEGG_Enrichmen_result.txt", quote = F, sep = "\t")

# Filter KEGG subgroup in data set

kegg_genesets <- result_C2 |>
  filter(grepl("KEGG", ID, ignore.case = TRUE))
view(kegg_genesets)
write.table(kegg_genesets, "KEGG_Enrichmen_result_geneset.txt", quote = F, sep = "\t")

# Plot

library(enrichplot)

# GO Plot

ggsave("GO_dotplot.png", width = 10, height = 8, dpi = 300)
dotplot(enricher_run_C5,
        showCategory = 20,
        x = "GeneRatio",
        color = "p.adjust",
        font.size = 10,
        title = "GO Biological Process Enrichment") +
  labs(
    x = "Gene Ratio (DEGs / Total Genes)",
    y = "GO Terms",
    subtitle = paste0("Top ", 20, " Enriched Pathways"),
    caption = "Analysis by clusterProfiler"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray50"),
    plot.caption = element_text(hjust = 1, size = 10, color = "gray50")
  )
dev.off()



# KEGG Plot

ggsave("KEGG_dotplot.png", width = 10, height = 8, dpi = 300)
dotplot(enricher_run_C2,
        showCategory = 20,
        x = "GeneRatio",
        color = "p.adjust",
        font.size = 8,
        title = "KEGG Biological Process Enrichment") +
  labs(
    x = "Gene Ratio (DEGs / Total Genes)",
    y = "KEGG Terms",
    subtitle = paste0("Top ", 20, " Enriched Pathways"),
    caption = "Analysis by clusterProfiler"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray50"),
    plot.caption = element_text(hjust = 1, size = 10, color = "gray50")
  )
dev.off()

#########
# KEGG Data Set Plot

kegg_dataset <- read.table("KEGG_Enrichmen_result_geneset.txt", header = T)
View(kegg_dataset)

kegg_dataset_new <- kegg_dataset |>
  separate(
    col = GeneRatio,
    into = c("Gene", "Total"),
    convert = TRUE
  ) |>
  mutate(
    GeneRatio = Gene/Total
  )
View(kegg_dataset_new)
?separate

top_20 <- kegg_dataset_new |>
  arrange(p.adjust) |>
  slice(1:20)
View(top_20)


ggsave("KEGG_Plot_2.png", width = 10, height = 8, dpi = 300)
ggplot(top_20,
       aes(x = GeneRatio,
           y = reorder(Description, GeneRatio))) +
  geom_point(aes(size = Count,
                 colour = p.adjust)) +
  scale_color_continuous(low = "red", high = "blue") +
  labs(x = "Gene Ratio",
       y = "KEGG Terms",
       colour = "Adjusted P",
       size = "Count") +
  theme_bw()
dev.off()







