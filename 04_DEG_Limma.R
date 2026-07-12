# Set Work Directory
setwd("E:/Brain Cancer/GSE_124814")
dir()

#Load Packages
library(tidyverse)

# Load Data
my_data <- read.table("Raw_Data.txt" , header = T)
meta_data <- read.table("MetaData.txt" , header = T)
type(my_data)

# Limma

library(limma)

# Step 1: Convert to Matrix

step_1 <- my_data |>
  as.matrix()
class(step_1)

# Step 2: Make group

step_2 <- factor(meta_data$Sample_type)
table(step_2)

# Step 3: Design matrix

design <- model.matrix(~0 + step_2)
colnames(design) <- levels(step_2)
design

# Step 4: Model

fit <- lmFit(step_1, design)

# Step 5: Contrast

contrast.matrix <- makeContrasts(
  Medulloblastoma - Normal,
  levels = design
)

# Step 6: Execute contrast

fit2 <- contrasts.fit(fit, contrast.matrix)

# Step 7: Empirical Bayes

fit2 <- eBayes(fit2)

# Step 8: DEG table

deg <- topTable(
  fit2,
  adjust.method = "BH",
  number = Inf
)
view(deg)
write.table(deg, "DEG.txt", quote = F, sep = "\t")


deg_fil <- deg[
  deg$adj.P.Val < 0.05 &
    abs(deg$logFC) > 1,
]
view(deg_fil)
write.table(deg_fil, "Filtered_DEG.txt", quote = F, sep = "\t")

















