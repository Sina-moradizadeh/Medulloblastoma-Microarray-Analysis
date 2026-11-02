# Set Work Directory
setwd("E:/Brain Cancer/GSE_124814")
dir()

#Load Packages

library(GEOquery)
library(tidyverse)

# Download With GEOquery
#raw_data <- getGEOSuppFiles("GSE124814")

# Load Data
raw_data <- read.table("GSE124814_HW_expr_matrix.tsv", header = T)
View(raw_data)

raw_data_1 <- raw_data |>
  column_to_rownames(var = "Gene_Symbol")
view(raw_data_1)

nrow(raw_data_1)
# n: 14883
ncol(raw_data_1)
# n: 1641

write.table(raw_data_1, "Raw_Data.txt", sep = "\t", quote = F)







