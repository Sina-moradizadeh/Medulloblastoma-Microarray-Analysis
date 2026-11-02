# Set Work Directory
setwd("E:/Brain Cancer/GSE_124814")
dir()

#Load Packages
library(tidyverse)

# Load Data
my_data <- read.table("Raw_Data.txt" , header = T)
view(my_data)


data_mean <- my_data |>
  summarise(across(everything(), mean, na.rm = TRUE))
View(data_mean)

range(my_data)
dim(my_data)

# Boxplot samples 200-250
pdf("02_Boxplot.pdf", width = 8, height = 6)
boxplot(my_data[,200:250],
        las = 2,
        main = "Raw Data Expression")
dev.off()


# Histogram
pdf("02_Histogram.pdf", width = 7, height = 6)
hist(my_data[, 500],
     breaks = 100,
     main = "Sample_500"
)

hist(my_data[, 10],
     breaks = 100,
     main = "Sample_10"
)

hist(my_data[, 800],
     breaks = 100,
     main = "Sample_800"
)

hist(my_data[, 1000],
     breaks = 100,
     main = "Sample_1000"
)
dev.off()



# Density Plot
pdf("02_Density.pdf", width = 8, height = 6)
plot(density(my_data[,1]))

for (i in 50:60) {
  lines(density(my_data[, i]),
  col = i + 1,
  lwd = 1.5)
}

dev.off()



#The expression matrix obtained from GEO had already been background corrected
#and normalized using the gRMA algorithm, as reported in original publication.
#Therefore, no additional normalization was applied prior to downstream analysis.




