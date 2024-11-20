# 2023-12-07
# Loose lab, School of Life Sciences, University of Nottingham
# This script generates a split violin plot detailing the effects of removing 
# alignment errors from the sequencing data during the baseline analysis.


# Set working directory to the specified path
setwd("C:/Users/mbzsjk/My_Freyja/big_tables")

# Install and load required packages
install.packages("hrbrthemes")
library(hrbrthemes)
library(viridis)
library(readxl)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(magrittr)
library(asht)

# Install introdataviz from GitHub for split violin plot
#library(devtools)  # Required for installing packages from GitHub
#devtools::install_github("psyteachr/introdataviz", force = TRUE)
library(introdataviz)

# Define filenames for data files
file1 <- "singletable.exeter_ont_NT002_baseline.tsv"
file2 <- "singletable.exeter_ont_NT002_baseline_align_errors_removed.tsv"
file3 <- "singletable.exeter_ill_NT002_baseline.tsv"
file4 <- "singletable.exeter_ill_NT002_baseline_align_errors_removed.tsv"

# Read data from the files into separate tables
table1 <- read.table(file1, sep="\t", header=TRUE)
table2 <- read.table(file2, sep="\t", header=TRUE)
table3 <- read.table(file3, sep="\t", header=TRUE)
table4 <- read.table(file4, sep="\t", header=TRUE)

# Add a 'Threshold' column to each table to differentiate baseline and optimal data
table1$Threshold <- "Initial"
table2$Threshold <- "Align errors removed"
table3$Threshold <- "Initial"
table4$Threshold <- "Align errors removed"

# Combine all tables into a single DataFrame
combined_table <- rbind(table1, table2, table3, table4)

# Add a 'dataset' column to differentiate ONT and ILLUMINA samples
combined_table <- combined_table %>% 
  mutate(dataset = ifelse(substr(sample, nchar(sample) - 2, nchar(sample)) == "ont",
                          "ONT",
                          ifelse(substr(sample, nchar(sample) - 2, nchar(sample)) == "ill",
                                 "ILL",
                                 NA)))

# Select relevant columns and reorder them
combined_table <- combined_table %>% 
  select(sample, sensitivity, precision, Jaccard_similarity, lineage, mean_cov, proportion, Threshold, dataset)

# Convert the data to a longer format for plotting
combined_table <- combined_table %>% 
  pivot_longer(cols = c(sensitivity, precision, Jaccard_similarity),
               names_to = "parameter",
               values_to = "value") 

# Uppercase the parameter names for consistency
combined_table$parameter <- toupper(combined_table$parameter)

# Define the order of parameters for plotting
parameter_order <- c("SENSITIVITY", "PRECISION", "JACCARD_SIMILARITY")
threshold_order <- c("Initial", "Align errors removed")
combined_table$parameter <- factor(combined_table$parameter, levels = parameter_order)
combined_table$Threshold <- factor(combined_table$Threshold, levels = threshold_order)

# Create the first plot
split_violin <- ggplot(combined_table, aes(x = dataset, y = value, fill = Threshold)) +
  introdataviz::geom_split_violin(alpha = 0.8, trim = TRUE, scale = "width") +
  facet_grid(. ~ parameter, scales = "free") +
  geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(x = NULL, y = NULL)







sens <- wsrTest(table1$sensitivity, table3$sensitivity)
prec <- wsrTest(table1$precision, table3$precision)
jacc <- wsrTest(table1$Jaccard_similarity, table3$Jaccard_similarity)

result <- tibble(
  dataset = "NT002",
  med_ONT_sens = median(table1$sensitivity),
  med_ILL_sens = median(table3$sensitivity),
  wil_sens = sens[["p.value"]],
  med_ONT_prec = median(table1$precision),
  med_ILL_prec = median(table3$precision),
  wil_prec = prec[["p.value"]],
  med_ONT_jacc = median(table1$Jaccard_similarity),
  med_ILL_jacc = median(table3$Jaccard_similarity),
  wil_jacc = jacc[["p.value"]]
)



