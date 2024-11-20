setwd("C:/Users/mbzsjk/My_Freyja/big_tables")

# Install and load required packages
#install.packages("hrbrthemes")
library(hrbrthemes)
library(viridis)
library(readxl)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(magrittr)
library(asht)
library(data.table)

# Install introdataviz from GitHub for split violin plot
#library(devtools)  # Required for installing packages from GitHub
#devtools::install_github("psyteachr/introdataviz", force = TRUE)
library(introdataviz)


# Read in big_tables
initial_table <- read.table("C:/Users/mbzsjk/My_Freyja/big_tables/baseline/Initial_analysis/bigtable.exeter_ill_vs_exeter_ont_NT002.mincov1.mpileup2snp.nostrandbiasfilter.varscan.tsv") %>% 
  mutate(tech = str_sub(sample, -3, )) %>% 
  mutate(threshold = "Initial") %>% 
  select(sample, sensitivity, precision, Jaccard_similarity, threshold, tech)

align_errors_rem_table <- read.table("C:/Users/mbzsjk/My_Freyja/big_tables/baseline/Align_errors_removed/bigtable.exeter_ill_vs_exeter_ont_NT002_baseline_nostrandbiasfilter_align_errors_removed.tsv") %>% 
  mutate(tech = str_sub(sample, -3, ))%>% 
  mutate(threshold = "Align errors removed")%>% 
  select(sample, sensitivity, precision, Jaccard_similarity, threshold, tech)

# Filter to single tech only
initial_ont <- initial_table %>% 
  filter(tech == "ont")

initial_ill <- initial_table %>% 
  filter(tech == "ill")

align_errors_ont <- align_errors_rem_table %>% 
  filter(tech == "ont")

align_errors_ill <- align_errors_rem_table %>% 
  filter(tech == "ill")


# Bind data into a single table
combined_table <- rbind(initial_ont, initial_ill, align_errors_ont, align_errors_ill)

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
combined_table$threshold <- factor(combined_table$threshold, levels = threshold_order)

# Create the first plot
split_violin <- ggplot(combined_table, aes(x = tech, y = value, fill = threshold)) +
  introdataviz::geom_split_violin(alpha = 0.8, trim = TRUE, scale = "width") +
  facet_grid(. ~ parameter, scales = "free") +
  geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(x = NULL, y = NULL) 






# Attempt to add p-values to graph.

split_violin + stat_pvalue_manual(test2, 
                                  label = "p", 
                                  step.increase = 0.05, 
                                  tip.length = 0.02,
                                  bracket.shorten = 0)




metric <- c("sensitivity", "precision", "Jaccard_similarity")
stat_pos <- c(1.1, 0.60, 0.50)

p_val_summ <- function(metric, pos) {
  
  wilcoxon <- wsrTest(initial_ont[[metric]], initial_ill[[metric]])
  
  result <- tibble(
    metric = metric,
    group1 = "initial_ont",
    group2 = "initial_ill",
    p = wilcoxon[["p.value"]],
    y.position = pos,
    threshold = "Initial"
    
  )
  
}

test <- mapply(p_val_summ, metric, stat_pos)
test2 <- as.data.frame(t(test)) %>% 
  rownames_to_column() %>% 
  select(-c(rowname))

