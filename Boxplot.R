# Clear environment and set working directory
rm(list = ls(all = TRUE))
setwd("your/path/Data")
 
# Load packages
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(tibble))

#Set date
current_date <- format(Sys.Date(), "%d-%m-%Y")

# Load data
depth = data.table::fread("depthPerPos.csv", header = TRUE, stringsAsFactors = FALSE) %>% as.data.frame()
sampleSheet = data.table::fread("samplesheet.csv", stringsAsFactors = FALSE) %>% as.data.frame()
pathMITO = read.csv("pathMITO.csv")

# Define original condition list
conditions <- c("sub30x", "sub20x", "sub40x", "sub10x", "seq10x")

# Expand sampleSheet to include both condition columns
sampleSheet_expanded <- sampleSheet %>%
  pivot_longer(cols = c(condition, condition2), names_to = "condition_type", values_to = "condition") %>%
  filter(!is.na(condition))

# Group samples and calculate mean depth per condition
group_conditions <- function(cond) {
  sampleSheet_expanded %>%
    filter(condition == cond) %>%
    pull(sampleID) %>%
    unique()
}

for (cond in conditions) {
  samples <- group_conditions(cond)
  depth[[paste0("mean_", cond)]] <- rowMeans(depth %>% select(all_of(samples)), na.rm = TRUE)
}

# Filter to pathological positions only
indication <- depth %>%
  filter(POS %in% pathMITO$Position) %>%
  select(POS, contains("mean_"))

# Reformat for plotting/statistics
depth.longer <- indication %>%
  pivot_longer(cols = -POS, names_to = "conditionID", values_to = "depth") %>%
  mutate(
    conditionID = gsub("mean_", "", conditionID),
    conditionID = factor(conditionID, levels = conditions),
    depthLog10 = log10(depth)
  )

# Count unique samples per condition across both condition columns
sample_counts <- sampleSheet %>%
  pivot_longer(cols = c(condition, condition2), names_to = "condition_type", values_to = "condition") %>%
  filter(!is.na(condition), condition %in% conditions) %>%
  distinct(sampleID, condition) %>%
  group_by(condition) %>%
  summarise(n = n(), .groups = "drop")

# Define custom display names
condition_labels <- c(
  sub40x = "Mitotic",
  sub30x = "Meiotic I",
  sub20x = "Meiotic II",
  sub10x = "Euploid",
  seq10x = "Aneuploid"
)

# Create final label with name and sample count
sample_counts$condition <- as.character(sample_counts$condition)
sample_counts <- sample_counts %>%
  mutate(display_label = paste0(condition_labels[condition], "\n(n=", n, ")"))

labels.coverage <- setNames(sample_counts$display_label, sample_counts$condition)

# Total number of samples
total_samples <- nrow(sampleSheet)

# Get only conditions with data (non-empty groups)
available_conditions <- unique(na.omit(depth.longer$conditionID))
available_conditions <- as.character(available_conditions)

# Ensure each condition has at least 2 samples
valid_conditions <- sample_counts %>%
  filter(n >= 2) %>%
  pull(condition)

# Generate all pairwise comparisons of valid conditions
comparisons_list <- combn(valid_conditions, 2, simplify = FALSE)


# Create boxplot
p_box <- ggplot(depth.longer, aes(x = conditionID, y = depth, fill = conditionID)) +
  geom_boxplot(outlier.size = 0.5, width = 0.6) +
  stat_compare_means(comparisons = comparisons_list,
                     method = "t.test", label = "p.format", size = 2) +
  scale_x_discrete(labels = labels.coverage) +
  ylab("Reads per Position") +
  xlab("Condition") +
  ggtitle(paste("Pathological Positions - Total Samples:", total_samples, "-", current_date)) +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = 1, size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

# Save plot
ggsave(p_box, file = "your/path/Output/boxplot_path.jpg",
       width = 120, height = 80, units = "mm", dpi = 600)
