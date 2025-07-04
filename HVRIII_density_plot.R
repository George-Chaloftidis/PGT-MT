###########################################################################################################################
# Author: George Chaloftidis
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# Script purpose: Visualize heteroplasmy percentages in the HVRIII region (chrM:438-574) for euploid and aneuploid samples
#                using density plots. Also performs a Kolmogorov–Smirnov (KS) test to assess statistical differences 
#                in heteroplasmy distributions between the two groups.

# Input:
# - Samplesheet: /your/path/Data/samplesheet.csv
# - Per-sample heteroplasmy files: /your/path/Output/heteroplasmy/<sampleID>_heteroplasmy.csv

# Output:
# - Density plot with KS p-value annotation:
#   * your/path/Output/Density_plot_HVRIII_by_group.png

###########################################################################################################################

rm(list=ls(all=T))
library(data.table)
library(ggplot2)
 
#set date
current_date <- format(Sys.Date(), "%d%m%Y")
# Paths
folder_path <- "/your/path/Output/heteroplasmy"
samplesheet_path <- "/your/path/Data/samplesheet.csv"
# Load sample sheet
samplesheet <- fread(samplesheet_path)
samplesheet <- samplesheet[condition %in% c("sub10x", "seq10x"), .(sampleID, condition)]
samplesheet[, Group := fifelse(condition == "sub10x", "Euploid",
fifelse(condition == "seq10x", "Aneuploid", NA_character_))]
# List all CSV files
csv_files <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)
# Initialize list to store heteroplasmy data
heteroplasmy_data <- list()
# Define position range
position_min <- 438
position_max <- 574
# Loop through files
for (file in csv_files) {
# Extract sampleID from filename (remove '_heteroplasmy.csv')
filename_base <- basename(file)
sample_id <- sub("_heteroplasmy\\.csv$", "", filename_base)
# Match with samplesheet
match_row <- samplesheet[sampleID == sample_id]
if (nrow(match_row) == 0) {
warning(paste("Sample", sample_id, "not found in samplesheet or doesn't meet condition filter. Skipping."))
next
}
group <- match_row$Group
# Load data
data <- fread(file, na.strings = c("NA", ""))
if (!"heteroplasmy" %in% names(data) || !"position" %in% names(data)) {
warning(paste("No 'heteroplasmy' or 'position' column in", file))
next
}
# Extract position values (e.g., "chrM:1_G/<NON_REF>" => 1)
data[, position_numeric := as.numeric(sub(".*:(\\d+)_.*", "\\1", position))]
# Filter positions within the range
data <- data[position_numeric >= position_min & position_numeric <= position_max]
if (nrow(data) == 0) {
warning(paste("No data left after filtering positions for sample", sample_id))
next
}
# Extract and clean heteroplasmy values
h_raw <- data[["heteroplasmy"]]
h_values <- suppressWarnings(as.numeric(as.character(h_raw)))
h_values <- h_values[!is.na(h_values)]
h_values <- pmin(h_values, 80)
if (length(h_values) > 0) {
heteroplasmy_data[[length(heteroplasmy_data) + 1]] <- data.frame(
heteroplasmy = h_values,
Group = group
)
}
}
# Combine all data
if (length(heteroplasmy_data) == 0) {
stop("No valid numeric heteroplasmy values found.")
}
df <- rbindlist(heteroplasmy_data)# KS Test
if (all(c("Euploid", "Aneuploid") %in% unique(df$Group))) {
ks_result <- ks.test(
df[Group == "Euploid", heteroplasmy],
df[Group == "Aneuploid", heteroplasmy]
)
ks_p_value <- ks_result$p.value
ks_p_value_text <- paste0("KS test p-value = ", signif(ks_p_value, 3))
} else {
ks_p_value_text <- "KS test: insufficient data"
}

# Plot
plot <- ggplot(df, aes(x = heteroplasmy, fill = Group)) +
geom_density(alpha = 0.5) +
labs(title = "Density Plot of Heteroplasmy Percentages by Group",
x = "Heteroplasmy (%)",
y = "Density",
caption = paste("Generated on", current_date)) +
theme_minimal() +
scale_fill_manual(values = c("Euploid" = "steelblue", "Aneuploid" = "darkred")) +
theme(
panel.background = element_rect(fill = "white", color = NA),
plot.background = element_rect(fill = "white", color = NA),
legend.background = element_rect(fill = "white", color = NA),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank()
) +
annotate("text", x = 10, y = 0.04, label = ks_p_value_text, hjust = 0, size = 5)

ggsave(filename = "your/path/Output/Density_plot_HVRIII_by_group.png",
       plot = plot, width = 20, height = 10, dpi = 300)
