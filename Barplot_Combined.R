# Load required libraries
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))

rm(list = ls(all = TRUE))

#set date
current_date <- format(Sys.Date(), "%d%m%Y")

# Set paths
folder_path <- "your/path/heteroplasmy"
samplesheet_path <- "your/path/Data/samplesheet.csv"

# Read the samplesheet
samplesheet <- read.csv(samplesheet_path, header = TRUE, stringsAsFactors = FALSE)

# List all CSV files
csv_files <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)
file_names <- list.files(path = folder_path, pattern = "\\.csv$", full.names = FALSE)

# Initialize list for all mutation records
all_mutation_data <- data.frame()

# Loop through files
for (i in seq_along(csv_files)) {
  file <- csv_files[i]
  file_name <- file_names[i]
  
  # Extract sample ID
  sample_id <- sub("_heteroplasmy$", "", tools::file_path_sans_ext(file_name))
  
  # Match sample to samplesheet
  sample_row <- samplesheet[samplesheet[[1]] == sample_id, ]
  if (nrow(sample_row) == 0) {
    message("Sample ID not found in samplesheet: ", sample_id)
    next
  }
  
  # Define condition
  condition <- as.character(sample_row[[3]])
  condition <- ifelse(condition == "seq10x", "Aneuploid",
                      ifelse(condition == "sub10x", "Euploid", condition))
  
  # Read heteroplasmy data (Heteroplasmy % is stored on the 12th column)
  data <- tryCatch({
    read.csv(file, header = TRUE, skip = 5)
  }, error = function(e) {
    message("Error reading file: ", file, "\n", e)
    return(NULL)
  })
  if (is.null(data)) next
  
  data <- data[!is.na(data[[12]]), ]
  data <- data[data[[12]] > 0, ]
  positions <- str_extract(data[[1]], "(?<=chrM:)\\d+")
  positions <- as.numeric(positions)
  
  if (!is.null(positions) && length(positions) > 0) {
    tmp_df <- data.frame(
      Mutation = positions,
      Condition = condition,
      stringsAsFactors = FALSE
    )
    all_mutation_data <- rbind(all_mutation_data, tmp_df)
  }
}

# Count mutations per condition
mutation_summary <- all_mutation_data %>%
  group_by(Mutation, Condition) %>%
  summarise(Count = n(), .groups = "drop")

# Set factor order for consistent coloring
mutation_summary$Condition <- factor(mutation_summary$Condition, levels = c("Euploid", "Aneuploid"))

# Plot combined bar chart
plot <- ggplot(mutation_summary, aes(x = as.factor(Mutation), y = Count, fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "SNP Frequencies by Condition",
       x = "Mutation (chrM Position)",
       y = "Count",
caption = paste("Generated on", current_date)) +
  theme_minimal() +
 theme(
  axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_rect(fill = "white", color = NA),
  plot.background = element_rect(fill = "white", color = NA),
  legend.background = element_rect(fill = "white", color = NA),
  axis.line = element_line(color = "black"),
  plot.title = element_text(hjust = 0.5)
) +
  scale_fill_manual(values = c("Euploid" = "forestgreen", "Aneuploid" = "darkorange")) +
  guides(fill = guide_legend(title = "Condition"))

# Save plot
ggsave(filename ="your/path/Output/Barplot_Combined.png",
       plot = plot, width = 30, height = 15, dpi = 300)
