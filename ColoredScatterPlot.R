###########################################################################################################################
# Author: George Chaloftidis
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# Script purpose: Visualization of heteroplasmic SNP distribution across mitochondrial positions per sample category

# Input: Per-sample heteroplasmy files (.csv), samplesheet (.csv)

# Output: Scatter plots showing heteroplasmy (%) by mitochondrial position for Euploid and Aneuploid samples (.png)

###########################################################################################################################

# Load required libraries
library(ggplot2)
library(dplyr)
library(stringr)

rm(list = ls(all = TRUE))

#set date
current_date <- format(Sys.Date(), "%d%m%Y")


# Set paths
folder_path <- "your/path/Output/heteroplasmy"
output_path <- "your/path/Output"
samplesheet_path <- file.path("your/path/Data/samplesheet.csv")

# Read the sample sheet
samplesheet <- read.csv(samplesheet_path, header = TRUE, stringsAsFactors = FALSE)

# Ensure it has at least 3 columns
if (ncol(samplesheet) < 3) stop("Samplesheet must have at least 3 columns.")

# Construct expected file paths using the pattern: sampleID_heteroplasmy.csv
samplesheet$FullPath <- file.path(folder_path, paste0(samplesheet[[1]], "_heteroplasmy.csv"))

# Filter to existing files
samplesheet <- samplesheet[file.exists(samplesheet$FullPath), ]

# Map categories
samplesheet$Category <- ifelse(samplesheet[[3]] == "sub10x", "Euploid",
                               ifelse(samplesheet[[3]] == "seq10x", "Aneuploid", NA))

# Remove rows with unknown categories
samplesheet <- samplesheet[!is.na(samplesheet$Category), ]

# Function to process and plot by category
process_category <- function(category_name, file_paths) {
  scatter_data <- data.frame(Position = numeric(), Value = numeric())
  
  for (file in file_paths) {
    cat("Reading:", file, "\n")
    
    tryCatch({
      data <- read.csv(file, header = TRUE)
      
      if (ncol(data) < 12) {
        cat("  Skipped: Less than 12 columns\n")
        next
      }
      
      data <- data[!is.na(data[[12]]), ]
      data <- data[data[[12]] > 0, ]
      positions <- as.numeric(str_extract(data[[1]], "(?<=chrM:)\\d+"))
      
      temp_df <- data.frame(Position = positions, Value = data[[12]])
      scatter_data <- rbind(scatter_data, temp_df)
    }, error = function(e) {
      cat("  Skipped:", conditionMessage(e), "\n")
    })
  }
  
  scatter_data <- scatter_data[!is.na(scatter_data$Position), ]
  scatter_data$ColorGroup <- with(scatter_data, ifelse(
    Position >= 57 & Position <= 372, "HVRII",
    ifelse(Position >= 438 & Position <= 574, "HVRIII",
           ifelse(Position >= 16024 & Position <= 16383, "HVRI", "Other"))))
  scatter_data$ColorGroup <- factor(scatter_data$ColorGroup,
                                    levels = c("HVRII", "HVRIII", "HVRI", "Other"))
  
  color_map <- c(
    "HVRII" = "#FF4500",
    "HVRIII" = "#6A5ACD",
    "HVRI" = "#FFD700",
    "Other" = "#20B2AA"
  )
  
  # Count the number of samples
  sample_count <- length(file_paths)
  
  plot <- ggplot(scatter_data, aes(x = Position, y = Value)) +
    geom_point(aes(color = ColorGroup), alpha = 0.7) +
    scale_color_manual(values = color_map) +
    ylim(0, 100) +
    labs(
      title = paste0("Scatter Plot of SNPs - ", category_name, " (n=", sample_count, ")"),
      x = "Position",
      y = "Heteroplasmy %",
      color = "Genomic Region",
caption = paste("Generated on", current_date)
    ) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = NA),
      axis.line = element_line(colour = "black"),
      plot.title = element_text(hjust = 0.5)
    )
  
  # Save the plot
  output_file <- file.path(output_path, paste0("ColoredScatterPlot_", category_name, ".png"))
  ggsave(filename = output_file, plot = plot, width = 10, height = 6, dpi = 300)
}

# Split and process each category
euploid_files <- samplesheet$FullPath[samplesheet$Category == "Euploid"]
aneuploid_files <- samplesheet$FullPath[samplesheet$Category == "Aneuploid"]

process_category("Euploid", euploid_files)
process_category("Aneuploid", aneuploid_files)
