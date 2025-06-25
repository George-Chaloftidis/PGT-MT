rm(list = ls(all = TRUE))

# Load libraries
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(gridBase))
suppressPackageStartupMessages(library(rtracklayer))

#Set Date
current_date <- format(Sys.Date(), "%Y%m%d")

# Load depth data
depth <- data.table::fread("your/path/Data/depthPerPos.csv",
                           header = TRUE, stringsAsFactors = FALSE) %>% as.data.frame()

# Load sample sheet
setwd("your/path/Data")
sampleSheet <- data.table::fread("samplesheet.csv", stringsAsFactors = FALSE)

# Get sample IDs by condition
sub10x <- sampleSheet %>% filter(condition == "sub10x") %>% pull(sampleID)
seq10x <- sampleSheet %>% filter(condition == "seq10x") %>% pull(sampleID)
validation <- sampleSheet %>% filter(condition == "validation") %>% pull(sampleID)

# Calculate mean depth per group
depth$mean_seq10x <- rowMeans(depth %>% select(all_of(seq10x)))
depth$mean_sub10x <- rowMeans(depth %>% select(all_of(sub10x)))
depth$mean_validation <- rowMeans(depth %>% select(all_of(validation)))

# Prepare depth matrix
depth2 <- depth %>% select(mean_validation, mean_sub10x, mean_seq10x, POS) %>% arrange(POS)
depth2$section <- factor(1)  # Ensure a valid factor for split

# Set color schemes for coverage heatmaps
viridis_colors <- as.vector(viridis(n = 256))
fill <- colorRamp2(breaks = seq(from = 0, to = 10000, by = 10000 / 255), colors = rev(viridis_colors))
fill2 <- colorRamp2(breaks = seq(from = 0, to = 10000, by = 10000 / 255), colors = rev(viridis_colors))

# Import GFF and filter to genes, rRNA, and tRNA
gff <- import("Homo_sapiens.GRCh38.113.chromosome.MT.gff3")
gff <- gff[!is.na(gff$type), ]
gff <- gff[gff$type %in% c("gene", "rRNA", "tRNA"), ]

# Feature color map
feature_colors <- c(
  gene = "#DC267F",
  rRNA = "#FFB000",
  tRNA = "#FE6100"   
)

# Assign color to each feature
gff$color <- feature_colors[as.character(gff$type)]
gff$color[is.na(gff$color)] <- "#FCBF49"

# Combine annotations
annotations <- data.frame(
  start = start(gff),
  end = end(gff),
  type = gff$type,
  name = mcols(gff)$Name,
  pos = (start(gff) + end(gff)) / 2,
  color = gff$color
)
annotations$section <- factor(1)
annotations$name <- gsub("^MT-", "", annotations$name)

# Begin circos plot
jpeg("your/path/Output/coverageCircos.jpg",
     width = 3, height = 3, units = "in", res = 600)

circos.clear()
circos.par(start.degree = 90)

# Set the track gap
set_track_gap(gap = 0)

# Draw the heatmap tracks
circos.heatmap(mat = depth2[, "mean_sub10x", drop = FALSE], col = fill,
               split = depth2$section, track.height = 0.09, cluster = FALSE)
circos.heatmap(mat = depth2[, "mean_seq10x", drop = FALSE], col = fill2,
               split = depth2$section, track.height = 0.09, cluster = FALSE)

# Draw the GFF annotation track
circos.track(
  ylim = c(0, 1), track.height = 0.08, bg.border = NA,
  panel.fun = function(x, y) {
    last_end <- 0  
    
    for (i in seq_len(nrow(annotations))) {
      # Add gray background for each feature (even if it is a gap)
      if (annotations$start[i] > last_end) {
        # Draw a gray rectangle for the gap between features
        circos.rect(
          xleft = last_end,
          xright = annotations$start[i],
          ybottom = 0,
          ytop = 1,
          col = "gray90",
          border = NA
        )
      }
      
      # Draw the feature itself (gene, rRNA, tRNA)
      if (!is.na(annotations$color[i])) {
        circos.rect(
          xleft = annotations$start[i],
          xright = annotations$end[i],
          ybottom = 0,
          ytop = 1,
          col = annotations$color[i],
          border = NA
        )
        
        # Draw gene names for genes and rRNA
        if (annotations$type[i] %in% c("gene", "rRNA")) {
          circos.text(
            x = annotations$pos[i],
            y = 0.5,
            labels = annotations$name[i],
            facing = "bending",
            niceFacing = TRUE,
            cex = 0.3
          )
        }
      }
      
      # Update the last_end to track the last position drawn
      last_end <- annotations$end[i]
    }
    
    # Fill any remaining space after the last feature with gray
    if (last_end < max(depth2$POS)) {
      circos.rect(
        xleft = last_end,
        xright = max(depth2$POS),
        ybottom = 0,
        ytop = 1,
        col = "gray90",
        border = NA

      
      )
    }
  }
)

# Clear and finalize the plot
circos.clear()
dev.off()

# Generate legend image
jpeg("your/path/Output/coverageCircos_legend.jpg",
     width = 3.5, height = 4, units = "in", res = 600)

plot.new()

# Create manual legends for GFF feature types
feature_colors_ext <- c(feature_colors, "D-loop" = "gray")

# Create updated legend
lgd_gff <- Legend(
  at = names(feature_colors_ext),
  title = "Features",
  legend_gp = gpar(fill = feature_colors_ext),
  labels_gp = gpar(fontsize = 8)
)

lgd1 <- Legend(at = c(0, 2500, 5000, 7500, 10000), title = "sub10x", col_fun = fill, title_position = "topcenter")
lgd2 <- Legend(at = c(0, 2500, 5000, 7500, 10000), title = "seq10x", col_fun = fill2, title_position = "topcenter")

draw(packLegend(lgd1, lgd2, lgd_gff))
dev.off()
