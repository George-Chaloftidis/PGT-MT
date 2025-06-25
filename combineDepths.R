rm(list=ls(all=T))
# loadPackagesx
library(dplyr)
library(data.table)

file_path = "your/path/Data"
  setwd(file_path)


# collect the data from all embryos & combine

toRun = data.table::fread("samplesheet.csv", stringsAsFactors = F) 

for(i in 1:nrow(toRun)){
  
  sample = toRun[i, "sampleID"]
  sample_path = paste0(toRun[i, "dataDir"], "/", sample, "/MT_depth")
  cat("Looking for file:", sample_path, "\n")
  
  if (dir.exists(sample_path)) {
    setwd(sample_path)
  } else {
    warning("Directory does not exist: ", sample_path)
    next  
  }

  
  
  file_path = paste0(sample, "_MT.txt")
  
  if (file.exists(file_path)) {
    data = data.table::fread(file_path)
  } else {
    warning("File does not exist or is unreadable: ", file_path)
    next
  }
  
  data = data.table::fread(paste0(sample, "_MT.txt"))
  
  colnames(data) = gsub(".bam", "", colnames(data))
    
  if(i == 1){
    mtCoverage = data
  }else{
    
    mtCoverage = full_join(mtCoverage, data, by = c("#CHROM", "POS"))
  }
  
}

# write out the data
write.csv(mtCoverage, "your/path/Output/depthPerPos.csv", row.names = F)
write.csv(mtCoverage, "your/path/Data/depthPerPos.csv", row.names = F)


