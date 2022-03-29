# A script to read in slurm output files of running samtools flagstat on bwa mem mappings of walleye RAPTURE data to the yellow perch genome. It returns a simple mean and standard deviation for mapping percentages.
library(tidyverse)

# List the slurm files
slurm_files <- list.files(path = "slurm_out/", pattern = "*.out", full.names = TRUE)

# Create an empty data frame
mapping_data <- tibble(percent = rep(NA, length(slurm_files)))

# Loop through the slurm files, splitting up row 6 (the line with the mapping % information) and taking only the number
for (file in 1:length(slurm_files)){
  
  print(slurm_files[file]) # Track where we are in the loop
  
  data <- read_csv(slurm_files[file], col_names = FALSE, show_col_types = FALSE) %>% 
    slice(6) %>% # Select row 6
    str_split(., "\\),\\(|\\)|\\(", simplify = TRUE) # Split things by parentheses
  
  data <- tibble(data[2]) %>% # Take the data and turn it back into a tibble
    str_split(., "\\%", simplify = TRUE) # Split things by the percent sign
  
  mapping_data[file, 1] <- as.numeric(data[,1]) # Take the number, turn it into a numeric value and put it into the output data frame
  
  rm(data)
}
rm(file, slurm_files)

# Calculate simple summary statistics of mapping percents
mean(mapping_data$percent)
sd(mapping_data$percent)
