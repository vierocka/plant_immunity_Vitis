# Load necessary packages
library(dplyr)
library(readr)
library(purrr)

# Set the directory containing your CSV files
csv_dir <- "Splicing_junctions/"  # update with your folder

# List all CSV files
SJ_files <- list.files(path = csv_dir, pattern = "\\.csv$", full.names = TRUE)

# Define a function to read and aggregate one BED file.
# Assumes:
#   - Column 2 contains gene names.
#   - Column 1 contains the score (value to sum).
read_and_aggregate <- function(file) {
  # Read the file; adjust parameters if your files have headers or different delimiters.
  df <- read.table(file, header = FALSE, stringsAsFactors = FALSE)
  # Rename columns so that V4 is gene and V5 is the score.
  colnames(df)[1:2] <- c("count", "gene")
  # Aggregate: sum score for each gene (in case a gene has multiple records)
  agg <- df %>%
    group_by(gene) %>%
    summarise(score = sum(count, na.rm = TRUE))
  # Rename the aggregated score column with the sample name (or file basename)
  sample_name <- tools::file_path_sans_ext(basename(file))
  agg <- agg %>% rename(!!sample_name := score)
  return(agg)
}

# Process each CSV file and store the results in a list
bed_list <- lapply(SJ_files, read_and_aggregate)

# Combine the aggregated data from all files by performing a full join by "gene"
combined_table <- reduce(bed_list, full_join, by = "gene")

# Replace NA values with 0 (i.e. genes missing in some BED files)
combined_table[is.na(combined_table)] <- 0

# View the final combined table
print(combined_table)
# write.table(combined_table, "data_files/SJ_allSamples.counts", sep="\t")
