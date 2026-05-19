# Run this script ONCE locally to convert data files into a DuckDB database.
# After running, commit vitis_data.duckdb to the repo / include it in shinyapps.io deployment.
# You do NOT need to re-run this unless the source files change.
# Then, continue with the shiny app.

library(duckdb)
library(readxl)

db_path <- file.path(dirname(rstudioapi::getSourceEditorContext()$path), "vitis_data.duckdb")
con <- dbConnect(duckdb(), db_path)

message("Loading metadata xlsx...")
metadata <- readxl::read_xlsx("26169genes_with_AthalHomologs_allIDs_exprPatterns.xlsx")
dbWriteTable(con, "metadata", metadata, overwrite = TRUE)
message(sprintf("  -> %d rows, %d columns", nrow(metadata), ncol(metadata)))

message("Loading Rlogs.csv...")
rlogs <- data.table::fread("Rlogs.csv", sep = "\t")
dbWriteTable(con, "rlogs", as.data.frame(rlogs), overwrite = TRUE)
message(sprintf("  -> %d rows, %d columns", nrow(rlogs), ncol(rlogs)))

message("Loading RawCounts.csv...")
raw_counts <- data.table::fread("RawCounts.csv", sep = "\t")
dbWriteTable(con, "raw_counts", as.data.frame(raw_counts), overwrite = TRUE)
message(sprintf("  -> %d rows, %d columns", nrow(raw_counts), ncol(raw_counts)))

message("Loading TranscrpModule_Spearman_Pval_Padj_FDR.csv...")
module_info <- read.csv("TranscrpModule_Spearman_Pval_Padj_FDR.csv")
dbWriteTable(con, "module_info", module_info, overwrite = TRUE)
message(sprintf("  -> %d rows, %d columns", nrow(module_info), ncol(module_info)))

# Create indexes for fast single-gene lookup (used on every row click)
message("Creating indexes...")
dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_rlogs_gene      ON rlogs(geneID)")
dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_module_gene     ON module_info(geneID)")
dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_metadata_gene   ON metadata(PN40024v4_ENSMBL)")

dbDisconnect(con, shutdown = TRUE)
message("Done. vitis_data.duckdb created at: ", db_path)
