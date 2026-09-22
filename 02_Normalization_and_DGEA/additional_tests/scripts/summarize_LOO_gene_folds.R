###############################################################################
# Summarize persisted gene-level leave-one-out DESeq2 results
#
# DGEA_leave_one_out_sensitivity.R writes one compressed RDS file for every
# contrast x dropped-sample fold (9 x 36 = 324). This script IMPORTS those
# files and creates auditable gene-level retention summaries without refitting.
###############################################################################

input_dir <- "02_Normalization_and_DGEA/additional_tests/results/LOO_sensitivity/gene_level_folds"
output_dir <- "02_Normalization_and_DGEA/additional_tests/results/LOO_sensitivity"
fold_files <- sort(list.files(input_dir, pattern = "\\.rds$", full.names = TRUE))
if (length(fold_files) != 324L) {
  stop("Expected 324 gene-level fold files; found ", length(fold_files), ".")
}

# Determine contrast identity from file content, not from filenames.
index_rows <- lapply(fold_files, function(path) {
  x <- readRDS(path)
  required <- c(
    "gene", "genotype", "timing", "dropped_sample", "dropped_role",
    "full_LFC", "LOO_LFC", "LFC_change", "sign_flip", "full_call",
    "LOO_call", "call_retained", "call_lost", "call_gained"
  )
  missing <- setdiff(required, names(x))
  if (length(missing)) stop(basename(path), " is missing: ", paste(missing, collapse = ", "))
  if (length(unique(x$genotype)) != 1L || length(unique(x$timing)) != 1L ||
      length(unique(x$dropped_sample)) != 1L) {
    stop(basename(path), " does not contain one contrast/fold.")
  }
  data.frame(
    path = path,
    genotype = unique(x$genotype),
    timing = unique(x$timing),
    dropped_sample = unique(x$dropped_sample),
    stringsAsFactors = FALSE
  )
})
index <- do.call(rbind, index_rows)
if (anyDuplicated(index[c("genotype", "timing", "dropped_sample")])) {
  stop("Duplicate persisted contrast/fold files detected.")
}

contrast_keys <- unique(index[c("genotype", "timing")])
contrast_keys <- contrast_keys[order(
  match(contrast_keys$genotype, c("Rpv12", "Rpv12+1", "Rpv12+1+3")),
  contrast_keys$timing
), ]

summary_rows <- vector("list", nrow(contrast_keys))
focal_rows <- list()
focal_genes <- c(
  "Vitvi17g00964", # SOBIR1
  "Vitvi09g04536", # GSO2
  "Vitvi09g04548", # LRR RLP
  "Vitvi01g04416", # LRR RLP
  "Vitvi09g01951", # RGI1
  "Vitvi09g04565", # RLP
  "Vitvi17g04216", # EDS1
  "Vitvi05g00637"  # ACD6
)

for (k in seq_len(nrow(contrast_keys))) {
  key <- contrast_keys[k, ]
  paths <- index$path[index$genotype == key$genotype & index$timing == key$timing]
  if (length(paths) != 36L) stop("A contrast does not have all 36 folds.")
  folds <- lapply(paths, readRDS)
  reference_genes <- folds[[1]]$gene
  if (!all(vapply(folds, function(x) identical(x$gene, reference_genes), logical(1)))) {
    stop("Gene order/universe differs among folds for ", key$genotype, " at ", key$timing)
  }

  full_call <- folds[[1]]$full_call
  full_lfc <- folds[[1]]$full_LFC
  if (!all(vapply(folds, function(x) identical(x$full_call, full_call), logical(1))) ||
      !all(vapply(folds, function(x) isTRUE(all.equal(x$full_LFC, full_lfc)), logical(1)))) {
    stop("Full-data baseline differs among folds for one contrast.")
  }

  retained <- rowSums(do.call(cbind, lapply(folds, `[[`, "call_retained")), na.rm = TRUE)
  lost <- rowSums(do.call(cbind, lapply(folds, `[[`, "call_lost")), na.rm = TRUE)
  gained <- rowSums(do.call(cbind, lapply(folds, `[[`, "call_gained")), na.rm = TRUE)
  sign_flips <- rowSums(do.call(cbind, lapply(folds, `[[`, "sign_flip")), na.rm = TRUE)
  delta_matrix <- do.call(cbind, lapply(folds, function(x) x$LFC_change))

  summary_rows[[k]] <- data.frame(
    gene = reference_genes,
    genotype = key$genotype,
    timing = key$timing,
    full_LFC = full_lfc,
    full_call = full_call,
    n_folds = length(folds),
    retained_call_folds = retained,
    call_retention_fraction = ifelse(full_call, retained / length(folds), NA_real_),
    lost_call_folds = lost,
    gained_call_folds = gained,
    sign_flip_folds = sign_flips,
    mean_LFC_change = rowMeans(delta_matrix, na.rm = TRUE),
    mean_abs_LFC_change = rowMeans(abs(delta_matrix), na.rm = TRUE),
    max_abs_LFC_change = apply(abs(delta_matrix), 1, max, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  for (i in seq_along(folds)) {
    focal <- folds[[i]][folds[[i]]$gene %in% focal_genes, , drop = FALSE]
    focal_rows[[length(focal_rows) + 1L]] <- focal
  }
}

gene_summary <- do.call(rbind, summary_rows)
summary_connection <- gzfile(
  file.path(output_dir, "LOO_DESeq2_gene_retention_summary.tsv.gz"), "wt"
)
write.table(gene_summary, summary_connection, sep = "\t", quote = FALSE,
            row.names = FALSE, na = "NA")
close(summary_connection)

focal <- do.call(rbind, focal_rows)
write.csv(focal,
          file.path(output_dir, "LOO_SOBIR1_network_full_vs_LOO_stability.csv"),
          row.names = FALSE)

called_summary <- gene_summary[gene_summary$full_call, ]
retention_distribution <- do.call(rbind, lapply(
  split(called_summary, interaction(called_summary$genotype,
                                    called_summary$timing, drop = TRUE)),
  function(x) data.frame(
    genotype = x$genotype[1], timing = x$timing[1],
    min = min(x$call_retention_fraction),
    q25 = unname(quantile(x$call_retention_fraction, 0.25)),
    median = median(x$call_retention_fraction),
    mean = mean(x$call_retention_fraction),
    q75 = unname(quantile(x$call_retention_fraction, 0.75)),
    max = max(x$call_retention_fraction),
    stringsAsFactors = FALSE
  )
))
rownames(retention_distribution) <- NULL
write.csv(retention_distribution,
          file.path(output_dir, "LOO_DESeq2_call_retention_distribution.csv"),
          row.names = FALSE)
writeLines(capture.output(sessionInfo()),
           file.path(output_dir, "summarize_gene_folds_sessionInfo.txt"))
message("Completed LOO gene-level summaries in: ", output_dir)
