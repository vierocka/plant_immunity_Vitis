###############################################################################
# Task 6: bootstrap stability analysis (200x resample the 36 samples WITH
# replacement, rerun the extracted gcna_module_builder.R pipeline), for the
# same 4 base conditions as 00_reproduce_baseline.R / 05_noise_perturbation.R.
#
# Resampling is of SAMPLES (columns), not genes: draw 36 sample indices with
# replacement, build a resampled expression matrix using those columns, and
# apply the SAME resampled indices to each anchor's trait vector (trait
# values are tied to sample identity, so they must be permuted identically
# to the expression columns, not independently).
#
# CHECKPOINTED / RESUMABLE, same design as 05_noise_perturbation.R (~6h
# estimated total - shorter than the noise loop since there is only one
# "noise_level" (the resampling itself), not 4).
###############################################################################

source("network_robustness/gcna_module_builder.R")

N_REPS <- 200
output_dir <- "network_robustness/06_bootstrap"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_expr <- function(path) {
  df <- read.table(path, header = TRUE, sep = "\t")
  mat <- as.matrix(df[, 2:37]); rownames(mat) <- df[, 1]
  mat
}
expr <- list(
  unprotected = read_expr("../data_files/Rlogs.csv"),
  protected   = read_expr("../data_files/Rlogs_ComBat_protected.csv")
)
panels <- load_de_panels(repo_root = "..")

jaccard <- function(a, b) if (length(union(a, b)) == 0) NA_real_ else length(intersect(a, b)) / length(union(a, b))

args <- commandArgs(trailingOnly = TRUE)
combat_targets <- if (length(args) >= 1) args[1] else names(expr)
panel_targets  <- if (length(args) >= 2) args[2] else names(panels)

for (combat_label in combat_targets) {
  for (panel_label in panel_targets) {
    condition_label <- paste0(combat_label, "_", panel_label)
    baseline_members <- readRDS(file.path("network_robustness/00_baseline",
                                           paste0("module_members_", condition_label, ".rds")))
    baseline_anchors <- names(baseline_members)
    base_mat <- expr[[combat_label]]
    anchor_genes <- panels[[panel_label]]$anchor_genes
    trait_fn <- panels[[panel_label]]$trait_fn

    checkpoint_file <- file.path(output_dir, sprintf("bootstrap_%s.csv", condition_label))
    done_reps <- integer(0)
    if (file.exists(checkpoint_file)) {
      existing <- read.csv(checkpoint_file, stringsAsFactors = FALSE)
      done_reps <- existing$rep
      message(condition_label, ": resuming, ", length(done_reps), "/", N_REPS, " reps already done.")
    }
    remaining <- setdiff(seq_len(N_REPS), done_reps)
    if (!length(remaining)) {
      message(condition_label, ": already complete, skipping.")
      next
    }

    for (r in remaining) {
      set.seed(20000L * r + nchar(condition_label))
      boot_idx <- sample.int(ncol(base_mat), ncol(base_mat), replace = TRUE)
      boot_mat <- base_mat[, boot_idx, drop = FALSE]
      boot_trait_fn <- function(gene) {
        t <- trait_fn(gene)
        if (is.null(t)) return(NULL)
        t[boot_idx]
      }

      res <- build_gcna_modules(boot_mat, anchor_genes, boot_trait_fn,
                                 r_cutoff = 0.817, min_partners = 3, verbose = FALSE)
      sig_anchors <- names(res$module_members)
      sizes <- vapply(res$module_members, length, integer(1))

      anchor_jac <- jaccard(baseline_anchors, sig_anchors)
      shared <- intersect(baseline_anchors, sig_anchors)
      partner_jac <- if (length(shared)) {
        mean(vapply(shared, function(g) jaccard(baseline_members[[g]], res$module_members[[g]]), numeric(1)))
      } else NA_real_

      row <- data.frame(rep = r, n_significant = length(sig_anchors),
                         median_module_size = if (length(sizes)) median(sizes) else NA_real_,
                         anchor_identity_jaccard = anchor_jac,
                         mean_partner_set_jaccard = partner_jac, stringsAsFactors = FALSE)
      write.table(row, checkpoint_file, sep = ",", row.names = FALSE,
                  col.names = !file.exists(checkpoint_file), append = file.exists(checkpoint_file))
    }
    message(condition_label, ": complete (", N_REPS, " reps).")
  }
}

############################ SUMMARY (only if all checkpoints complete) ########
all_files <- list.files(output_dir, pattern = "^bootstrap_.*\\.csv$", full.names = TRUE)
summary_rows <- lapply(all_files, function(f) {
  d <- read.csv(f, stringsAsFactors = FALSE)
  if (nrow(d) < N_REPS) return(NULL)
  condition <- sub("^bootstrap_(.*)\\.csv$", "\\1", basename(f))
  data.frame(condition = condition, n_reps = nrow(d),
             mean_n_significant = mean(d$n_significant), sd_n_significant = sd(d$n_significant),
             ci95_low = quantile(d$n_significant, 0.025), ci95_high = quantile(d$n_significant, 0.975),
             mean_anchor_jaccard = mean(d$anchor_identity_jaccard, na.rm = TRUE),
             mean_partner_jaccard = mean(d$mean_partner_set_jaccard, na.rm = TRUE),
             stringsAsFactors = FALSE)
})
summary_rows <- summary_rows[!vapply(summary_rows, is.null, logical(1))]
if (length(summary_rows)) {
  summary_df <- do.call(rbind, summary_rows)
  write.csv(summary_df, file.path(output_dir, "bootstrap_summary.csv"), row.names = FALSE)
  message("\n=== Bootstrap summary (conditions with all ", N_REPS, " reps complete) ===")
  print(summary_df, digits = 3)
}
message("Done (or partially done - rerun to resume any incomplete checkpoints).")
