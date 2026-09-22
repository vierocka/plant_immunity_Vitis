###############################################################################
# Task 5: noise-perturbation robustness loop (0.5/1/2/5% per-gene Gaussian
# jitter, 200 reps/level/condition), wrapping the extracted, unit-tested
# gcna_module_builder.R around the SAME 4 base conditions as
# 00_reproduce_baseline.R. Noise model: N(0, level * SD_gene) added
# independently per gene. Measures Jaccard stability of module membership vs the r=0.817,
# no-noise baseline already verified to reproduce 155/147 modules.
#
# CHECKPOINTED / RESUMABLE: each (condition, noise_level) writes one
# checkpoint CSV, appended after every repetition. On restart, already-
# completed repetition indices are skipped - required given the ~24h
# estimated total runtime (150,000+ times longer than a typical interactive
# session; must tolerate being interrupted and resumed).
#
# Accepts optional CLI args to allow partial/targeted runs, e.g.:
#   Rscript 05_noise_perturbation.R unprotected 3553
# reruns just one condition (all 4 noise levels) - default (no args) runs
# all 4 conditions.
###############################################################################

source("network_robustness/gcna_module_builder.R")

N_REPS <- 200
NOISE_LEVELS <- c(0.005, 0.01, 0.02, 0.05)
output_dir <- "network_robustness/05_noise_perturbation"
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
    gene_sd <- apply(base_mat, 1, sd)
    anchor_genes <- panels[[panel_label]]$anchor_genes
    trait_fn <- panels[[panel_label]]$trait_fn

    for (noise in NOISE_LEVELS) {
      checkpoint_file <- file.path(output_dir, sprintf("noise_%s_level%.3f.csv", condition_label, noise))
      done_reps <- integer(0)
      if (file.exists(checkpoint_file)) {
        existing <- read.csv(checkpoint_file, stringsAsFactors = FALSE)
        done_reps <- existing$rep
        message(condition_label, " noise=", noise, ": resuming, ", length(done_reps), "/", N_REPS, " reps already done.")
      }
      remaining <- setdiff(seq_len(N_REPS), done_reps)
      if (!length(remaining)) {
        message(condition_label, " noise=", noise, ": already complete, skipping.")
        next
      }

      for (r in remaining) {
        set.seed(10000L * r + round(noise * 1000) + nchar(condition_label))
        jitter <- matrix(rnorm(length(base_mat), mean = 0, sd = rep(noise * gene_sd, ncol(base_mat))),
                          nrow = nrow(base_mat), dimnames = dimnames(base_mat))
        noisy_mat <- base_mat + jitter

        res <- build_gcna_modules(noisy_mat, anchor_genes, trait_fn,
                                   r_cutoff = 0.817, min_partners = 3, verbose = FALSE)
        sig_anchors <- names(res$module_members)

        anchor_jac <- jaccard(baseline_anchors, sig_anchors)
        shared <- intersect(baseline_anchors, sig_anchors)
        partner_jac <- if (length(shared)) {
          mean(vapply(shared, function(g) jaccard(baseline_members[[g]], res$module_members[[g]]), numeric(1)))
        } else NA_real_

        row <- data.frame(rep = r, n_significant = length(sig_anchors),
                           anchor_identity_jaccard = anchor_jac,
                           mean_partner_set_jaccard = partner_jac, stringsAsFactors = FALSE)
        write.table(row, checkpoint_file, sep = ",", row.names = FALSE,
                    col.names = !file.exists(checkpoint_file), append = file.exists(checkpoint_file))
      }
      message(condition_label, " noise=", noise, ": complete (", N_REPS, " reps).")
    }
  }
}

############################ SUMMARY (only if all checkpoints complete) ########
all_files <- list.files(output_dir, pattern = "^noise_.*\\.csv$", full.names = TRUE)
summary_rows <- lapply(all_files, function(f) {
  d <- read.csv(f, stringsAsFactors = FALSE)
  if (nrow(d) < N_REPS) return(NULL)
  info <- regmatches(basename(f), regexec("^noise_(.*)_level([0-9.]+)\\.csv$", basename(f)))[[1]]
  data.frame(condition = info[2], noise_level = as.numeric(info[3]),
             n_reps = nrow(d),
             mean_n_significant = mean(d$n_significant), sd_n_significant = sd(d$n_significant),
             mean_anchor_jaccard = mean(d$anchor_identity_jaccard, na.rm = TRUE),
             mean_partner_jaccard = mean(d$mean_partner_set_jaccard, na.rm = TRUE),
             stringsAsFactors = FALSE)
})
summary_rows <- summary_rows[!vapply(summary_rows, is.null, logical(1))]
if (length(summary_rows)) {
  summary_df <- do.call(rbind, summary_rows)
  write.csv(summary_df, file.path(output_dir, "noise_perturbation_summary.csv"), row.names = FALSE)
  message("\n=== Noise-perturbation summary (conditions with all ", N_REPS, " reps complete) ===")
  print(summary_df, digits = 3)
}
message("Done (or partially done - rerun to resume any incomplete checkpoints).")
