###############################################################################
# Network checks on the CURRENT network without the module-trait test (2026-10-22).
# Modules = every DE anchor (9,459 canonical DESeq2 genes) with >= 3 partners at r >= cutoff
# (rlog, condition-protected ComBat, all 26,169 genes); no trait-association filter.
# Replaces, for the manuscript: FDR-"significant" module counts and the trait-shuffle permutation
# (circular: the trait comes from the same data), and the 3,553-gene-panel numbers formerly
# quoted for STRING enrichment (script 12) and PC1 variance (script 11).
# Outputs (15_candidate_modules/): module sizes per cutoff (Figure S7B), member lists (Figure S7C,
# read by string-db/scripts/SF6_network_robustness_rebuild.R), PC1 variance and size-matched
# STRING score test for the r = 0.8077 modules and for protected WGCNA modules.
# Run from 06_PCNWA/:  Rscript network_robustness/15_candidate_modules_no_trait_filter.R
###############################################################################
Sys.setenv(OPENBLAS_NUM_THREADS = "1")
source("network_robustness/gcna_module_builder.R")
out <- "network_robustness/15_candidate_modules"; dir.create(out, showWarnings = FALSE)
set.seed(20261022)

expr_df <- read.delim("../data_files/Rlogs_ComBat_protected.csv", check.names = FALSE)
expr <- as.matrix(expr_df[, -1]); rownames(expr) <- as.character(expr_df[[1]])
stopifnot(!anyNA(expr), ncol(expr) == 36, nrow(expr) == 26169)
anchors <- as.character(read.delim("../05_transcriptional_dynamics/canonical_DESeq2/tables/canonical_primary_DEG_direction_matrix.tsv")$gene)
stopifnot(length(anchors) == 9459)
r_main <- read.csv("combat_protected/results/GCNA_rebuild_9459panel_FDR_datadriven_rcutoff/r_cutoff_derivation.csv")$r_cutoff_top0.5pct
cutoffs <- c("0.757" = 0.757, "0.8077" = r_main, "0.8839" = 0.8839)   # top 1%, 0.5%, 0.1% of all pairwise r

# ---- modules at the three cutoffs (stage A once, stage B per cutoff) ----
ac <- compute_anchor_correlation(expr, anchors)
modules <- lapply(cutoffs, function(rc) lapply(extract_partners_pc1(ac, expr, r_cutoff = rc, min_partners = 3), `[[`, "partners"))
rm(ac); invisible(gc())
for (k in names(modules)) saveRDS(modules[[k]], file.path(out, paste0("module_members_r", k, ".rds")), compress = "xz")
sizes <- do.call(rbind, lapply(names(modules), function(k)
  data.frame(r = k, anchor = names(modules[[k]]), n_partners = lengths(modules[[k]]))))
write.csv(sizes, file.path(out, "module_sizes_by_cutoff.csv"), row.names = FALSE)
print(aggregate(n_partners ~ r, sizes, function(x) c(n = length(x), median = median(x), max = max(x))))

# ---- PC1 variance explained (same definition as script 11), r = 0.8077 modules vs protected WGCNA ----
pve <- function(genes) { pc <- prcomp(t(expr[genes, , drop = FALSE]), center = TRUE, scale. = FALSE); pc$sdev[1]^2 / sum(pc$sdev^2) }
m8077 <- Map(c, names(modules[["0.8077"]]), modules[["0.8077"]])          # module gene set = anchor + partners
wcol <- readRDS("WGCNA/results/WGCNA_classical_comparison_protected/WGCNA_module_colors_power8.rds")
wgcna <- split(names(wcol), wcol)[setdiff(unique(wcol), "grey")]
pv <- rbind(data.frame(method = "correlation_modules_r0.8077", module = names(m8077), size = lengths(m8077),
                       pve_pc1 = unlist(mclapply(m8077, pve, mc.cores = default_mc_cores()))),
            data.frame(method = "WGCNA_power8_protected", module = names(wgcna), size = lengths(wgcna),
                       pve_pc1 = vapply(wgcna, pve, numeric(1))))
write.csv(pv, file.path(out, "module_pc1_variance_explained.csv"), row.names = FALSE)
print(aggregate(pve_pc1 ~ method, pv, function(x) c(mean = mean(x), median = median(x), n = length(x))))

# ---- size-matched STRING score test (setup = script 12, sections 1-2; raw-data paths updated) ----
setup <- readLines("network_robustness/12_string_score_permutation.R")
setup <- setup[seq_len(grep("^#+ 3\\. LOAD REAL MODULES", setup) - 1)]
setup <- gsub('"string-db/29760', '"string-db/raw_data/29760', setup)
eval(parse(text = setup))
all_sets <- list(correlation_modules_r0.8077 = m8077, WGCNA_power8_protected = wgcna)
sz <- unlist(lapply(all_sets, lengths))
bin_edges <- unique(round(exp(seq(log(max(2, min(sz))), log(max(sz)), length.out = 30))))
nearest <- function(n) bin_edges[which.min(abs(bin_edges - n))]
null_by_bin <- mclapply(setNames(bin_edges, bin_edges), function(n) lapply(seq_len(K_DRAWS), function(k) {
  v <- logical(N); v[sample.int(N, n)] <- TRUE; module_stats_from_membership(v) }), mc.cores = default_mc_cores(), mc.set.seed = TRUE)
res <- do.call(rbind, lapply(names(all_sets), function(grp) do.call(rbind, mclapply(names(all_sets[[grp]]), function(m) {
  g <- all_sets[[grp]][[m]]; obs <- module_stats_from_genes(g); nul <- null_by_bin[[as.character(nearest(length(g)))]]
  row <- data.frame(method = grp, module = m, size = length(g), n_mapped = obs$n_mapped)
  for (nm in names(score_tables)) {
    o <- obs[[nm]]["mean_score"]; nm_ <- vapply(nul, function(d) d[[nm]]["mean_score"], numeric(1)); nm_ <- nm_[!is.na(nm_)]
    row[[paste0(nm, "_p")]] <- if (is.na(o) || !length(nm_)) NA_real_ else (1 + sum(nm_ >= o)) / (1 + length(nm_))
  }
  row }, mc.cores = default_mc_cores()))))
write.csv(res, file.path(out, "module_string_score_permutation_results.csv"), row.names = FALSE)
summ <- do.call(rbind, lapply(names(score_tables), function(nm) do.call(rbind, lapply(names(all_sets), function(grp) {
  p <- res[res$method == grp, paste0(nm, "_p")]; p <- p[!is.na(p)]
  data.frame(score_type = nm, method = grp, n_modules_tested = length(p), frac_p_lt_0.05 = mean(p < 0.05)) }))))
write.csv(summ, file.path(out, "module_string_score_permutation_summary.csv"), row.names = FALSE)
print(summ, digits = 3)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
