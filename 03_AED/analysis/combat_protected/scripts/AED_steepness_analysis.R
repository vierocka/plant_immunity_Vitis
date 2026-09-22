###############################################################################
# Steepness of the AED trajectories -- rate of change of divergence per unit
# of ordinal stage, for both the raw AggrDiv and the null-relative z-score.
#
# FRAMING: every real 3-vs-3 biological comparison in
# this project's AED tests lands in the right tail of its own permutation
# null (p_emp always saturates at the design's floor, ~2/(n_null+1)) -- this
# is not an artifact to explain away, it is the expected signature of real,
# consistent triplicate structure: a genuine biological grouping should
# reliably out-diverge random same-pool resamplings, and it does, every
# time. What DOES vary test-to-test, and is not saturated, is HOW FAR into
# the tail the observed value sits -- the z-score / fold-over-null-median
# computed in AED_zscore_comparability_check.R. A small value there (z close
# to the bulk of the null) is the informative "this looks more like noise /
# unstructured global stress than a distinct, right-tail-shifting biological
# program" signal the user described. Steepness = how fast that
# distance-into-the-tail (or, secondarily, raw magnitude) grows from one
# measurement to the next -- a positive, accelerating slope is the direct
# quantitative form of "the shift becomes more pronounced at later times".
###############################################################################

output_dir <- "03_AED/analysis/combat_protected/tables"
d <- read.csv(file.path(output_dir, "AED_zscore_comparability_all_tests.csv"), stringsAsFactors = FALSE)

mock_order <- c(12, 24, 48, 96, 120)   # within_mock's own-hpi order -> ordinal position
inoc_order <- c(24, 48, 96)            # within_inoculated's own-hpi order -> ordinal position (baseline=12h)
stage_of <- function(dataset, family, contrast) {
  # NOTE: dataset=="OwnStudy" is NOT specific enough on its own -- it also
  # matches cross_genotype_vs_Susceptible (also OwnStudy) and would silently
  # swallow it into the wrong branch (caught 2026-08-17: every
  # cross_genotype contrast except "6hpi_vs_0hpi" literally was mapping to
  # stage=2 regardless of its real hpi, because this branch matched first --
  # visible as d_stage=0 / Inf slopes once printed). Must check the specific
  # within_genotype_temporal family, not just the dataset.
  if (family == "within_Syrah_background") as.integer(sub("^DEV([0-9]+)_vs_DEV01$", "\\1", contrast))
  else if (dataset == "OwnStudy" && grepl("^within_(Susceptible|Rpv12(_1(_3)?)?)$", family)) ifelse(contrast == "6hpi_vs_0hpi", 1, 2)
  else if (family == "within_mock") match(as.integer(sub("^([0-9]+)h_vs_0h$", "\\1", contrast)), mock_order)
  else if (family == "within_inoculated") match(as.integer(sub("^([0-9]+)h_vs_12h$", "\\1", contrast)), inoc_order)
  else if (family == "cross_condition_inoculated_vs_mock") as.integer(sub("hpi$", "", contrast))  # real hpi kept, like cross_cultivar's real T-index
  else if (family == "cross_genotype_vs_Susceptible") as.integer(sub(".*_([0-9]+)hpi$", "\\1", contrast))  # real hpi kept
  else as.integer(sub("^T([0-9]+).*$", "\\1", contrast))
}
d$stage <- mapply(stage_of, d$dataset, d$family, d$contrast)
# cross_genotype_vs_Susceptible needs one series PER GENOTYPE (Rpv12/Rpv12_1/
# Rpv12_1_3 each have their own 0/6/24hpi trajectory vs Susceptible) --
# every other family is already one series per family.
cgeno_series <- function(family, contrast) {
  if (family != "cross_genotype_vs_Susceptible") return(NA_character_)
  paste0("OwnStudy_cross_genotype_", sub("_[0-9]+hpi$", "", contrast))
}
d$series <- mapply(cgeno_series, d$family, d$contrast)
d$series <- ifelse(is.na(d$series), ifelse(d$dataset == "OwnStudy", paste0("OwnStudy_", d$family), d$family), d$series)

# Every WITHIN-cultivar series' trajectory implicitly starts at (stage=0,
# value=0) -- the baseline vs itself (AggrDiv==0 by the run_aed() self-
# check) -- added here so the first interval's steepness is computed on the
# same footing as every later one. NOT added for cross_cultivar_MV102_vs_MV32:
# that family's stage=0 is a REAL measured test (MV102 vs MV32 at T0, both
# real samples, AggrDiv=1.036), not a trivial self-comparison -- it is
# already present in `d` and must not be overwritten with a synthetic 0.
# cross_cultivar_MV102_vs_MV32 and cross_condition_inoculated_vs_mock never
# get a synthetic baseline -- both are real cross-condition/cross-cultivar
# comparisons at every stage already, not trivial self-comparisons.
cross_type_families <- c("cross_cultivar_MV102_vs_MV32", "cross_condition_inoculated_vs_mock", "cross_genotype_vs_Susceptible")
cross_type_series <- unique(d$series[d$family %in% cross_type_families])  # resolves per-genotype series names too
series_needing_baseline <- unique(d$series[!d$family %in% cross_type_families])
baseline_rows <- do.call(rbind, lapply(series_needing_baseline, function(s) {
  r <- d[d$series == s, ][1, ]
  r$stage <- 0; r$AggrDiv <- 0; r$z_score <- 0; r$contrast <- "(own baseline)"
  r
}))
d2 <- rbind(d[, c("series","dataset","family","contrast","stage","AggrDiv","z_score")],
            baseline_rows[, c("series","dataset","family","contrast","stage","AggrDiv","z_score")])

steepness <- do.call(rbind, lapply(unique(d2$series), function(s) {
  dd <- d2[d2$series == s, ]
  dd <- dd[order(dd$stage), ]
  if (nrow(dd) < 2) return(NULL)
  data.frame(
    series = s,
    from_stage = head(dd$stage, -1), to_stage = tail(dd$stage, -1),
    to_contrast = tail(dd$contrast, -1),
    d_stage = diff(dd$stage),
    d_AggrDiv = diff(dd$AggrDiv), slope_AggrDiv_per_stage = diff(dd$AggrDiv) / diff(dd$stage),
    d_zscore = diff(dd$z_score), slope_zscore_per_stage = diff(dd$z_score) / diff(dd$stage),
    stringsAsFactors = FALSE
  )
}))
write.csv(steepness, file.path(output_dir, "AED_steepness_per_interval.csv"), row.names = FALSE)
cat("=== Steepness per interval (both raw AggrDiv/stage and z-score/stage) ===\n")
print(steepness, digits = 3)

cat("\n=== Mean steepness per series, ALL intervals (z-score/stage) -- MISLEADING, see correction below ===\n")
mean_steep <- aggregate(slope_zscore_per_stage ~ series, data = steepness, mean)
mean_steep <- mean_steep[order(-mean_steep$slope_zscore_per_stage), ]
print(mean_steep, digits = 3)

# CORRECTION: the first interval of every within-cultivar series runs from
# the synthetic (stage=0, z=0) self-comparison baseline to the first real
# measurement. That jump is close to the largest possible move (0 -> ~2-3)
# for almost ANY real triplicate effect, essentially by construction (a real
# 3-vs-3 grouping beats a random pool draw immediately) -- it swamps the
# per-series mean and answers "is there a real effect at all" (yes,
# everywhere), not "does the shift get MORE pronounced at later times"
# (the user's actual question). Excluding it isolates genuine
# acceleration/deceleration AFTER the initial response.
# cross_cultivar_MV102_vs_MV32 never had a synthetic baseline added (its
# stage=0 IS a real MV102-vs-MV32 test, not a trivial self-comparison), so
# all 3 of its intervals are kept.
steepness$is_first_interval <- with(steepness, ave(seq_along(from_stage), series, FUN = function(i) i == min(i)))
steepness_post <- steepness[steepness$is_first_interval == 0 | steepness$series %in% cross_type_series, ]
mean_steep_post <- aggregate(slope_zscore_per_stage ~ series, data = steepness_post, mean)
mean_steep_post <- mean_steep_post[order(-mean_steep_post$slope_zscore_per_stage), ]
write.csv(mean_steep_post, file.path(output_dir, "AED_steepness_post_initial_jump_mean.csv"), row.names = FALSE)

cat("\n=== CORRECTED: mean steepness EXCLUDING the trivial from-baseline first jump ===\n")
cat("(this is the metric that actually answers 'does divergence keep growing at later times')\n")
print(mean_steep_post, digits = 3)
cat("\nPositive & large = divergence keeps growing after the initial response (structured, progressive).\nNear-zero = plateaus after the initial jump. Negative = partially REVERTS toward baseline after\nthe initial response -- also a real, informative pattern, not just 'no effect'.\n")

############################ PLOT: STEEPNESS PER INTERVAL, PER SERIES ##########
suppressPackageStartupMessages(library(Cairo))
CairoPDF(file.path(output_dir, "AED_steepness_plot.pdf"), width = 11, height = 7)
par(mfrow = c(1, 2), mar = c(8, 5, 5, 1))

series_order1 <- mean_steep$series
cols1 <- colorRampPalette(c("firebrick", "gray70", "steelblue"))(length(series_order1))
names(cols1) <- series_order1
barplot(mean_steep$slope_zscore_per_stage, names.arg = mean_steep$series, las = 2, cex.names = 0.6,
        col = cols1[mean_steep$series], ylab = "Mean z-score gained per ordinal stage step",
        main = "ALL intervals (MISLEADING)\ndominated by the trivial 0-to-first-measurement jump\n-- answers 'is there any real effect', not steepness",
        cex.main = 0.68)
abline(h = 0, lty = 2)

series_order2 <- mean_steep_post$series
cols2 <- colorRampPalette(c("firebrick", "gray70", "steelblue"))(length(series_order2))
names(cols2) <- series_order2
barplot(mean_steep_post$slope_zscore_per_stage, names.arg = mean_steep_post$series, las = 2, cex.names = 0.6,
        col = cols2[mean_steep_post$series], ylab = "Mean z-score gained per ordinal stage step (post-initial-jump)",
        main = "CORRECTED: excludes trivial first jump\n-- this is the real 'does divergence keep growing\nat later times' answer",
        cex.main = 0.68)
abline(h = 0, lty = 2)
dev.off()
message("\nWrote AED_steepness_per_interval.csv, AED_steepness_post_initial_jump_mean.csv, AED_steepness_plot.pdf")
