###############################################################################
# Poisson (log-link, exposure-offset) model of STRING-confirmed edge counts:
# raw precision percentages (n_overlap / n_method_edges)
# are not directly comparable between the hub-anchored method (330k-1.36M
# candidate edges) and WGCNA's all-pairs-within-module proxy (14.6-17.4M
# edges) without accounting for that ~20-40x difference in "exposure" - a
# naive comparison lets WGCNA's much larger haystack look artificially
# similar or different by chance alone (a "how much are we testing"
# concern, distinct from raw precision).
#
# MODEL: n_overlap ~ Poisson(mu), log(mu) = log(n_method_edges) + predictors
# i.e. n_method_edges enters as a fixed-coefficient-1 OFFSET, not a covariate
# - this is the standard rate-ratio / standardized-incidence-ratio design for
# comparing two count totals measured over two different "exposures" (here:
# number of candidate edges each method actually proposed). The fitted
# coefficient on `method` is then a log(rate ratio) directly comparable
# between methods regardless of network size, with a Wald p-value.
#
# THREE ANALYSES:
#   A. hub-anchored (at the published r=0.817) vs WGCNA, per STRING reference,
#      per ComBat setting - the most granular, directly interpretable test.
#   B. Same, pooled across ComBat (combat as a covariate) - one rate ratio +
#      p-value per reference.
#   C. Trend test, hub-anchored method only: does precision (rate) increase
#      significantly with r_cutoff (0.75->0.85)? Turns the earlier
#      descriptive "rose modestly and monotonically" claim into a formal
#      slope + p-value.
#
# Dispersion is checked (Pearson chi-sq / df) and quasipoisson used instead
# of plain Poisson wherever overdispersion is detected, which inflates
# standard errors appropriately rather than reporting falsely narrow p-values.
###############################################################################

overlap <- read.csv("network_robustness/string_reference/hub_method_vs_wgcna_vs_string_overlap.csv",
                     stringsAsFactors = FALSE)
overlap$r_cutoff <- suppressWarnings(as.numeric(sub("r>=", "", overlap$parameter)))

fit_poisson_offset <- function(data, formula_rhs) {
  f <- as.formula(paste("n_overlap ~", formula_rhs))
  m_pois <- glm(f, data = data, family = poisson(link = "log"),
                offset = log(n_method_edges))
  disp <- sum(residuals(m_pois, type = "pearson")^2) / df.residual(m_pois)
  if (disp > 1.5) {
    m <- glm(f, data = data, family = quasipoisson(link = "log"),
             offset = log(n_method_edges))
    family_used <- sprintf("quasipoisson (dispersion=%.2f)", disp)
  } else {
    m <- m_pois
    family_used <- sprintf("poisson (dispersion=%.2f, no overdispersion correction needed)", disp)
  }
  list(model = m, family_used = family_used)
}

references <- unique(overlap$reference)

############################ A. single fully-pooled headline model #############
# NOTE: a per-(combat x reference) cell model was tried first and is not
# statistically usable - with exactly 1 hub-anchored point and 1 WGCNA point
# per cell, a method-effect model there is fully saturated (0 residual df),
# so no dispersion or CI can be estimated. Pooling across combat AND
# reference (as covariates) gives real residual df and one headline,
# directly interpretable rate ratio + p-value for "is the hub-anchored
# method's STRING-precision different from WGCNA's, at r=0.817, adjusting
# for ComBat setting and which STRING reference is used".
message("=== A. hub-anchored (r=0.817) vs WGCNA, single pooled model (all references + both ComBat settings as covariates) ===")
sub <- overlap[(overlap$method == "hub_anchored" & overlap$r_cutoff == 0.817) | overlap$method == "WGCNA", ]
sub$method <- relevel(factor(sub$method), ref = "WGCNA")
fit <- fit_poisson_offset(sub, "method + combat + reference")
co <- summary(fit$model)$coefficients
rr <- exp(co["methodhub_anchored", "Estimate"])
ci <- exp(co["methodhub_anchored", "Estimate"] + c(-1.96, 1.96) * co["methodhub_anchored", "Std. Error"])
pval <- co["methodhub_anchored", grepl("^Pr", colnames(co))]
resultsA_df <- data.frame(
  model = "pooled (method + combat + reference)", family = fit$family_used,
  rate_ratio_hub_vs_WGCNA = rr, ci_low = ci[1], ci_high = ci[2], p_value = pval,
  stringsAsFactors = FALSE
)
print(resultsA_df, digits = 3)

############################ B. per-reference, ComBat pooled as covariate ######
message("\n=== B. hub-anchored (r=0.817) vs WGCNA, per reference, ComBat pooled ===")
resultsB <- list()
for (ref in references) {
  sub <- overlap[overlap$reference == ref &
                  ((overlap$method == "hub_anchored" & overlap$r_cutoff == 0.817) |
                   overlap$method == "WGCNA"), ]
  sub$method <- relevel(factor(sub$method), ref = "WGCNA")
  fit <- fit_poisson_offset(sub, "method + combat")
  co <- summary(fit$model)$coefficients
  rr <- exp(co["methodhub_anchored", "Estimate"])
  ci <- exp(co["methodhub_anchored", "Estimate"] + c(-1.96, 1.96) * co["methodhub_anchored", "Std. Error"])
  pval <- co["methodhub_anchored", grepl("^Pr", colnames(co))]
  resultsB[[ref]] <- data.frame(
    reference = ref, family = fit$family_used,
    rate_ratio_hub_vs_WGCNA = rr, ci_low = ci[1], ci_high = ci[2], p_value = pval,
    stringsAsFactors = FALSE
  )
}
resultsB_df <- do.call(rbind, resultsB); rownames(resultsB_df) <- NULL
print(resultsB_df, digits = 3)

############################ C. r_cutoff trend test, hub-anchored only #########
message("\n=== C. Precision trend vs r_cutoff (hub-anchored method only), per reference ===")
resultsC <- list()
for (ref in references) {
  sub <- overlap[overlap$method == "hub_anchored" & overlap$reference == ref, ]
  fit <- fit_poisson_offset(sub, "r_cutoff + combat")
  co <- summary(fit$model)$coefficients
  slope <- co["r_cutoff", "Estimate"]
  rr_per_0.01 <- exp(slope * 0.01)  # rate ratio per +0.01 increase in r_cutoff
  pval <- co["r_cutoff", grepl("^Pr", colnames(co))]
  resultsC[[ref]] <- data.frame(
    reference = ref, family = fit$family_used,
    rate_ratio_per_0.01_rcutoff = rr_per_0.01, p_value_trend = pval,
    stringsAsFactors = FALSE
  )
}
resultsC_df <- do.call(rbind, resultsC); rownames(resultsC_df) <- NULL
print(resultsC_df, digits = 3)

output_dir <- "network_robustness/string_reference"
write.csv(resultsA_df, file.path(output_dir, "poisson_rateratio_hub_vs_wgcna_by_combat.csv"), row.names = FALSE)
write.csv(resultsB_df, file.path(output_dir, "poisson_rateratio_hub_vs_wgcna_pooled.csv"), row.names = FALSE)
write.csv(resultsC_df, file.path(output_dir, "poisson_precision_trend_vs_rcutoff.csv"), row.names = FALSE)
message("\nSaved: poisson_rateratio_hub_vs_wgcna_by_combat.csv, ",
        "poisson_rateratio_hub_vs_wgcna_pooled.csv, poisson_precision_trend_vs_rcutoff.csv")
