###############################################################################
# T-index cross-cultivar match validation, MV102 (RUN1/RPV1 carrier) vs MV32
# (non-carrier), Shi et al. 2024 (Plants 13:2095) reprocessing.
#
# WHY THIS SCRIPT EXISTS: the T0-T6 labels attached to each sample are
# per-genotype developmental/physiological stage indices (paper's own words,
# Materials & Methods 4.1/4.4: staging is "guided by considerations of
# relative growth, sugars, and organic acids"), NOT calendar-synchronized
# sampling dates. Assuming "MV102-T3" and "MV32-T3" are the same physiological
# ripening stage just because they share a label is an unverified assumption,
# not a fact -- exactly the kind of thing to test rather than trust. This
# script tests it directly against Table S2 (single-berry malic acid,
# tartaric acid, and glucose+fructose/"sugar" HPLC measurements per sample),
# using sugar (G+F, mmol/L) as the primary staging metric (matches the
# paper's own staging criterion) and malic acid as an independent
# cross-check (malic acid degrades monotonically during ripening, an
# orthogonal biochemical clock).
#
# LOGIC (deliberately adversarial to the "T-index = T-index" default
# assumption, per instruction to use negative logic):
#   1. Null/default assumption to disprove: nominal T-index pairing is no
#      better than a random pairing of MV102 stages to MV32 stages.
#   2. Global test: brute-force all 7! = 5040 possible 1-to-1 alignments of
#      MV102's 7 stages to MV32's 7 stages; rank the nominal (identity)
#      alignment's total squared log-sugar-difference among all 5040.
#   3. Per-pair test (what actually matters for deciding which T-pairs are
#      usable in the cross-cultivar AED script): for each MV102 T-index,
#      find its nearest MV32 stage in log-sugar space. A nominal pair
#      MV102-Ti/MV32-Ti is only called VALIDATED if (a) MV32-Ti is in fact
#      the nearest neighbor to MV102-Ti (not some other MV32 stage), AND
#      (b) the relative sugar difference is below a pre-set 30% tolerance.
#      Anything failing either test is flagged NOT VALIDATED and excluded
#      from the cross-cultivar AED script by default.
#   4. Malic acid is used only as a corroboration check on the sugar-based
#      call, not as an independent veto -- reported alongside, not fused
#      into a combined score (avoids post-hoc metric-shopping).
###############################################################################

suppressPackageStartupMessages(library(dplyr))

script_dir <- "."
supp_dir <- file.path(script_dir, "plants-13-02095_supplements")
output_dir <- file.path(script_dir, "AED_check_results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

traits <- read.csv(file.path(supp_dir, "Table S2_Single berry traits.csv"), stringsAsFactors = FALSE)

# Keep only MV102/MV32, parse genotype + T-index + replicate from Sample.
mv <- traits[traits$Genotype %in% c("MV102", "MV32"), ]
mv$T_index <- as.integer(sub(".*T([0-9]+)_.*", "\\1", mv$Sample))
mv$replicate <- as.integer(sub(".*_([0-9]+)$", "\\1", mv$Sample))
stopifnot(!anyNA(mv$T_index), !anyNA(mv$replicate))

# Sanity check: expect T0-T6 (7 stages) x up to 3 reps per genotype, nothing
# outside that range -- if this fails, the regex above is wrong, not the data.
stopifnot(all(mv$T_index %in% 0:6))

per_stage <- mv %>%
  group_by(Genotype, T_index) %>%
  summarise(
    n_samples = n(),
    mean_sugar_GF = mean(G.F..mmol.L.),
    mean_malic = mean(Malic.acid..mEq.L.),
    .groups = "drop"
  )

mv102 <- per_stage[per_stage$Genotype == "MV102", ]
mv32 <- per_stage[per_stage$Genotype == "MV32", ]
mv102 <- mv102[order(mv102$T_index), ]
mv32 <- mv32[order(mv32$T_index), ]
stopifnot(nrow(mv102) == 7, nrow(mv32) == 7)  # T0-T6 both sides, no gaps

cat("=== Per-stage sugar (G+F) and malic acid means ===\n")
cat("MV102:\n"); print(as.data.frame(mv102))
cat("MV32:\n");  print(as.data.frame(mv32))

############################ GLOBAL PERMUTATION TEST ###########################
# Null hypothesis to knock down: nominal (identity) T-index alignment is not
# meaningfully better, in total squared log-sugar difference, than a random
# alignment of MV102's 7 stages to MV32's 7 stages.
log_sugar_102 <- log(mv102$mean_sugar_GF)
log_sugar_32 <- log(mv32$mean_sugar_GF)

# all_perms elements are permutations of POSITIONS 1:7 (log_sugar_32 is
# already ordered by T_index 0..6, so position i <-> T_index i-1). Keeping
# this 1-based-position-only (no value/index mixing) to avoid off-by-one bugs.
all_perms <- combinat::permn(1:7)
sse_for_perm <- function(perm_positions) sum((log_sugar_102 - log_sugar_32[perm_positions])^2)
all_sse <- vapply(all_perms, sse_for_perm, numeric(1))
identity_sse <- sse_for_perm(1:7)
identity_rank <- sum(all_sse <= identity_sse)  # 1 = best possible alignment
identity_percentile <- mean(all_sse <= identity_sse)
best_perm_positions <- all_perms[[which.min(all_sse)]]
best_perm_Tindex <- best_perm_positions - 1L  # position -> T-index label

cat("\n=== Global permutation test (n=5040 possible 1:1 alignments) ===\n")
cat(sprintf(
  "Identity (nominal T=T) alignment SSE=%.4f, rank %d/5040 (percentile %.4f -- lower is better)\n",
  identity_sse, identity_rank, identity_percentile
))
cat("Best-fitting alignment found (MV102 T-index -> MV32 T-index):\n")
print(data.frame(MV102_T = 0:6, best_match_MV32_T = best_perm_Tindex))
cat(sprintf("Best possible SSE=%.4f (vs identity's %.4f)\n", min(all_sse), identity_sse))

############################ PER-PAIR NEAREST-NEIGHBOR VALIDATION ##############
rel_diff <- function(a, b) abs(a - b) / mean(c(a, b))

nn_check <- do.call(rbind, lapply(seq_len(7), function(i) {
  t_idx <- mv102$T_index[i]
  sugar_102 <- mv102$mean_sugar_GF[i]
  diffs_sugar <- abs(log(sugar_102) - log(mv32$mean_sugar_GF))
  nn_sugar_T <- mv32$T_index[which.min(diffs_sugar)]
  nominal_is_nn_sugar <- (nn_sugar_T == t_idx)
  reldiff_sugar_nominal <- rel_diff(sugar_102, mv32$mean_sugar_GF[mv32$T_index == t_idx])

  malic_102 <- mv102$mean_malic[i]
  diffs_malic <- abs(log(malic_102) - log(mv32$mean_malic))
  nn_malic_T <- mv32$T_index[which.min(diffs_malic)]
  nominal_is_nn_malic <- (nn_malic_T == t_idx)

  validated <- nominal_is_nn_sugar && (reldiff_sugar_nominal < 0.30)

  data.frame(
    T_index = t_idx,
    MV102_sugar = sugar_102, MV32_sugar = mv32$mean_sugar_GF[mv32$T_index == t_idx],
    nearest_MV32_T_by_sugar = nn_sugar_T,
    nominal_is_nearest_by_sugar = nominal_is_nn_sugar,
    relative_sugar_diff_at_nominal_pair = round(reldiff_sugar_nominal, 3),
    nearest_MV32_T_by_malic = nn_malic_T,
    nominal_is_nearest_by_malic = nominal_is_nn_malic,
    VALIDATED_for_crosscultivar_AED = validated,
    stringsAsFactors = FALSE
  )
}))

cat("\n=== Per-T-index validation (primary=sugar nearest-neighbor + <30% rel diff; malic acid = corroboration only) ===\n")
print(nn_check)

n_validated <- sum(nn_check$VALIDATED_for_crosscultivar_AED)
cat(sprintf(
  "\n%d/7 nominal T-index pairs VALIDATED for cross-cultivar AED: %s\n",
  n_validated, paste(nn_check$T_index[nn_check$VALIDATED_for_crosscultivar_AED], collapse = ", ")
))
cat(sprintf(
  "%d/7 REJECTED (physiologically desynchronized at the same nominal label): %s\n",
  7 - n_validated, paste(nn_check$T_index[!nn_check$VALIDATED_for_crosscultivar_AED], collapse = ", ")
))

write.csv(nn_check, file.path(output_dir, "Tindex_crosscultivar_validation_MV102_vs_MV32.csv"), row.names = FALSE)
writeLines(
  as.character(nn_check$T_index[nn_check$VALIDATED_for_crosscultivar_AED]),
  file.path(output_dir, "Tindex_validated_for_crosscultivar_AED.txt")
)

message("\nDone. Validation table + validated-T-index list written to: ", output_dir)
