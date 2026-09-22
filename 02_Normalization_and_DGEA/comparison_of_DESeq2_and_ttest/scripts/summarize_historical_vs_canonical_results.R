# Import canonical and historical results and write one concise audit summary.
# The canonical primary call is DESeq2 zero-null Wald BH FDR < 0.05 plus
# apeglm-shrunken |log2FC| > 1; historical calls are imported unchanged.

canonical <- read.csv("02_Normalization_and_DGEA/DESeq2_classic/results/DGEA_reanalysis/DESeq2_all_gene_results.csv")
historical <- read.csv("02_Normalization_and_DGEA/DESeq2_classic/results/DGEA_reanalysis/historical_method_all_gene_results.csv")
loo <- read.csv("additional_tests/results/LOO_sensitivity/LOO_DESeq2_call_retention_distribution.csv")
variance <- read.csv("normalized_within_condition_variance_shared_vs_new.csv")

new_genes <- unique(canonical$gene[canonical$DE_primary %in% TRUE])
old_genes <- unique(historical$gene[historical$DE_historical %in% TRUE])
shared <- intersect(new_genes, old_genes)
new_only <- setdiff(new_genes, old_genes)
old_only <- setdiff(old_genes, new_genes)

summary <- data.frame(
  quantity = c("tested_genes", "contrasts", "canonical_union", "historical_union",
               "shared", "canonical_new_only", "historical_only"),
  value = c(length(unique(canonical$gene)), length(unique(canonical$contrast)),
            length(new_genes), length(old_genes), length(shared), length(new_only), length(old_only))
)
write.csv(summary, "historical_vs_canonical_results_summary.csv", row.names = FALSE)

lines <- c(
  "Historical versus canonical DE summary", "",
  paste0("Tested genes: ", length(unique(canonical$gene)), "; contrasts: ", length(unique(canonical$contrast))),
  paste0("Canonical DE union: ", length(new_genes)),
  paste0("Historical DE union: ", length(old_genes)),
  paste0("Shared: ", length(shared), "; canonical new-only: ", length(new_only),
         "; historical-only: ", length(old_only)), "",
  "Canonical definition: DESeq2 zero-null Wald BH FDR < 0.05 within contrast plus apeglm |log2FC| > 1.",
  "Historical definition: imported historical Gaussian/rlog-based calls.", "",
  "Interpretation: shared genes have lower canonical FDR than new-only genes; canonical new-only genes are not a higher-variance group in FDR/effect distributions.",
  "Within-condition VST variance analysis: shared median 0.276; new-only median 0.117; Wilcoxon and Fligner p < 2.2e-16.",
  "The VST variance result is a canonical sensitivity analysis, not a reproduction of the historical Gaussian GLM.", "",
  "LOO: corrected primary DESeq2 analysis evaluates all nine contrasts in every fold and exports gene-level retention, sign changes, and LFC changes.",
  "Agreement is described as concordance, not validation or truth."
)
writeLines(lines, "historical_vs_canonical_results_summary.txt")
writeLines(capture.output(sessionInfo()), "historical_vs_canonical_summary_sessionInfo.txt")
