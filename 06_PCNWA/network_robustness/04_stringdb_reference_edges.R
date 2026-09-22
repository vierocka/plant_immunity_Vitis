###############################################################################
# Build STRING-derived reference edge sets for benchmarking the hub-anchored
# co-transcriptional network against 3 independent "positive interaction"
# definitions (all 3 kept, rather than committing to just one):
#
#   ref_detailed_any  - any edge listed in 29760.protein.links.detailed.v12.0.
#                        filtered.txt (already pre-filtered to require >0
#                        evidence in at least one non-textmining channel:
#                        neighborhood/fusion/cooccurence/coexpression/
#                        experimental/database).
#   ref_detailed_700  - same file, combined_score > 700 (STRING "high
#                        confidence").
#   ref_physical_any  - any edge listed in
#                        29760.protein.physical.links.v12.0.txt (physical PPI
#                        evidence only - the closest thing to "ground truth"
#                        available for this non-model species).
#   ref_coexpression_above_median - detailed-filtered file, column 6
#                        ("coexpression" evidence channel) strictly above its
#                        own median (50, computed across all 8,319,580
#                        pre-filtered rows incl. zeros). IMPORTANT CAVEAT:
#                        STRING's coexpression channel is mined
#                        from public expression compendia spanning many
#                        unrelated conditions/tissues/experiments, and for a
#                        less-studied species like V. vinifera may partly be
#                        transferred via orthology from other plants. A high
#                        score therefore reflects "frequently co-expressed
#                        across many, mostly UNRELATED conditions in
#                        general", not "co-expressed specifically during
#                        downy mildew infection". This makes it the reference
#                        set conceptually closest to what our own network
#                        measures - so overlap here should be read as a
#                        SOFTER, not stronger, validation than the physical-
#                        PPI reference: it can't distinguish "true infection-
#                        specific co-regulation STRING happens to have also
#                        seen elsewhere" from "coincidental general
#                        co-expression unrelated to immunity", and real
#                        infection-specific co-regulation with no cross-
#                        condition STRING signal will show as a false
#                        negative here, not a method flaw.
#
# IMPORTANT CAVEATS (documented here, repeated in the README):
#  1. STRING is PPI/functional-association evidence, not co-expression
#     evidence. A co-transcriptional module edge that is absent from STRING
#     is NOT necessarily a false positive - STRING coverage for a non-model
#     crop is incomplete, and co-expression captures real regulatory
#     relationships (e.g. shared transcription factor, same pathway) that
#     are not physical interactions at all. These reference sets bound
#     "how much of my network coincides with independently-derived
#     interaction evidence", not "how many of my edges are wrong".
#  2. ID mapping (06_PCNWA/string-db/results/26169genes_conversions_gene_protein_IDs.tsv) is
#     incomplete: 1804 / 26169 genes (6.9%) have no UniProt ID and are
#     therefore INVISIBLE to STRING comparison entirely - any edge touching
#     them can't be evaluated and is excluded (not counted as
#     positive or negative).
#  3. 4057 UniProt IDs are shared by >1 Vitvi gene model (paralogs/isoforms
#     collapsed onto one STRING protein). A pair of Vitvi genes that BOTH
#     map to the same UniProt ID is ambiguous (STRING cannot distinguish
#     them) and is excluded from evaluation, flagged in ambiguous_string_ids.
#
# STRING files list every edge in BOTH directions (verified: exactly 2x per
# canonical pair) - canonicalized here via pmin/pmax pair keys.
###############################################################################

suppressPackageStartupMessages(library(data.table))

output_dir <- "network_robustness/string_reference"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
stringdb_dir <- "string-db/raw_data"

############################ 1. GENE <-> STRING ID MAPPING #####################
conv <- read.table(file.path("string-db/results/26169genes_conversions_gene_protein_IDs.tsv"),
                    sep = "\t", header = TRUE, quote = "", comment.char = "")
conv <- conv[!is.na(conv$Uniprot_proteinID) & conv$Uniprot_proteinID != "", ]
conv$string_id <- paste0("29760.", conv$Uniprot_proteinID)

ambiguous_string_ids <- unique(conv$string_id[duplicated(conv$string_id) | duplicated(conv$string_id, fromLast = TRUE)])
message(nrow(conv), " genes have a UniProt ID; ", length(ambiguous_string_ids),
        " STRING IDs are shared by >1 gene (ambiguous, will be excluded from edge evaluation).")

gene_to_string <- setNames(conv$string_id, conv$ENSMBL_geneID)
# Drop ambiguous mappings from the lookup used for edge evaluation - a Vitvi
# gene whose STRING ID is shared with another gene cannot be unambiguously
# tested against STRING at the gene-pair level.
gene_to_string_unambiguous <- gene_to_string[!gene_to_string %in% ambiguous_string_ids]
message(length(gene_to_string_unambiguous), " / ", length(gene_to_string),
        " gene->STRING mappings are unambiguous and usable for edge evaluation.")

pair_key <- function(a, b) paste(pmin(a, b), pmax(a, b), sep = "||")

#' Map a set of Vitvi gene-gene pairs (two character vectors, same length) to
#' STRING pair keys, dropping pairs where either gene lacks an unambiguous
#' STRING mapping.
#' @return character vector of pair keys (may be shorter than input).
genes_to_string_pairkeys <- function(gene_a, gene_b) {
  sa <- gene_to_string_unambiguous[gene_a]
  sb <- gene_to_string_unambiguous[gene_b]
  keep <- !is.na(sa) & !is.na(sb)
  if (!any(keep)) return(character(0))
  pair_key(sa[keep], sb[keep])
}

############################ 2. LOAD STRING REFERENCE FILES ####################
message("Loading STRING physical links...")
physical <- fread(file.path(stringdb_dir, "29760.protein.physical.links.v12.0.txt"),
                   sep = " ", header = TRUE)
message("Loading STRING detailed (pre-filtered, non-textmining-only) links...")
detailed <- fread(file.path(stringdb_dir, "29760.protein.links.detailed.v12.0.filtered.txt"),
                   sep = " ", header = FALSE,
                   col.names = c("protein1", "protein2", "neighborhood", "fusion",
                                 "cooccurence", "coexpression", "experimental",
                                 "database", "textmining", "combined_score"))

canonicalize <- function(dt) {
  dt[, pk := pair_key(protein1, protein2)]
  unique(dt, by = "pk")
}
physical <- canonicalize(physical)
detailed <- canonicalize(detailed)
message("Canonical (undirected, deduplicated) edges: physical=", nrow(physical),
        ", detailed-filtered=", nrow(detailed))

############################ 3. BUILD 3 POSITIVE-EDGE REFERENCE SETS ###########
ref_physical_any <- physical$pk
ref_detailed_any <- detailed$pk
ref_detailed_700 <- detailed$pk[detailed$combined_score > 700]
coexpression_median <- median(detailed$coexpression)
ref_coexpression_above_median <- detailed$pk[detailed$coexpression > coexpression_median]

message("Reference edge-set sizes: physical_any=", length(ref_physical_any),
        ", detailed_any=", length(ref_detailed_any),
        ", detailed_700=", length(ref_detailed_700),
        ", coexpression_above_median (>", coexpression_median, ")=",
        length(ref_coexpression_above_median))

saveRDS(list(
  ref_physical_any = ref_physical_any,
  ref_detailed_any = ref_detailed_any,
  ref_detailed_700 = ref_detailed_700,
  ref_coexpression_above_median = ref_coexpression_above_median
), file.path(output_dir, "string_reference_pairkeys.rds"))

saveRDS(gene_to_string_unambiguous, file.path(output_dir, "gene_to_string_unambiguous.rds"))
saveRDS(ambiguous_string_ids, file.path(output_dir, "ambiguous_string_ids.rds"))
# also keep combined_score per detailed pair for ROC sweeps that want to plot
# STRING confidence rather than a fixed threshold.
saveRDS(detailed[, .(pk, combined_score)], file.path(output_dir, "detailed_pk_scores.rds"))

############################ 4. STRUCTURAL VALIDATION CHECK ####################
# Round-trip sanity check: take genes that DO have an unambiguous STRING
# mapping, restrict the raw STRING file to pairs where BOTH endpoints are in
# that mapped universe, and confirm every such row is found by our
# genes_to_string_pairkeys() + reference-set lookup (must be true by
# construction, but a mismatch here would reveal a pair-key or mapping bug).
mapped_string_ids <- unique(gene_to_string_unambiguous)
string_to_gene <- setNames(names(gene_to_string_unambiguous), gene_to_string_unambiguous)

sample_rows <- physical[protein1 %in% mapped_string_ids & protein2 %in% mapped_string_ids][sample(.N, min(200, .N))]
sample_gene_a <- string_to_gene[sample_rows$protein1]
sample_gene_b <- string_to_gene[sample_rows$protein2]
recovered_pk <- genes_to_string_pairkeys(sample_gene_a, sample_gene_b)
check1 <- all(recovered_pk %in% ref_physical_any)
message("Validation check 1 (round-trip: mapped gene pairs from a random STRING-physical ",
        "sample are found in ref_physical_any): ", ifelse(check1, "PASS", "FAIL"))
if (!check1) stop("STRING reference-edge round-trip validation FAILED - do not trust downstream ROC/Jaccard results.")

# Second check: a pair of genes known NOT to be in the physical file (two
# genes chosen at random, extremely unlikely to collide with a real STRING
# physical edge given ~459k canonical physical edges vs C(24000,2) possible
# pairs) should correctly evaluate as absent.
set.seed(1)
random_genes <- sample(names(gene_to_string_unambiguous), 4)
neg_pk <- genes_to_string_pairkeys(random_genes[1:2], random_genes[3:4])
check2 <- length(neg_pk) == 0 || !any(neg_pk %in% ref_physical_any)
message("Validation check 2 (random unrelated gene pair correctly absent from ref_physical_any): ",
        ifelse(check2, "PASS", "FAIL (extremely unlikely collision, or a bug - investigate)"))

writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
message("Completed STRING reference-edge build. Outputs in: ", output_dir)
