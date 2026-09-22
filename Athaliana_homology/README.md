# Athaliana_homology

## Aim

Annotate protein-coding genes in *Vitis vinifera* with their *Arabidopsis
thaliana* homologs, since Arabidopsis is the reference plant species with
by far the most extensively studied plant immunity mechanisms — this gives
Vitis genes (including all network/GCNA hub and module genes) an
interpretable functional annotation.

Output feeds `data_files/26169genes_with_AthalHomologs_allIDs_exprPatterns_TAIR10ids.csv`,
used throughout `06_PCNWA/` for gene ID conversion and functional annotation.

## Scripts

Two independent homology approaches, run separately:

- `qVvin_refAthalProt_blP.sh` + `Athal_homologs.sh` — direct BLASTP
  of all Vitis PN40024 v4 proteins against the TAIR10 Arabidopsis proteome
  (NCBI RefSeq assembly GCF_000001735.3), best hit per gene with percent
  identity.
- `call_homologs_Athal_Vvit.sh` — orthologous-group (COG) based homology
  via STRING-db v12.0 mappings, using the Vitis (taxid 29760) and
  Arabidopsis (taxid 3702) COG assignments.

`Athal_tair_proteinID_homologsSearch.sh` is a narrower, older lookup against
a single-genotype (Rpv1 vs. wild-type cultivar) gene panel from an earlier
analysis iteration, using the STRING-db output above; kept for reference,
not part of the current 26,169-gene annotation pipeline.

These are standalone shell scripts intended to be run manually, referencing
a local reference-data folder outside this repository
(`~/Dropbox/MendelUni_Vinselect/reference/`).
