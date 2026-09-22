#!/bin/bash -l
#SBATCH --cpus-per-task=24
#SBATCH --mem=68gb
#SBATCH --time=2:00:00
#SBATCH --output=froussios2019_prep_reference.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=youremail@here
#SBATCH --account=xxx

###############################################################################
# Reference download + STAR index build ONLY, split out of
# froussios2019_star_featurecounts_ramses.sh so it can run in parallel with
# froussios2019_download_ramses.sh (the fastq download) rather than serially
# after it. Submitting the full STAR/featureCounts script before the fastq
# download finishes would build this same reference/index fine, then abort
# immediately (set -euo pipefail + explicit missing-fastq check) on the first
# sample whose fastq isn't downloaded yet -- wasting the SLURM allocation.
# This script is idempotent (same existence checks) and safe to run either
# before or after the main script -- whichever runs second just skips
# everything here since it's already done.
###############################################################################

set -euo pipefail

REF_DIR="/scratch/USERNAME/Athal/TAIR10ref"
GENOME_FA="${REF_DIR}/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
ANNOT_GFF3="${REF_DIR}/Arabidopsis_thaliana.TAIR10.63.gff3"
STAR_INDEX="${REF_DIR}/STAR_index_TAIR10"

mkdir -p "$REF_DIR"

############################ 0. REFERENCE (verified URLs, 2026-08-17) ##########
if [ ! -f "$GENOME_FA" ]; then
  curl -sS -o "${GENOME_FA}.gz" \
    "https://ftp.ebi.ac.uk/ensemblgenomes/pub/plants/release-63/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz"
  gunzip "${GENOME_FA}.gz"
fi
if [ ! -f "$ANNOT_GFF3" ]; then
  curl -sS -o "${ANNOT_GFF3}.gz" \
    "https://ftp.ebi.ac.uk/ensemblgenomes/pub/plants/release-63/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.63.gff3.gz"
  gunzip "${ANNOT_GFF3}.gz"
fi

# NEGATIVE-LOGIC CHECK (same as the main script): confirm gene_id attribute
# exists before featureCounts would need it downstream.
# NOTE (2026-08-18): originally piped through `| head -1`, which -- under
# `set -o pipefail` -- reports pipeline failure whenever `head` closes the
# pipe early and awk gets SIGPIPE while still writing the rest of the file,
# REGARDLESS of what grep found downstream. Confirmed on Ramses: the check
# fired "ASSUMPTION VIOLATED" while its own debug line printed a row that
# plainly contains gene_id=AT1G01010. Fixed by having awk stop itself after
# the first match (`{print; exit}`) instead of relying on `head` to truncate
# a still-writing pipe.
FIRST_GENE_LINE=$(awk -F'\t' '$3=="gene"{print; exit}' "$ANNOT_GFF3")
if [[ "$FIRST_GENE_LINE" != *gene_id* ]]; then
  echo "ASSUMPTION VIOLATED: no 'gene_id' attribute found on gene-level GFF3 records." >&2
  echo "$FIRST_GENE_LINE" >&2
  exit 1
fi

GENOME_LEN=$(grep -v '^>' "$GENOME_FA" | tr -d '\n' | wc -c)
echo "TAIR10 toplevel genome length: ${GENOME_LEN} bp"
SAINDEX_CHECK=$(python3 -c "import math; print(min(14, int(math.log2(${GENOME_LEN})/2 - 1)))")
echo "Formula-recommended --genomeSAindexNbases: ${SAINDEX_CHECK} (script uses 12 below -- confirm these match)"

############################ 1. STAR GENOME INDEX ###############################
module unload compiler/GCC 2>/dev/null || true
module add bio/STAR/2.7.10b-GCC-11.3.0

if [ ! -f "${STAR_INDEX}/SAindex" ]; then
  mkdir -p "$STAR_INDEX"
  STAR --runThreadN 24 \
       --runMode genomeGenerate \
       --genomeDir "$STAR_INDEX" \
       --genomeFastaFiles "$GENOME_FA" \
       --sjdbGTFfile "$ANNOT_GFF3" \
       --sjdbGTFtagExonParentTranscript Parent \
       --sjdbGTFfeatureExon exon \
       --sjdbOverhang 100 \
       --genomeSAindexNbases 12
fi

echo "Done. Reference + STAR index ready in ${REF_DIR}."
echo "froussios2019_star_featurecounts_ramses.sh will skip both steps above (already-present checks) and start directly at per-sample alignment once the fastq download finishes."
