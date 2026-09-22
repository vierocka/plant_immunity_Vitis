#!/bin/bash -l
#SBATCH --cpus-per-task=24
#SBATCH --mem=64gb
#SBATCH --time=6:00:00
#SBATCH --output=chitarrini2020_counting_%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=youremail@here
#SBATCH --account=xxx

###############################################################################
# Gene-level read counting for Chitarrini et al. 2020 (33 runs), one combined
# featureCounts call across all sample BAMs -> one count matrix (columns =
# samples). Run after chitarrini2020_mapping_ramses.sh.
###############################################################################

set -euo pipefail

REF_DIR="/scratch/USERNAME/Vitis/GGref"
WORK_DIR="/scratch/USERNAME/Vitis/Chitarrini2020_RNAseq"
READS="/scratch/USERNAME/Vitis/"
RUN_LIST="${WORK_DIR}/chitarrini_run_list.tsv"
ANNOT_GFF3="${REF_DIR}/Vitis_vinifera.PN40024.v4.56.gff3"

cd "$READS"

module unload compiler/GCC 2>/dev/null || true
module add bio/Subread/2.1.1-GCC-13.2.0

BAMS=$(tail -n +2 "$RUN_LIST" | while IFS=$'\t' read -r run sample url1 url2; do
  echo -n "${run}_Aligned.sortedByCoord.out.bam "
done)

featureCounts \
  -T 24 \
  -p -B -C \
  -t gene \
  -g gene_id \
  -a "$ANNOT_GFF3" \
  -o "${WORK_DIR}/chitarrini2020_all_samples.counts.tsv" \
  $BAMS

echo "Done. Combined gene-level count matrix: ${WORK_DIR}/chitarrini2020_all_samples.counts.tsv"
