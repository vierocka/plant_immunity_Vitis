#!/bin/bash -l
#SBATCH --cpus-per-task=32
#SBATCH --mem=64gb
#SBATCH --time=6:00:00
#SBATCH --output=froussios2019_counting_%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=youremail@here
#SBATCH --account=xxx

###############################################################################
# Gene-level read counting for Froussios et al. 2019 (17 runs), one combined
# featureCounts call across all sample BAMs. Run after
# froussios2019_mapping_ramses.sh.
###############################################################################

set -euo pipefail

REF_DIR="/scratch/USERNAME/Athal/TAIR10ref"
WORK_DIR="/scratch/USERNAME/Athal"
RUN_LIST="${WORK_DIR}/froussios2019_run_list.tsv"
ANNOT_GFF3="${REF_DIR}/Arabidopsis_thaliana.TAIR10.63.gff3"

cd "$WORK_DIR"

module unload compiler/GCC 2>/dev/null || true
module add bio/Subread/2.1.1-GCC-13.2.0

BAMS=$(tail -n +2 "$RUN_LIST" | while IFS=$'\t' read -r run sample exp mix url1 url2 md5_1 md5_2; do
  echo -n "${WORK_DIR}/${run}_${sample}/${run}_Aligned.sortedByCoord.out.bam "
done)

featureCounts \
  -T 32 \
  -p -B -C \
  -t gene \
  -g gene_id \
  -a "$ANNOT_GFF3" \
  -o "${WORK_DIR}/froussios2019_samples_from1_to14_perGene.counts.tsv" \
  $BAMS

echo "Done. Combined gene-level count matrix: ${WORK_DIR}/froussios2019_samples_from1_to14_perGene.counts.tsv"
