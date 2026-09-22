#!/bin/bash -l
#SBATCH --cpus-per-task=24
#SBATCH --mem=64gb
#SBATCH --time=6:00:00
#SBATCH --output=chitarrini2020_featurecounts_only_%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=youremail@here
#SBATCH --account=xxx

###############################################################################
# Standalone featureCounts-only run for Chitarrini2020 (2026-09-07). Trimming
# and STAR mapping already completed for 29 of 33 samples (user confirmed via
# `ls trimmomatic/ | grep -c _2P.fq.gz` = 29, "I am fine with 29 samples" -
# not chasing the remaining 4 down further). This script does ONLY the
# featureCounts step, on whatever BAMs actually exist - no re-running of
# trimming/mapping, no genome-index step, no ulimit fix needed (not relevant
# to featureCounts). Picks up paired and F-only BAMs; R-only will be entirely
# absent for this run (the R-only STAR branch was added to
# chitarrini2020_trimmomatic_star_fc_ramses.sh AFTER this job was already
# submitted/running - see project_bmc_bioinformatics_batch_simulation_paper.md
# 2026-09-04 note) and will just be skipped with a warning per sample, exactly
# as the existence-guard in the full script already handles - not an error.
###############################################################################

set -euo pipefail

REF_DIR="/scratch/USERNAME/Vitis/GGref"
WORK_DIR="/scratch/USERNAME/Vitis/Chitarrini2020"
RUN_LIST="/scratch/USERNAME/Vitis/Chitarrini2020/chitarrini_run_list.tsv"
ANNOT_GFF3="${REF_DIR}/Vitis_vinifera.PN40024.v4.56.gff3"

STAR_PAIRED_DIR="${WORK_DIR}/STARmap1"
STAR_FONLY_DIR="${WORK_DIR}/STARmap1_F"
STAR_RONLY_DIR="${WORK_DIR}/STARmap1_R"

mkdir -p "${STAR_PAIRED_DIR}/fC" "${STAR_FONLY_DIR}/fC" "${STAR_RONLY_DIR}/fC"
cd "$WORK_DIR"

module unload compiler/GCC 2>/dev/null || true
module add bio/Subread/2.1.1-GCC-13.2.0

tail -n +2 "$RUN_LIST" | while IFS=$'\t' read -r run sample url1 url2; do
  if [ ! -f "${STAR_PAIRED_DIR}/${run}_Aligned.sortedByCoord.out.bam" ]; then
    echo "SKIP paired featureCounts for ${run}: no BAM (mapping not completed for this run)." >&2
  elif [ ! -f "${STAR_PAIRED_DIR}/fC/${run}.count" ]; then
    featureCounts -T 24 -p -a "$ANNOT_GFF3" -t gene -g gene_id \
      -o "${STAR_PAIRED_DIR}/fC/${run}.count" \
      "${STAR_PAIRED_DIR}/${run}_Aligned.sortedByCoord.out.bam"
  fi
  if [ ! -f "${STAR_FONLY_DIR}/${run}_Aligned.sortedByCoord.out.bam" ]; then
    echo "SKIP F-only featureCounts for ${run}: no BAM (mapping not completed for this run)." >&2
  elif [ ! -f "${STAR_FONLY_DIR}/fC/${run}.count" ]; then
    featureCounts -T 24 -a "$ANNOT_GFF3" -t gene -g gene_id \
      -o "${STAR_FONLY_DIR}/fC/${run}.count" \
      "${STAR_FONLY_DIR}/${run}_Aligned.sortedByCoord.out.bam"
  fi
  if [ ! -f "${STAR_RONLY_DIR}/${run}_Aligned.sortedByCoord.out.bam" ]; then
    echo "SKIP R-only featureCounts for ${run}: no BAM (R-only branch not run for this job - see script header)." >&2
  elif [ ! -f "${STAR_RONLY_DIR}/fC/${run}.count" ]; then
    featureCounts -T 24 -a "$ANNOT_GFF3" -t gene -g gene_id \
      -o "${STAR_RONLY_DIR}/fC/${run}.count" \
      "${STAR_RONLY_DIR}/${run}_Aligned.sortedByCoord.out.bam"
  fi
done

N_TOTAL=$(( $(wc -l < "$RUN_LIST") - 1 ))
N_PAIRED=$(ls "${STAR_PAIRED_DIR}/fC/"*.count 2>/dev/null | wc -l)
N_FONLY=$(ls "${STAR_FONLY_DIR}/fC/"*.count 2>/dev/null | wc -l)
N_RONLY=$(ls "${STAR_RONLY_DIR}/fC/"*.count 2>/dev/null | wc -l)
echo "Of ${N_TOTAL} total runs in the run list: paired counts=${N_PAIRED}, F-only counts=${N_FONLY}, R-only counts=${N_RONLY}."
echo "Done. Paired counts: ${STAR_PAIRED_DIR}/fC/. F-only counts: ${STAR_FONLY_DIR}/fC/. R-only counts: ${STAR_RONLY_DIR}/fC/."
