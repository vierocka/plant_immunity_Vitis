#!/bin/bash -l
#SBATCH --cpus-per-task=8
#SBATCH --mem=16gb
#SBATCH --time=24:00:00
#SBATCH --output=chitarrini2020_filtering_%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=youremail@here
#SBATCH --account=xxx

###############################################################################
# Read filtering for Chitarrini et al. 2020 (ENA PRJEB28042, 33 runs). One
# quality filter only: average per-read quality >= Q30 (fastp), every other
# fastp filter explicitly disabled -- STAR has no native raw-quality
# filter, so this is the minimal addition needed before alignment, not a
# trimming step.
###############################################################################

set -euo pipefail

WORK_DIR="/scratch/USERNAME/Vitis/Chitarrini2020_RNAseq"
READS="/scratch/USERNAME/Vitis/"
RUN_LIST="${WORK_DIR}/chitarrini_run_list.tsv"

cd "$READS"

module unload compiler/GCC 2>/dev/null || true
module add bio/fastp/1.0.1-GCC-13.3.0   # CHECK EXACT MODULE NAME/VERSION: `module spider fastp` first.

for ((i=2; i<35; i++)); do
  run=$(sed -n ''$i'p' "$RUN_LIST" | cut -f3 | cut -d"/" -f7 | cut -d"_" -f1)

  fastp \
    -i "${run}_1.fastq.gz" -I "${run}_2.fastq.gz" \
    -o "${run}_1.q30.fastq.gz" -O "${run}_2.q30.fastq.gz" \
    --average_qual 30 \
    --disable_length_filtering \
    -w 8 \
    -j "${run}_fastp.json" -h "${run}_fastp.html"
    # --disable_quality_filtering turns off fastp's OWN default per-base
    # quality filter -- --average_qual 30 is a distinct filter, applied
    # regardless of that flag.
done

echo "Done. Q30-filtered fastqs in ${READS}."
