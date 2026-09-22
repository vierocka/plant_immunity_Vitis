#!/bin/bash -l
#SBATCH --cpus-per-task=8
#SBATCH --mem=16gb
#SBATCH --time=24:00:00
#SBATCH --output=froussios2019_filtering_%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=youremail@here
#SBATCH --account=xxx

###############################################################################
# Read filtering for Froussios et al. 2019 (ArrayExpress E-MTAB-5446 / ENA
# ERP021226, 17 runs). One quality filter only: average per-read quality
# >= Q30 (fastp), every other fastp filter explicitly disabled -- STAR has
# no native raw-quality filter, so this is the minimal addition needed
# before alignment, not a trimming step.
###############################################################################

set -euo pipefail

WORK_DIR="/scratch/USERNAME/Athal"
FASTQ_DIR="${WORK_DIR}/fastq"                        # from froussios2019_download_ramses.sh
RUN_LIST="${WORK_DIR}/froussios2019_run_list.tsv"

cd "$WORK_DIR"

module unload compiler/GCC 2>/dev/null || true
module add bio/fastp/1.0.1-GCC-13.3.0   # CHECK EXACT MODULE NAME/VERSION: `module spider fastp` first.

tail -n +2 "$RUN_LIST" | while IFS=$'\t' read -r run sample exp mix url1 url2 md5_1 md5_2; do
  sample_dir="${WORK_DIR}/${run}_${sample}"
  mkdir -p "$sample_dir"
  cd "$sample_dir"

  [ -f "${FASTQ_DIR}/${run}_1.fastq.gz" ] || { echo "Missing ${FASTQ_DIR}/${run}_1.fastq.gz -- run froussios2019_download_ramses.sh first" >&2; exit 1; }
  [ -f "${FASTQ_DIR}/${run}_2.fastq.gz" ] || { echo "Missing ${FASTQ_DIR}/${run}_2.fastq.gz -- run froussios2019_download_ramses.sh first" >&2; exit 1; }
  [ -f "${run}_1.fastq.gz" ] || ln -s "${FASTQ_DIR}/${run}_1.fastq.gz" "${run}_1.fastq.gz"
  [ -f "${run}_2.fastq.gz" ] || ln -s "${FASTQ_DIR}/${run}_2.fastq.gz" "${run}_2.fastq.gz"

  if [ ! -f "${run}_1.q30.fastq.gz" ]; then
    fastp \
      -i "${run}_1.fastq.gz" -I "${run}_2.fastq.gz" \
      -o "${run}_1.q30.fastq.gz" -O "${run}_2.q30.fastq.gz" \
      --average_qual 30 \
      --disable_adapter_trimming \
      --disable_length_filtering \
      --disable_quality_filtering \
      -w 8 \
      -j "${run}_fastp.json" -h "${run}_fastp.html"
  fi
  cd "$WORK_DIR"
done

echo "Done. Q30-filtered fastqs in each \${WORK_DIR}/<run>_<sample>/ subfolder."
