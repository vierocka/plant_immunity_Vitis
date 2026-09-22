#!/bin/bash -l
#SBATCH --cpus-per-task=12
#SBATCH --mem=48gb
#SBATCH --time=24:00:00
#SBATCH --output=sra_download_%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=youremail@here
#SBATCH --account=xxx

###############################################################################
# Download all 83 runs backing "Vitis rotundifolia Genes Introgressed with
# RUN1 and RPV1: Poor Recombination and Impact on V. vinifera Berry
# Transcriptome" (Plants, 2024) via sra-tools (prefetch + fasterq-dump), for
# reprocessing with this project's own STAR/featureCounts pipeline.
#
# 2026-08-14: originally parallelized with GNU parallel
# (tools/parallel/20240322-GCCcore-13.2.0), but EVERY SRA-Toolkit build on
# Ramses (3.0.10-gompi-2023a, 3.0.10-gompi-2023a-Java-8, 3.0.5-gompi-2022b)
# needs a GCCcore version (12.3.0 / 12.2.0) that Lmod refuses to load
# alongside parallel's GCCcore/13.2.0 -- confirmed directly, no combination
# of the two modules loads cleanly. Switched to `xargs -P` instead: it's a
# base coreutils/findutils tool, no module required, so it can't conflict
# with SRA-Toolkit's toolchain at all. Same 12-way concurrency as before.
#
# 12 concurrent fasterq-dump extractions need proportionally more scratch
# DISK headroom than one-at-a-time -- budget for several runs' temp files
# live simultaneously, not just one.
#
# Three BioProjects, confirmed via ENA portal API 2026-08-14:
#   PRJNA862686   Syrah, single-berry developmental time course, 33 runs
#   PRJNA1118503  G5 (microvine parental line), 8 runs
#   PRJNA1121615  MV102 (RUN1/RPV1 carrier) vs MV032 (non-carrier) during
#                 ripening, 42 runs (21 MV102 + 21 MV032)
#
# sra_accession_list.txt (this folder) = all 83 run accessions, one per line.
# all_runs_combined.tsv = same 83 runs with bioproject/sample_alias/library_
# name/experiment_title, for mapping SRR IDs back to condition/replicate.
###############################################################################

set -euo pipefail

module purge
module add bio/SRA-Toolkit/3.0.10-gompi-2023a-Java-8

WORK_DIR="/scratch/USERNAME/Vitis/RUN1_RPV1_microvine_2024"
ACC_LIST="${WORK_DIR}/sra_accession_list.txt"
FASTQ_DIR="${WORK_DIR}/fastq"
JOBLOG="${WORK_DIR}/xargs_joblog.tsv"

mkdir -p "$WORK_DIR" "$FASTQ_DIR"
cd "$WORK_DIR"

export WORK_DIR FASTQ_DIR JOBLOG

process_one() {
    local RUN="$1"
    local START
    START=$(date +%s)

    if [[ -f "${FASTQ_DIR}/${RUN}_1.fastq.gz" && -f "${FASTQ_DIR}/${RUN}_2.fastq.gz" ]]; then
        echo "SKIP ${RUN}: already downloaded"
        return 0
    fi
    # clear any stale partial output from an interrupted prior run (the
    # SRR29330621/622 case: old .gz sitting next to a fresh re-extraction)
    rm -f "${FASTQ_DIR}/${RUN}_1.fastq" "${FASTQ_DIR}/${RUN}_2.fastq" \
          "${FASTQ_DIR}/${RUN}_1.fastq.gz" "${FASTQ_DIR}/${RUN}_2.fastq.gz"

    echo "=== ${RUN} start $(date) ==="

    local STATUS=0
    prefetch --max-size 100G -O "$WORK_DIR" "$RUN" \
        && fasterq-dump \
               --split-files \
               --threads 1 \
               --outdir "$FASTQ_DIR" \
               "${WORK_DIR}/${RUN}/${RUN}.sra" \
        && gzip -f "${FASTQ_DIR}/${RUN}_1.fastq" "${FASTQ_DIR}/${RUN}_2.fastq" \
        || STATUS=$?

    rm -rf "${WORK_DIR}/${RUN}"

    # flock keeps the 12 concurrent workers from interleaving joblog lines
    { flock -x 200
      echo -e "${RUN}\t${STATUS}\t$(( $(date +%s) - START ))s" >> "$JOBLOG"
    } 200>>"${JOBLOG}.lock"

    echo "=== ${RUN} done $(date), exit ${STATUS} ==="
}
export -f process_one

echo -e "run_accession\texit_status\tseconds" > "$JOBLOG"

xargs -P "${SLURM_CPUS_PER_TASK}" -I{} bash -c 'process_one "$@"' _ {} < "$ACC_LIST"

echo "Done. Fastqs in ${FASTQ_DIR}. Joblog (check exit_status != 0 for failures): ${JOBLOG}"
