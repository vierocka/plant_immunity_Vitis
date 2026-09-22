#!/bin/bash -l
#SBATCH --cpus-per-task=8
#SBATCH --mem=16gb
#SBATCH --time=24:00:00
#SBATCH --output=froussios2019_download_%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=youremail@here
#SBATCH --account=xxx
# ^ same caveat as the other RAMSES scripts in this project: verify the actual
#   account/partition name for this cluster before submitting.

###############################################################################
# Download all 17 runs backing Froussios, Schurch, Mackinnon, Gierlinski, Duc,
# Simpson, Barton (2019), "How well do RNA-Seq differential gene expression
# tools perform in a complex eukaryote? A case study in Arabidopsis thaliana,"
# Bioinformatics 35(18):3372-3377 (doi:10.1093/bioinformatics/btz089).
# ArrayExpress E-MTAB-5446 / ENA study ERP021226.
#
# Independent isogenic-background test for the 03_AED "12% of genes explain
# 75% of divergence" question: 17 WT Col-0 seedlings (14 days old), same lab,
# same growth conditions, RNA-Seq -- used the SAME WAY as the 03_AED noise
# negative control (real, same-genotype replicate variation with NO cross-
# genotype/cross-background confound), but from real biology instead of
# synthetic Gaussian jitter.
#
# VERIFIED DIRECTLY against the ENA/BioStudies APIs on 2026-08-17 (not taken
# from the paper text alone) -- see Froussios_2019/NOTES.md for the full
# provenance trail and exact queries used:
#   - Organism: Arabidopsis thaliana; Ecotype: Col-0; Genotype: wild type.
#   - 17 runs, paired-end, Illumina HiSeq 2000, ~77-130M read pairs/sample
#     (deep), read length 2x101bp (base_count/read_count/2 = ~101).
#   - IMPORTANT CAVEAT found during verification, not mentioned in the
#     05_noise_negative_control design discussion: the 17 samples are NOT one
#     flat replicate set. They split into 3 sequencing batches ("Experiment"
#     1/2/3, n=7/7/3) further crossed with 2 ERCC spike-in mixes within each
#     batch (this is a DE-tool-benchmarking dataset built around a known-
#     ground-truth spike-in design, not a pure replicate-variability study).
#     The spike-in mix affects only synthetic ERCC transcripts, not
#     endogenous Arabidopsis genes, and our STAR/featureCounts run below
#     targets the TAIR10 genome+annotation only (no ERCC sequences added to
#     the reference), so spike-in reads simply won't be counted -- endogenous
#     gene counts across all 17 samples remain a valid same-genotype,
#     same-condition replicate set. Kept as explicit columns
#     (experiment_batch, ercc_mix) in froussios2019_run_list.tsv regardless,
#     in case a batch-structure check is wanted later (mirrors this project's
#     habit of treating batch as a first-class column, cf. Batch1/BatchOrigin
#     in the Vitis scripts).
#
# WHY DIRECT ENA DOWNLOAD, NOT sra_download_ramses.sh's prefetch+fasterq-dump:
# these are ERR (ENA-native) accessions with fastq already hosted directly on
# the ENA FTP/HTTPS mirror, published together with per-file MD5 checksums --
# re-deriving fastq via the SRA-Toolkit round-trip (SRA->fastq) would be
# slower and adds a re-encoding step for no benefit here. Direct download also
# gives free integrity verification against ENA's own MD5s (the "negative
# logic" check requested: don't just assume a completed download is correct
# -- prove it against an independent, pre-published checksum).
###############################################################################

set -euo pipefail

WORK_DIR="/scratch/USERNAME/Athal"
RUN_LIST="${WORK_DIR}/froussios2019_run_list.tsv"   # run_accession, sample_title, experiment_batch, ercc_mix, url_R1, url_R2, md5_R1, md5_R2
FASTQ_DIR="${WORK_DIR}/fastq"
JOBLOG="${WORK_DIR}/download_joblog.tsv"

mkdir -p "$WORK_DIR" "$FASTQ_DIR"
cd "$WORK_DIR"
# Copy froussios2019_run_list.tsv here before submitting (scp'd alongside this script).

export FASTQ_DIR JOBLOG

download_one() {
    local RUN="$1" SAMPLE="$2" URL1="$3" URL2="$4" MD5_1="$5" MD5_2="$6"
    local START
    START=$(date +%s)
    local OUT1="${FASTQ_DIR}/${RUN}_1.fastq.gz"
    local OUT2="${FASTQ_DIR}/${RUN}_2.fastq.gz"

    if [[ -f "$OUT1" && -f "$OUT2" ]] && \
       echo "${MD5_1}  ${OUT1}" | md5sum -c --status - 2>/dev/null && \
       echo "${MD5_2}  ${OUT2}" | md5sum -c --status - 2>/dev/null; then
        echo "SKIP ${RUN} (${SAMPLE}): already downloaded and MD5-verified"
        return 0
    fi
    # clear any stale/partial output from an interrupted prior attempt
    rm -f "$OUT1" "$OUT2"

    echo "=== ${RUN} (${SAMPLE}) start $(date) ==="
    local STATUS=0
    curl -sS --retry 5 --retry-delay 10 -o "$OUT1" "https://${URL1}" \
        && curl -sS --retry 5 --retry-delay 10 -o "$OUT2" "https://${URL2}" \
        || STATUS=$?

    # NEGATIVE-LOGIC CHECK: don't trust a curl exit-0 alone -- verify against
    # ENA's own published MD5 before counting this run as done. A corrupted or
    # truncated download that curl still reports success for (e.g. a dropped
    # connection resumed into a garbled file) would otherwise pass silently.
    if [[ $STATUS -eq 0 ]]; then
        if ! echo "${MD5_1}  ${OUT1}" | md5sum -c --status - || \
           ! echo "${MD5_2}  ${OUT2}" | md5sum -c --status -; then
            echo "MD5 MISMATCH for ${RUN} -- removing corrupted output" >&2
            rm -f "$OUT1" "$OUT2"
            STATUS=99
        fi
    fi

    { flock -x 200
      echo -e "${RUN}\t${SAMPLE}\t${STATUS}\t$(( $(date +%s) - START ))s" >> "$JOBLOG"
    } 200>>"${JOBLOG}.lock"
    echo "=== ${RUN} done $(date), exit ${STATUS} ==="
}
export -f download_one

echo -e "run_accession\tsample_title\texit_status\tseconds" > "$JOBLOG"

tail -n +2 "$RUN_LIST" | \
  xargs -P "${SLURM_CPUS_PER_TASK}" -L 1 -I{} bash -c '
    IFS=$'"'"'\t'"'"' read -r run sample exp mix url1 url2 md5_1 md5_2 <<< "{}"
    download_one "$run" "$sample" "$url1" "$url2" "$md5_1" "$md5_2"
  '

# NEGATIVE-LOGIC CHECK: assert the expected sample count, don't just assume
# "the loop finished" means "all 17 samples are present" -- this exact class
# of silent gap (fewer files than expected) already bit the Shi2024 (74/83)
# and Chitarrini2020 (29/33) reprocessing in this project.
N_OK=$(awk -F'\t' 'NR>1 && $3==0' "$JOBLOG" | wc -l)
echo "Runs downloaded and MD5-verified: ${N_OK} / 17"
if [[ "$N_OK" -ne 17 ]]; then
    echo "WARNING: expected all 17 runs to succeed -- check ${JOBLOG} for non-zero exit_status rows before proceeding to alignment." >&2
fi

echo "Done. Fastqs in ${FASTQ_DIR}. Joblog: ${JOBLOG}"
