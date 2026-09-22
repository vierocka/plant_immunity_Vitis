#!/bin/bash -l
#SBATCH --cpus-per-task=8
#SBATCH --mem=16gb
#SBATCH --time=24:00:00
#SBATCH --output=chitarrini2020_download_%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=youremail@here
#SBATCH --account=xxx

###############################################################################
# Standalone download script for Chitarrini et al. 2020 (ENA PRJEB28042, 33
# runs) - NEW (2026-09-07). No such script existed before in this project;
# the raw fastqs backing both chitarrini2020_star_featurecounts_ramses.sh
# (fastp version) and chitarrini2020_trimmomatic_star_fc_ramses.sh
# (Trimmomatic version) were apparently fetched by some other/manual process
# not captured as a script - which is exactly why ERR2987455 went missing
# without anyone noticing until the Trimmomatic run hit it directly (2026-09-07).
#
# VERIFIED DIRECTLY against the ENA portal API (2026-09-07): all 33 runs in
# chitarrini_run_list.tsv are STILL PRESENT on ENA with valid fastq_ftp URLs -
# ERR2987455 specifically confirmed present, NOT withdrawn. So a missing local
# file on RAMSES is a download gap, not a genuine data-loss/withdrawal case -
# do not assume the latter without checking first.
#
# Same MD5-verified, skip-if-already-downloaded, negative-logic-checked
# pattern as froussios2019_download_ramses.sh - built here specifically
# because that pattern was missing for this dataset and directly caused the
# 2026-09-07 mid-run failure.
###############################################################################

set -euo pipefail

WORK_DIR="/scratch/USERNAME/Vitis/Chitarrini2020"
FASTQ_DIR="${WORK_DIR}"   # matches READS_RAW in chitarrini2020_trimmomatic_star_fc_ramses.sh (files land directly in WORK_DIR, not a fastq/ subdir)
RUN_LIST="${WORK_DIR}/chitarrini_run_list_with_md5.tsv"   # scp alongside this script - see build note below
JOBLOG="${WORK_DIR}/download_joblog.tsv"

mkdir -p "$WORK_DIR"
cd "$WORK_DIR"
# BUILD NOTE: chitarrini_run_list.tsv (already in this repo folder) has no MD5
# column. Before submitting, either (a) scp the MD5-augmented version fetched
# locally on 2026-09-07 (ask for it - saved alongside this script when it was
# written), or (b) regenerate directly on RAMSES if internet access allows:
#   curl -sS "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB28042&result=read_run&fields=run_accession,fastq_ftp,fastq_md5&format=tsv&limit=0" -o chitarrini_run_list_with_md5.tsv

export FASTQ_DIR JOBLOG

download_one() {
    local RUN="$1" URL1="$2" URL2="$3" MD5_1="$4" MD5_2="$5"
    local OUT1="${FASTQ_DIR}/${RUN}_1.fastq.gz"
    local OUT2="${FASTQ_DIR}/${RUN}_2.fastq.gz"

    if [[ -f "$OUT1" && -f "$OUT2" ]] && \
       echo "${MD5_1}  ${OUT1}" | md5sum -c --status - 2>/dev/null && \
       echo "${MD5_2}  ${OUT2}" | md5sum -c --status - 2>/dev/null; then
        echo "SKIP ${RUN}: already downloaded and MD5-verified"
        return 0
    fi
    rm -f "$OUT1" "$OUT2"

    local STATUS=0
    curl -sS --retry 5 --retry-delay 10 -o "$OUT1" "https://${URL1}" \
        && curl -sS --retry 5 --retry-delay 10 -o "$OUT2" "https://${URL2}" \
        || STATUS=$?

    if [[ $STATUS -eq 0 ]]; then
        if ! echo "${MD5_1}  ${OUT1}" | md5sum -c --status - || \
           ! echo "${MD5_2}  ${OUT2}" | md5sum -c --status -; then
            echo "MD5 MISMATCH for ${RUN} -- removing corrupted output" >&2
            rm -f "$OUT1" "$OUT2"
            STATUS=99
        fi
    fi
    { flock -x 200; echo -e "${RUN}\t${STATUS}" >> "$JOBLOG"; } 200>>"${JOBLOG}.lock"
}
export -f download_one

echo -e "run_accession\texit_status" > "$JOBLOG"

# plain bash job control, not xargs -I{} - avoids a field-mangling bug:
# bash's plain echo does not expand \t, and nested quoting through xargs
# corrupted fields further. Use parameter expansion + & + wait -n instead.
MAX_JOBS=8
tail -n +2 "$RUN_LIST" | \
  while IFS=$'\t' read -r run ftp md5; do
    url1="${ftp%%;*}"; url2="${ftp#*;}"
    md5_1="${md5%%;*}"; md5_2="${md5#*;}"
    download_one "$run" "$url1" "$url2" "$md5_1" "$md5_2" &
    while (( $(jobs -rp | wc -l) >= MAX_JOBS )); do wait -n || true; done
  done
wait

N_OK=$(awk -F'\t' 'NR>1 && $2==0' "$JOBLOG" | wc -l)
N_EXPECTED=$(( $(wc -l < "$RUN_LIST") - 1 ))
echo "Runs downloaded and MD5-verified: ${N_OK} / ${N_EXPECTED}"
if [[ "$N_OK" -ne "$N_EXPECTED" ]]; then
  echo "WARNING: expected all ${N_EXPECTED} runs to succeed -- check ${JOBLOG} for non-zero exit_status rows before proceeding to trimming/alignment." >&2
fi

echo "Done. Fastqs in ${FASTQ_DIR}. Joblog: ${JOBLOG}"
