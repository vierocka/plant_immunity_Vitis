#!/bin/bash -l
#SBATCH --cpus-per-task=24
#SBATCH --mem=128gb
#SBATCH --time=24:00:00
#SBATCH --output=froussios2019_add_samples15_16_17_%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=youremail@here
#SBATCH --account=xxx

###############################################################################
# Fill the Experiment-3 gap: samples 15/16/17 (ERR1811902/ERR1811903/
# ERR1811904) were never reliably carried through this project's own
# Trimmomatic->STAR->featureCounts pipeline the way samples 1-14 (ExpA/ExpB)
# were -- ERR1811903 in particular was a known download gap noted in
# froussios2019_trimmomatic_star_fc_ramses.sh ("user is OK with ExpA and
# ExpB" at the time). This script closes that gap end-to-end: download via
# sra-tools (prefetch+fasterq-dump) -> Trimmomatic PE -> STAR (paired/
# STARmap1, F-only/STARmap1_F, R-only/STARmap1_R) -> featureCounts, writing
# into the SAME shared directories the main 14-sample pipeline uses, so the
# result is indistinguishable from a sample that went through the main
# script -- not a parallel/separate output tree.
#
# CHECKED BEFORE WRITING THIS (2026-09-14), not assumed:
#   - ERR1811902 (sample 15): local Dropbox has PARTIAL downstream outputs
#     only -- STAR_logs/*/ERR1811902_Log.final.out, trimm_reports/
#     ERR1811902_summary.txt, and fC/STAR_{both,Fonly,Ronly}/ERR1811902.count
#     (dated 2026-09-08, non-empty). User confirmed (2026-09-14) the
#     upstream data -- raw fastq, trimmed fastq, BAMs -- was DELETED from the
#     cluster (RAMSES scratch). So the local .count files are the only
#     surviving artifact; they are NOT reproducible from current scratch
#     state and every stage below (fastq/trim/BAM/count, all checked against
#     RAMSES scratch paths, not local Dropbox) WILL redo sample 15 from
#     scratch -- correct and necessary here, not wasted work.
#   - Do NOT let this script's step 4 (featureCounts) silently overwrite the
#     existing local fC/STAR_*/ERR1811902.count files when synced back --
#     those are the only current record of the original run. Diff the new
#     RAMSES output against them first (see the note at the end of this
#     script) before replacing, per this project's extend-don't-overwrite
#     convention.
#   - ERR1811903/ERR1811904 (samples 16/17): confirmed genuinely absent --
#     no fastq, no trimmomatic output, no BAM, no count file anywhere in
#     this project's local tree, and (per the original script's own note)
#     no successful download on the cluster either.
#
# WHY sra-tools HERE, NOT froussios2019_download_ramses.sh's direct-ENA-FTP
# method: that script's whole rationale (see its own header) was avoiding
# the SRA round-trip because direct ENA fastq + published MD5 gives free
# integrity verification. That reasoning still holds in general -- but if
# the ENA-FTP path is what silently produced the sample-16/17 gap in the
# first place (transient host issue, path typo, etc.), retrying the exact
# same method risks repeating the same silent failure. sra-tools resolves
# ERR accessions through the SRA/ENA cross-archive mirror independently of
# the ENA FTP tree, so it's a genuinely different code path, not just a
# retry -- reasonable to reach for here even though it's the non-default
# choice project-wide.
#
# CAVEAT: fasterq-dump re-derives fastq from the .sra object -- the result
# is NOT byte-identical to ENA's own hosted fastq.gz, so the published
# md5_R1/md5_R2 in froussios2019_run_list.tsv CANNOT be used to verify these
# downloads (this is a real re-encoding, not a corruption risk to check the
# same way). Verified instead by: (a) R1/R2 read-count equality (fasterq-dump
# --split-files on a genuinely paired run must produce matched mate counts --
# a mismatch means something is wrong with the accession's layout or the
# extraction), (b) a sanity floor against the ~77-130M read-pairs/sample
# range independently verified for this whole dataset in NOTES.md.
#
# CHECK BEFORE FIRST SUBMIT: exact SRA-Toolkit module name/version on RAMSES
# (`module spider SRA-Toolkit` or `module spider sra-tools`) -- set
# SRA_MODULE below. Same caveat as every other RAMSES script in this
# project: module names are RAMSES-specific and unverified from where this
# was written.
###############################################################################

set -euo pipefail
ulimit -n 65536 2>/dev/null || ulimit -n 10000 2>/dev/null || true

REF_DIR="/scratch/USERNAME/Athal/TAIR10ref"
WORK_DIR="/scratch/USERNAME/Athal"
FASTQ_DIR="${WORK_DIR}/fastq"                 # same dir the ENA-download script uses -- output lands here indistinguishably
SRA_TMP_DIR="${WORK_DIR}/sra_tmp"             # scratch space for fasterq-dump's intermediate uncompressed fastq + prefetch cache
GENOME_FA="${REF_DIR}/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
ANNOT_GFF3="${REF_DIR}/Arabidopsis_thaliana.TAIR10.63.gff3"
STAR_INDEX="${REF_DIR}/STAR_index_TAIR10"      # this project's OWN index (sjdbOverhang 100) -- NOT the paperrepro scripts' separate sjdbOverhang99 index

TRIM_DIR="${WORK_DIR}/trimmomatic"
STAR_PAIRED_DIR="${WORK_DIR}/STARmap1"
STAR_FONLY_DIR="${WORK_DIR}/STARmap1_F"
STAR_RONLY_DIR="${WORK_DIR}/STARmap1_R"

mkdir -p "$FASTQ_DIR" "$SRA_TMP_DIR" "$TRIM_DIR" \
         "${STAR_PAIRED_DIR}/fC" "${STAR_FONLY_DIR}/fC" "${STAR_RONLY_DIR}/fC"
cd "$WORK_DIR"

[ -f "$GENOME_FA" ] || { echo "ERROR: ${GENOME_FA} not found -- this script assumes the reference already exists from the main pipeline." >&2; exit 1; }
[ -f "$ANNOT_GFF3" ] || { echo "ERROR: ${ANNOT_GFF3} not found -- this script assumes the reference already exists from the main pipeline." >&2; exit 1; }

############################ THE THREE SAMPLES (from froussios2019_run_list.tsv / NOTES.md) ####
# run_accession  sample_title  experiment_batch  ercc_mix
RUNS=(
  "ERR1811902 Sample_15 3 1"
  "ERR1811903 Sample_16 3 2"
  "ERR1811904 Sample_17 3 1"
)

############################ MODULES ############################################
# UPDATE to the exact name/version from `module spider SRA-Toolkit` before first submit
load_sra()         { module unload compiler/GCC 2>/dev/null || true; module add bio/SRA-Toolkit/3.0.10-gompi-2023a; }
load_trimmomatic() { module unload compiler/GCC 2>/dev/null || true; module add module add lang/Java/17.0.6; module add bio/Trimmomatic/0.39-Java-17; }
load_star()        { module unload compiler/GCC 2>/dev/null || true; module add bio/STAR/2.7.10b-GCC-11.3.0; }
load_subread()     { module unload compiler/GCC 2>/dev/null || true; module add bio/Subread/2.1.1-GCC-13.2.0; }

############################ 0. STAR GENOME INDEX (reuse if already built by the main pipeline) ####
load_star
if [ ! -f "${STAR_INDEX}/SAindex" ]; then
  mkdir -p "$STAR_INDEX"
  STAR --runThreadN 24 --runMode genomeGenerate \
       --genomeDir "$STAR_INDEX" \
       --genomeFastaFiles "$GENOME_FA" \
       --sjdbGTFfile "$ANNOT_GFF3" \
       --sjdbGTFtagExonParentTranscript Parent \
       --sjdbGTFfeatureExon exon \
       --sjdbOverhang 100 \
       --genomeSAindexNbases 12
fi

############################ 1. DOWNLOAD via sra-tools (skip-if-already-present) ####
for row in "${RUNS[@]}"; do
  read -r run sample exp mix <<< "$row"
  OUT1="${FASTQ_DIR}/${run}_1.fastq.gz"
  OUT2="${FASTQ_DIR}/${run}_2.fastq.gz"

  if [ -f "$OUT1" ] && [ -f "$OUT2" ]; then
    echo "SKIP download ${run} (${sample}): fastq already present at ${FASTQ_DIR}"
    continue
  fi

  echo "=== ${run} (${sample}) download start $(date) ==="
  load_sra
  command -v prefetch >/dev/null 2>&1 || { echo "ERROR: prefetch not found after loading ${SRA_MODULE} -- check module name (module spider SRA-Toolkit / sra-tools)." >&2; exit 1; }

  # prefetch refuses to run while a lock file exists, even one left behind by
  # a killed/interrupted prior attempt (a real failure mode here: this
  # exact script has already been aborted once mid-run by an unrelated
  # module-loading bug). A stale lock with no genuinely complete .sra next
  # to it just wedges every future resubmit until removed by hand -- so
  # clear it here automatically rather than requiring a manual `rm -rf` on
  # the cluster each time this happens.
  if [ -f "${SRA_TMP_DIR}/${run}/${run}.sra.lock" ]; then
    echo "NOTE: stale lock found for ${run} (leftover from an interrupted prior attempt) -- removing ${SRA_TMP_DIR}/${run} before retrying." >&2
    rm -rf "${SRA_TMP_DIR}/${run}"
  fi

  prefetch --max-size 100G --output-directory "$SRA_TMP_DIR" "$run"

  fasterq-dump --split-files --threads "${SLURM_CPUS_PER_TASK:-24}" \
    --temp "$SRA_TMP_DIR" --outdir "$SRA_TMP_DIR" \
    "${SRA_TMP_DIR}/${run}/${run}.sra"

  [ -f "${SRA_TMP_DIR}/${run}_1.fastq" ] && [ -f "${SRA_TMP_DIR}/${run}_2.fastq" ] || \
    { echo "ERROR: fasterq-dump did not produce paired output for ${run} -- check whether this accession is genuinely paired-end (expected, per NOTES.md)." >&2; exit 1; }

  # NEGATIVE-LOGIC CHECK: mate read-count equality and a sanity floor against
  # the independently-verified 77-130M read-pairs/sample range (NOTES.md) --
  # don't just trust a zero exit status from fasterq-dump.
  N1=$(( $(wc -l < "${SRA_TMP_DIR}/${run}_1.fastq") / 4 ))
  N2=$(( $(wc -l < "${SRA_TMP_DIR}/${run}_2.fastq") / 4 ))
  if [ "$N1" -ne "$N2" ]; then
    echo "ERROR: ${run} mate read counts differ (R1=${N1}, R2=${N2}) -- pairing is broken, do not proceed with this sample." >&2
    exit 1
  fi
  if [ "$N1" -lt 50000000 ]; then
    echo "WARNING: ${run} has only ${N1} read pairs -- well below the ~77-130M/sample range verified for this dataset in NOTES.md. Proceeding, but flag this before trusting downstream results." >&2
  fi
  echo "${run}: ${N1} read pairs extracted, mate counts match."

  gzip -c "${SRA_TMP_DIR}/${run}_1.fastq" > "$OUT1"
  gzip -c "${SRA_TMP_DIR}/${run}_2.fastq" > "$OUT2"

  # free scratch space -- these intermediates are fully reproducible from the
  # gzipped output above and the .sra cache is large (tens of GB across 3 runs)
  rm -rf "${SRA_TMP_DIR}/${run}" "${SRA_TMP_DIR}/${run}_1.fastq" "${SRA_TMP_DIR}/${run}_2.fastq"
  echo "=== ${run} download done $(date) ==="
done

############################ 2. Trimmomatic PE (identical params to froussios2019_trimmomatic_star_fc_ramses.sh) ####
for row in "${RUNS[@]}"; do
  read -r run sample exp mix <<< "$row"
  R1="${FASTQ_DIR}/${run}_1.fastq.gz"
  R2="${FASTQ_DIR}/${run}_2.fastq.gz"
  [ -f "$R1" ] && [ -f "$R2" ] || { echo "WARNING: missing fastq for ${run} -- download step above must have failed. Skipping." >&2; continue; }

  if [ ! -f "${TRIM_DIR}/${run}_1P.fq.gz" ]; then
    load_trimmomatic
    ADAPTERS="${EBROOTTRIMMOMATIC}/adapters/TruSeq3-PE.fa"
    java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads 24 -trimlog "${TRIM_DIR}/${run}.log" -summary "${TRIM_DIR}/${run}_summary.txt" \
      -quiet -validatePairs \
      "$R1" "$R2" \
      "${TRIM_DIR}/${run}_1P.fq.gz" "${TRIM_DIR}/${run}_1U.fq.gz" \
      "${TRIM_DIR}/${run}_2P.fq.gz" "${TRIM_DIR}/${run}_2U.fq.gz" \
      ILLUMINACLIP:${ADAPTERS}:2:28:8 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
  fi
done

############################ 3. STAR -- paired, F-only, R-only (identical params to the main pipeline / chitarrini2020_ronly_star_fc_ramses.sh) ####
for row in "${RUNS[@]}"; do
  read -r run sample exp mix <<< "$row"
  [ -f "${TRIM_DIR}/${run}_1P.fq.gz" ] || { echo "WARNING: no trimmed reads for ${run} -- skipping STAR." >&2; continue; }

  if [ ! -f "${STAR_PAIRED_DIR}/${run}_Aligned.sortedByCoord.out.bam" ]; then
    load_star
    STAR --runThreadN 24 --genomeDir "$STAR_INDEX" \
      --outFilterType BySJout --outFilterMultimapNmax 10 \
      --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
      --outFilterMismatchNmax 12 --outFilterMismatchNoverReadLmax 0.04 \
      --alignIntronMin 20 --outFileNamePrefix "${STAR_PAIRED_DIR}/${run}_" \
      --outReadsUnmapped Fastx --limitBAMsortRAM 12000000000 \
      --outBAMsortingBinsN 20 \
      --outSAMtype BAM SortedByCoordinate --readFilesCommand gunzip -c \
      --readFilesIn "${TRIM_DIR}/${run}_1P.fq.gz" "${TRIM_DIR}/${run}_2P.fq.gz"
  fi

  if [ ! -f "${STAR_FONLY_DIR}/${run}_Aligned.sortedByCoord.out.bam" ]; then
    load_star
    STAR --runThreadN 24 --genomeDir "$STAR_INDEX" \
      --outFilterType BySJout --outFilterMultimapNmax 10 \
      --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
      --outFilterMismatchNmax 12 --outFilterMismatchNoverReadLmax 0.04 \
      --alignIntronMin 20 --outFileNamePrefix "${STAR_FONLY_DIR}/${run}_" \
      --outReadsUnmapped Fastx --limitBAMsortRAM 12000000000 \
      --outBAMsortingBinsN 20 \
      --outSAMtype BAM SortedByCoordinate --readFilesCommand gunzip -c \
      --readFilesIn "${TRIM_DIR}/${run}_1U.fq.gz"
  fi

  if [ ! -f "${STAR_RONLY_DIR}/${run}_Aligned.sortedByCoord.out.bam" ]; then
    load_star
    STAR --runThreadN 24 --genomeDir "$STAR_INDEX" \
      --outFilterType BySJout --outFilterMultimapNmax 10 \
      --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
      --outFilterMismatchNmax 12 --outFilterMismatchNoverReadLmax 0.04 \
      --alignIntronMin 20 --outFileNamePrefix "${STAR_RONLY_DIR}/${run}_" \
      --outReadsUnmapped Fastx --limitBAMsortRAM 12000000000 \
      --outBAMsortingBinsN 20 \
      --outSAMtype BAM SortedByCoordinate --readFilesCommand gunzip -c \
      --readFilesIn "${TRIM_DIR}/${run}_2U.fq.gz"
  fi
done

############################ 4. featureCounts -- paired, F-only, R-only ####
load_subread
for row in "${RUNS[@]}"; do
  read -r run sample exp mix <<< "$row"

  if [ ! -f "${STAR_PAIRED_DIR}/${run}_Aligned.sortedByCoord.out.bam" ]; then
    echo "WARNING: no paired BAM for ${run} -- skipping featureCounts (both)." >&2
  elif [ ! -f "${STAR_PAIRED_DIR}/fC/${run}.count" ]; then
    featureCounts -T 24 -p -a "$ANNOT_GFF3" -t gene -g gene_id \
      -o "${STAR_PAIRED_DIR}/fC/${run}.count" \
      "${STAR_PAIRED_DIR}/${run}_Aligned.sortedByCoord.out.bam"
  fi

  if [ ! -f "${STAR_FONLY_DIR}/${run}_Aligned.sortedByCoord.out.bam" ]; then
    echo "WARNING: no F-only BAM for ${run} -- skipping featureCounts (F-only)." >&2
  elif [ ! -f "${STAR_FONLY_DIR}/fC/${run}.count" ]; then
    featureCounts -T 24 -a "$ANNOT_GFF3" -t gene -g gene_id \
      -o "${STAR_FONLY_DIR}/fC/${run}.count" \
      "${STAR_FONLY_DIR}/${run}_Aligned.sortedByCoord.out.bam"
  fi

  if [ ! -f "${STAR_RONLY_DIR}/${run}_Aligned.sortedByCoord.out.bam" ]; then
    echo "WARNING: no R-only BAM for ${run} -- skipping featureCounts (R-only)." >&2
  elif [ ! -f "${STAR_RONLY_DIR}/fC/${run}.count" ]; then
    featureCounts -T 24 -a "$ANNOT_GFF3" -t gene -g gene_id \
      -o "${STAR_RONLY_DIR}/fC/${run}.count" \
      "${STAR_RONLY_DIR}/${run}_Aligned.sortedByCoord.out.bam"
  fi
done

############################ 5. NEGATIVE-LOGIC CHECK: assert all 3 runs x 3 count sets exist, don't assume ####
N_EXPECTED=${#RUNS[@]}
N_BOTH=0; N_F=0; N_R=0
for row in "${RUNS[@]}"; do
  read -r run sample exp mix <<< "$row"
  [ -s "${STAR_PAIRED_DIR}/fC/${run}.count" ] && N_BOTH=$((N_BOTH+1))
  [ -s "${STAR_FONLY_DIR}/fC/${run}.count" ] && N_F=$((N_F+1))
  [ -s "${STAR_RONLY_DIR}/fC/${run}.count" ] && N_R=$((N_R+1))
done
echo "Expected samples: ${N_EXPECTED}. both=${N_BOTH} Fonly=${N_F} Ronly=${N_R}."
if [ "$N_BOTH" -ne "$N_EXPECTED" ] || [ "$N_F" -ne "$N_EXPECTED" ] || [ "$N_R" -ne "$N_EXPECTED" ]; then
  echo "WARNING: mismatch -- check which run(s)/stage(s) failed above before trusting this batch as complete." >&2
fi

echo "Done. Fastq: ${FASTQ_DIR}/. Trimmed: ${TRIM_DIR}/. Counts: ${STAR_PAIRED_DIR}/fC/, ${STAR_FONLY_DIR}/fC/, ${STAR_RONLY_DIR}/fC/."
echo "Next: scp the new .count/.count.summary files back to a STAGING location first (e.g. fC_rerun_sample15_16_17/), NOT directly over"
echo "the existing local fC/STAR_{both,Fonly,Ronly}/ERR1811902.count -- diff read totals against the current (2026-09-08) local copy before"
echo "replacing it, since that file is currently the only surviving record of the original (pre-deletion) sample-15 run."
