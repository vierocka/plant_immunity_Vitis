#!/bin/bash -l
#SBATCH --cpus-per-task=24
#SBATCH --mem=192gb
#SBATCH --time=72:00:00
#SBATCH --output=froussios2019_trimmomatic_STAR_featureCounts_%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=youremail@here
#SBATCH --account=xxx

###############################################################################
# Reprocessing of Froussios et al. 2019 (Bioinformatics 35(18):3372-3377,
# ArrayExpress E-MTAB-5446 / ENA ERP021226, 17 runs, WT Arabidopsis Col-0,
# HiSeq2000 2x101bp) with THIS PROJECT'S OWN main pipeline -- Trimmomatic ->
# STAR -> featureCounts, same tools/parameters as the local scripts
# (~/Dropbox/MendelUni_Vinselect/scripts/trimmomatic.sh,
# star_alignm_Vvinif_Pviticola_r1.sh, star1_onlyF_unpaired.sh,
# featureCounts.sh) -- NOT the fastp-Q30-only pipeline used previously in
# froussios2019_star_featurecounts_ramses.sh.
#
# WHY THIS REPLACES THE fastp VERSION: fastp's Q30-only paired mode DROPS THE
# ENTIRE PAIR if either mate fails, unlike Trimmomatic's native _1P/_2P
# (both-survived) vs _1U/_2U (singleton-survived) separation -- so the fastp
# version cannot produce the forward-only-survived-read-fraction diagnostic
# used throughout real_data_Fonly_degradation_DE_investigation.md (this was
# already flagged as a blocking gap there: "Froussios_2019 checked... NOT
# directly usable for the F/R-only breakdown as-is"). This script is that
# flagged TODO, done via a full Trimmomatic rerun rather than patching fastp.
#
# CAVEAT carried over from froussios2019_download_ramses.sh: these 17 runs
# are NOT one flat single-batch replicate set -- they split into 3 sequencing
# batches ("Experiment" 1/2/3, n=7/7/3, see experiment_batch column in
# froussios2019_run_list.tsv) further crossed with 2 ERCC spike-in mixes.
# Useful in its own right (a real, independently-documented multi-batch case
# to test the degradation/batch framework against) but NOT a clean single-
# batch dataset if that specific property is wanted for the "more datasets"
# search -- keep the batch column when analyzing results here.
#
# Deliberately mirrors the local Vitis scripts' STAR filter set (--outFilter
# MultimapNmax 10 etc.), a deliberate departure from froussios2019_star_
# featurecounts_ramses.sh's --outFilterMultimapNmax 1 uniquely-mapped-only
# convention, for parity with the main pipeline.
###############################################################################

set -euo pipefail

# FIX (2026-09-04): STAR's SortedByCoordinate BAM output opens many temp files
# in parallel across sort bins x threads (--outBAMsortingBinsN default 50 x
# --runThreadN 24 here) - the identical Chitarrini2020 script hit "could not
# create output file .../BAMsort/18/42, FATAL ERROR" on RAMSES from this;
# STAR's own message points at ulimit -n (not RAM despite how the error
# looks). Raise the open-file limit and cap the bin count as a belt-and-
# suspenders fix, applied here pre-emptively since this script shares the
# exact same STAR invocation pattern. If a run already failed with this
# error, also manually remove its stale "<sample>__STARtmp" directory before
# resubmitting - STAR does not clean that up on a failed run.
ulimit -n 65536 2>/dev/null || ulimit -n 10000 2>/dev/null || true

REF_DIR="/scratch/USERNAME/Athal/TAIR10ref"
WORK_DIR="/scratch/USERNAME/Athal"
FASTQ_DIR="/scratch/USERNAME/Athal/fastq"     # from froussios2019_download_ramses.sh
RUN_LIST="/scratch/USERNAME/Athal/froussios2019_run_list.tsv"
GENOME_FA="${REF_DIR}/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
ANNOT_GFF3="${REF_DIR}/Arabidopsis_thaliana.TAIR10.63.gff3"
STAR_INDEX="${REF_DIR}/STAR_index_TAIR10"
ADAPTERS="${WORK_DIR}/TruSeq3-PE.fa"   # scp from ~/anaconda3/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa
                                        # before submitting -- OR check `module spider Trimmomatic`
                                        # for a RAMSES-bundled adapters dir and point here instead.

TRIM_DIR="${WORK_DIR}/trimmomatic"
STAR_PAIRED_DIR="${WORK_DIR}/STARmap1"
STAR_FONLY_DIR="${WORK_DIR}/STARmap1_F"

mkdir -p "$WORK_DIR" "$TRIM_DIR" "$STAR_PAIRED_DIR" "$STAR_FONLY_DIR" "${STAR_PAIRED_DIR}/fC" "${STAR_FONLY_DIR}/fC"
cd "$WORK_DIR"

############################ MODULES ############################################
# CHECK EXACT MODULE NAMES/VERSIONS before first submit: `module spider Trimmomatic`.
load_trimmomatic() { module unload compiler/GCC 2>/dev/null || true; module add bio/Trimmomatic/0.39-Java-17; }
load_star()        { module unload compiler/GCC 2>/dev/null || true; module add bio/STAR/2.7.10b-GCC-11.3.0; }
load_subread()      { module unload compiler/GCC 2>/dev/null || true; module add bio/Subread/2.1.1-GCC-13.2.0; }

############################ 0. STAR GENOME INDEX (reuse if already built by the fastp-version script) ##
load_star
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

############################ 1. PER-SAMPLE: Trimmomatic PE -> STAR (paired + F-only) ####
tail -n +2 "$RUN_LIST" | while IFS=$'\t' read -r run sample exp mix url1 url2 md5_1 md5_2; do
  R1="${FASTQ_DIR}/${run}_1.fastq.gz"
  R2="${FASTQ_DIR}/${run}_2.fastq.gz"
  # FIX (2026-09-07, same bug/fix as chitarrini2020_trimmomatic_star_fc_ramses.sh):
  # a missing raw fastq used to `exit 1` here, aborting the WHOLE per-sample loop -
  # meaning featureCounts (only run after this entire loop finishes) never ran even
  # for samples that had already completed successfully. Confirmed against ENA
  # (2026-09-07): ERR1811903 is NOT withdrawn, still has a valid fastq_ftp URL -
  # this is a local download gap, not a data-availability problem, but user is
  # OK proceeding without it ("I am ok with ExpA and ExpB"). Now: log and skip
  # just this one run, let the rest of the batch (and featureCounts) complete.
  if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
    echo "WARNING: missing raw fastq for ${run} -- run froussios2019_download_ramses.sh first (or check for a transient download gap). Skipping this run, continuing with the rest." >&2
    continue
  fi

  if [ ! -f "${TRIM_DIR}/${run}_1P.fq.gz" ]; then
    load_trimmomatic
    trimmomatic PE -threads 24 -trimlog "${TRIM_DIR}/${run}.log" -summary "${TRIM_DIR}/${run}_summary.txt" \
      -quiet -validatePairs \
      "$R1" "$R2" \
      "${TRIM_DIR}/${run}_1P.fq.gz" "${TRIM_DIR}/${run}_1U.fq.gz" \
      "${TRIM_DIR}/${run}_2P.fq.gz" "${TRIM_DIR}/${run}_2U.fq.gz" \
      ILLUMINACLIP:${ADAPTERS}:2:28:8 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
  fi

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
done

############################ 2. featureCounts -- paired and F-only, separately ####
load_subread

tail -n +2 "$RUN_LIST" | while IFS=$'\t' read -r run sample exp mix url1 url2 md5_1 md5_2; do
  # FIX (2026-09-07): a run skipped earlier for missing raw fastq has no BAM here -
  # featureCounts on a nonexistent file would fail and, under set -e, abort the
  # WHOLE remaining loop (the exact same failure class as the exit-1 bug fixed
  # above, just one step later). Skip with a warning instead.
  if [ ! -f "${STAR_PAIRED_DIR}/${run}_Aligned.sortedByCoord.out.bam" ]; then
    echo "WARNING: no paired BAM for ${run} (likely skipped earlier for missing fastq) -- skipping featureCounts for this run." >&2
  elif [ ! -f "${STAR_PAIRED_DIR}/fC/${run}.count" ]; then
    featureCounts -T 24 -p -a "$ANNOT_GFF3" -t gene -g gene_id \
      -o "${STAR_PAIRED_DIR}/fC/${run}.count" \
      "${STAR_PAIRED_DIR}/${run}_Aligned.sortedByCoord.out.bam"
  fi
  if [ ! -f "${STAR_FONLY_DIR}/${run}_Aligned.sortedByCoord.out.bam" ]; then
    echo "WARNING: no F-only BAM for ${run} (likely skipped earlier for missing fastq) -- skipping featureCounts for this run." >&2
  elif [ ! -f "${STAR_FONLY_DIR}/fC/${run}.count" ]; then
    featureCounts -T 24 -a "$ANNOT_GFF3" -t gene -g gene_id \
      -o "${STAR_FONLY_DIR}/fC/${run}.count" \
      "${STAR_FONLY_DIR}/${run}_Aligned.sortedByCoord.out.bam"
  fi
done

N_EXPECTED=$(( $(wc -l < "$RUN_LIST") - 1 ))
N_PAIRED=$(ls "${STAR_PAIRED_DIR}/fC/"*.count 2>/dev/null | wc -l)
N_FONLY=$(ls "${STAR_FONLY_DIR}/fC/"*.count 2>/dev/null | wc -l)
echo "Expected samples: ${N_EXPECTED}. Paired count files: ${N_PAIRED}. F-only count files: ${N_FONLY}."
if [ "$N_PAIRED" -ne "$N_EXPECTED" ] || [ "$N_FONLY" -ne "$N_EXPECTED" ]; then
  echo "WARNING: mismatch -- check which run(s) failed before trusting these results." >&2
fi

echo "Done. Paired counts: ${STAR_PAIRED_DIR}/fC/. F-only counts: ${STAR_FONLY_DIR}/fC/."
