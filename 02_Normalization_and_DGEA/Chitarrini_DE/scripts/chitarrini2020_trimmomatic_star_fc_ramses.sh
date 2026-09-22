#!/bin/bash -l
#SBATCH --cpus-per-task=24
#SBATCH --mem=192gb
#SBATCH --time=72:00:00
#SBATCH --output=chitarrini2020_trimmomatic_STAR_featureCounts_%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=youremail@here
#SBATCH --account=xxx

###############################################################################
# Reprocessing of Chitarrini et al. 2020 (Sci Rep 10:12193, ENA PRJEB28042, 33
# runs: 0/12/24/48/96/120 hpi x mock/inoculated x 3 reps, HiSeq2500 2x50bp)
# with THIS PROJECT'S OWN main pipeline -- Trimmomatic -> STAR -> featureCounts,
# same tools/parameters as the local scripts (~/Dropbox/MendelUni_Vinselect/
# scripts/trimmomatic.sh, star_alignm_Vvinif_Pviticola_r1.sh, star1_onlyF_
# unpaired.sh, featureCounts.sh) -- NOT the fastp-Q30-only pipeline used in
# chitarrini2020_star_featurecounts_ramses.sh.
#
# WHY THIS REPLACES THE fastp VERSION: the fastp pipeline (Q30 filter only,
# no adapter/length trimming) does not separate reads into paired-survived
# vs forward-only-survived output the way Trimmomatic's PE mode natively
# does (_1P/_2P vs _1U/_2U) -- so it cannot produce the F-only-survived-read-
# fraction diagnostic used throughout real_data_Fonly_degradation_DE_
# investigation.md (Findings 1/4/10/13/14 in that memory file). Trimmomatic's
# _1U output IS that quantity directly, with zero extra processing.
#
# Two independent goals this reprocessing serves:
#  1. Genuine methodological parity with the main Vitis analysis, so any
#     F-only-vs-DE association found here is comparable (same trimming/
#     mapping/counting rules), not an artifact of a different pipeline.
#  2. Chitarrini2020 is a candidate independent-validation dataset for the
#     "does the degradation-vs-spurious-DE pattern reproduce in an unrelated,
#     nobody-previously-checked real dataset" question (see the memory file's
#     "Independent validation idea" TODO) -- SINGLE lab, SINGLE library prep
#     batch as far as the paper reports (verify this assumption against the
#     paper/ENA metadata before treating it as a clean single-batch case).
#
# Deliberately mirrors the local Vitis scripts' STAR filter set exactly
# (--outFilterMultimapNmax 10, not the fastp-script's --outFilterMultimapNmax
# 1 uniquely-mapped-only convention) -- this is a DELIBERATE departure from
# chitarrini2020_star_featurecounts_ramses.sh, made for parity with the main
# pipeline, not an oversight.
#
# Reference: reuses the same PN40024 v4 STAR index built for the fastp-based
# rerun (REF_DIR below) -- if that index does not already exist, build it
# first (block 1, currently commented out, sjdbOverhang 49 = 50bp reads - 1,
# matches this dataset's 2x50bp reads).
###############################################################################

set -euo pipefail

# FIX (2026-09-04): STAR's SortedByCoordinate BAM output opens many temp files
# in parallel across sort bins x threads (--outBAMsortingBinsN default 50 x
# --runThreadN 24 here) - hit "could not create output file .../BAMsort/18/42,
# FATAL ERROR" on RAMSES, STAR's own message points at ulimit -n (not RAM
# despite how the error looks). Raise the open-file limit and cap the bin
# count as a belt-and-suspenders fix. If a run already failed with this error,
# also manually remove its stale "<sample>__STARtmp" directory before
# resubmitting - STAR does not clean that up on a failed run, and a leftover
# one can cause a new immediate error even after this fix.
ulimit -n 65536 2>/dev/null || ulimit -n 10000 2>/dev/null || true

REF_DIR="/scratch/USERNAME/Vitis/GGref"
WORK_DIR="/scratch/USERNAME/Vitis/Chitarrini2020"
READS_RAW="/scratch/USERNAME/Vitis/Chitarrini2020"   # raw fastqs from the earlier download step
RUN_LIST="/scratch/USERNAME/Vitis/Chitarrini2020/chitarrini_run_list.tsv"
GENOME_FA="${REF_DIR}/Vitis_vinifera.PN40024.v4.dna.toplevel.fa"
ANNOT_GFF3="${REF_DIR}/Vitis_vinifera.PN40024.v4.56.gff3"
STAR_INDEX="${REF_DIR}/STAR_index_PN40024v4"
ADAPTERS="${WORK_DIR}/TruSeq3-PE.fa"   # scp from ~/anaconda3/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa
                                        # before submitting -- OR check `module spider Trimmomatic`
                                        # for a RAMSES-bundled adapters dir and point here instead.

TRIM_DIR="${WORK_DIR}/trimmomatic"
STAR_PAIRED_DIR="${WORK_DIR}/STARmap1"
STAR_FONLY_DIR="${WORK_DIR}/STARmap1_F"
STAR_RONLY_DIR="${WORK_DIR}/STARmap1_R"
FC_DIR="${WORK_DIR}/fC"

mkdir -p "$WORK_DIR" "$TRIM_DIR" "$STAR_PAIRED_DIR" "$STAR_FONLY_DIR" "$STAR_RONLY_DIR" "$FC_DIR" \
  "${STAR_PAIRED_DIR}/fC" "${STAR_FONLY_DIR}/fC" "${STAR_RONLY_DIR}/fC"
cd "$WORK_DIR"

############################ MODULES ############################################
# CHECK EXACT MODULE NAMES/VERSIONS before first submit: `module spider Trimmomatic`.
load_trimmomatic() { module unload compiler/GCC 2>/dev/null || true; module add module add lang/Java/17.0.6; module add bio/Trimmomatic/0.39-Java-13; }
load_star()        { module unload compiler/GCC 2>/dev/null || true; module add bio/STAR/2.7.10b-GCC-11.3.0; }
load_subread()      { module unload compiler/GCC 2>/dev/null || true; module add bio/Subread/2.1.1-GCC-13.2.0; }

############################ 0. STAR GENOME INDEX (build once, reuse if present) ##
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
       --sjdbOverhang 49 \
       --genomeSAindexNbases 13
fi

############################ 1. PER-SAMPLE: Trimmomatic PE -> STAR (paired + F-only) ####
tail -n +2 "$RUN_LIST" | while IFS=$'\t' read -r run sample url1 url2; do
  R1="${READS_RAW}/${run}_1.fastq.gz"
  R2="${READS_RAW}/${run}_2.fastq.gz"
  # FIX (2026-09-07): a missing raw fastq used to `exit 1` here, aborting the
  # WHOLE per-sample loop for every remaining run - meaning featureCounts (which
  # only runs after this entire loop finishes) never ran even for the samples
  # that had already completed trimming+mapping successfully. This actually
  # happened: ERR2987455 was missing locally (confirmed 2026-09-07 NOT
  # withdrawn from ENA - a local download gap, see chitarrini2020_download_
  # ramses.sh, not a data-loss case) and stopped the run partway through with
  # featureCounts never reached. Now: log and skip just this one run, let the
  # rest of the batch complete - the completeness check at the end of this
  # script already reports any gap (expected vs actual count file numbers).
  if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
    echo "WARNING: missing raw fastq for ${run} -- run chitarrini2020_download_ramses.sh first (or check for a transient download gap). Skipping this run, continuing with the rest." >&2
    continue
  fi

  # skip-if-done: resumable per sample (feedback_resumability_default convention)
  if [ ! -f "${TRIM_DIR}/${run}_1P.fq.gz" ]; then
    load_trimmomatic
    java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads 24 -trimlog "${TRIM_DIR}/${run}.log" -summary "${TRIM_DIR}/${run}_summary.txt" \
      -quiet -validatePairs \
      "$R1" "$R2" \
      "${TRIM_DIR}/${run}_1P.fq.gz" "${TRIM_DIR}/${run}_1U.fq.gz" \
      "${TRIM_DIR}/${run}_2P.fq.gz" "${TRIM_DIR}/${run}_2U.fq.gz" \
      ILLUMINACLIP:${ADAPTERS}:2:28:8 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
  fi

  # paired (_1P/_2P) alignment -- same filters as star_alignm_Vvinif_Pviticola_r1.sh
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

  # F-only (_1U) alignment -- same filters as star1_onlyF_unpaired.sh, single-end
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

  # R-only (_2U) alignment -- NEW (2026-09-04): no existing script in this project
  # has ever mapped reverse-mate-only-survived reads before, only forward-only.
  # Added so R-only fraction (and its own 5'->3' positional distribution) is
  # available per the positional-feature-extraction schema, added before this
  # script's first submission rather than retrofitted after wasting compute.
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

############################ 2. featureCounts -- paired, F-only, and R-only, separately ####
load_subread

tail -n +2 "$RUN_LIST" | while IFS=$'\t' read -r run sample url1 url2; do
  # FIX (2026-09-07): a run skipped earlier for missing raw fastq has no BAM here -
  # featureCounts on a nonexistent file would fail and, under set -e, abort the
  # WHOLE remaining loop. Skip with a warning instead (same fix as Froussios's
  # equivalent loop, same root cause).
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
  if [ ! -f "${STAR_RONLY_DIR}/${run}_Aligned.sortedByCoord.out.bam" ]; then
    echo "WARNING: no R-only BAM for ${run} (likely skipped earlier for missing fastq) -- skipping featureCounts for this run." >&2
  elif [ ! -f "${STAR_RONLY_DIR}/fC/${run}.count" ]; then
    featureCounts -T 24 -a "$ANNOT_GFF3" -t gene -g gene_id \
      -o "${STAR_RONLY_DIR}/fC/${run}.count" \
      "${STAR_RONLY_DIR}/${run}_Aligned.sortedByCoord.out.bam"
  fi
done

# NEGATIVE-LOGIC CHECK: assert expected sample count in all three count-file sets --
# same class of silent gap already caught for Shi2024 (74/83) and the original
# Chitarrini2020 reprocessing (29/33).
N_EXPECTED=$(( $(wc -l < "$RUN_LIST") - 1 ))
N_PAIRED=$(ls "${STAR_PAIRED_DIR}/fC/"*.count 2>/dev/null | wc -l)
N_FONLY=$(ls "${STAR_FONLY_DIR}/fC/"*.count 2>/dev/null | wc -l)
N_RONLY=$(ls "${STAR_RONLY_DIR}/fC/"*.count 2>/dev/null | wc -l)
echo "Expected samples: ${N_EXPECTED}. Paired count files: ${N_PAIRED}. F-only count files: ${N_FONLY}. R-only count files: ${N_RONLY}."
if [ "$N_PAIRED" -ne "$N_EXPECTED" ] || [ "$N_FONLY" -ne "$N_EXPECTED" ] || [ "$N_RONLY" -ne "$N_EXPECTED" ]; then
  echo "WARNING: mismatch -- check which run(s) failed before trusting these results." >&2
fi

echo "Done. Paired counts: ${STAR_PAIRED_DIR}/fC/. F-only counts: ${STAR_FONLY_DIR}/fC/. R-only counts: ${STAR_RONLY_DIR}/fC/."
