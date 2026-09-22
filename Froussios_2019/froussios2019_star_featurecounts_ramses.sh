#!/bin/bash -l
#SBATCH --cpus-per-task=24
#SBATCH --mem=128gb
#SBATCH --time=72:00:00
#SBATCH --output=froussios2019_STAR_featureCounts.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=youremail@here
#SBATCH --account=xxx

###############################################################################
# Reprocessing of Froussios, Schurch, Mackinnon, Gierlinski, Duc, Simpson,
# Barton (2019), Bioinformatics 35(18):3372-3377 (ArrayExpress E-MTAB-5446 /
# ENA study ERP021226), 17 runs, WT Arabidopsis thaliana Col-0 seedlings,
# Illumina HiSeq 2000 2x101bp -- with THIS PROJECT'S OWN pipeline (STAR +
# featureCounts against Ensembl Plants TAIR10), for direct comparability with
# the 03_AED / Chitarrini2020 / Shi2024 AggrDiv and gene-concentration
# methodology (same statistic, same normalization family, same decomposition).
#
# Purpose: independent, genuinely isogenic (Col-0 wild type, one lab, one
# growth batch of protocols) same-condition replicate set, to test whether
# the ~12% "genes needed for 75% of divergence" concentration figure seen
# BOTH in the main study's real cross-genotype/within-genotype AED tests AND
# in the 03_AED pure-Gaussian-noise negative control (AED_noise_negative_
# control.R, itself run on Susceptible's own triplicate -- i.e. already one
# genotype, no cross-background confound) is a generic property of this
# divergence decomposition under ANY replicate variation, or whether it is
# specifically shaped by the non-isogenic breeding-line backgrounds used
# throughout the rest of this project.
#
# ALL METADATA BELOW VERIFIED DIRECTLY against the ENA/BioStudies APIs on
# 2026-08-17 (organism, ecotype, genotype, layout, instrument, read counts) --
# see Froussios_2019/NOTES.md. Ensembl Plants release 63 (current as of
# 2026-08-17, confirmed via `curl https://rest.ensembl.org/info/eg_version`)
# used for the reference; re-verify the release number if run much later.
#
# Otherwise deliberately mirrors chitarrini2020_star_featurecounts_ramses.sh
# exactly (same fastp Q30-only filter, same STAR uniquely-mapped-only +
# two-pass settings, same single combined featureCounts call) so that any
# difference in downstream results traces to the biology/dataset, not a
# pipeline discrepancy.
###############################################################################

set -euo pipefail

REF_DIR="/scratch/USERNAME/Athal/TAIR10ref"
WORK_DIR="/scratch/USERNAME/Athal"
FASTQ_DIR="${WORK_DIR}/fastq"                        # from froussios2019_download_ramses.sh
RUN_LIST="${WORK_DIR}/froussios2019_run_list.tsv"     # run_accession, sample_title, experiment_batch, ercc_mix, url_R1, url_R2, md5_R1, md5_R2
GENOME_FA="${REF_DIR}/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
ANNOT_GFF3="${REF_DIR}/Arabidopsis_thaliana.TAIR10.63.gff3"
STAR_INDEX="${REF_DIR}/STAR_index_TAIR10"

mkdir -p "$WORK_DIR" "$REF_DIR"
cd "$WORK_DIR"
# Copy froussios2019_run_list.tsv here before submitting (scp'd alongside this script);
# fastq/ should already be populated by froussios2019_download_ramses.sh.

############################ 0. REFERENCE (download once, verified URLs) #######
# if [ ! -f "$GENOME_FA" ]; then
#  curl -sS -o "${GENOME_FA}.gz" \
#    "https://ftp.ebi.ac.uk/ensemblgenomes/pub/plants/release-63/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz"
#  gunzip "${GENOME_FA}.gz"
# fi
# if [ ! -f "$ANNOT_GFF3" ]; then
#  curl -sS -o "${ANNOT_GFF3}.gz" \
#    "https://ftp.ebi.ac.uk/ensemblgenomes/pub/plants/release-63/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.63.gff3.gz"
#  gunzip "${ANNOT_GFF3}.gz"
# fi

############################ MODULES (loaded per-tool-call, not once up front) ##
# STAR (2.7.10b) is built against GCC/11.3.0; fastp (1.0.1) is built against
# GCC/13.3.0. Lmod is hierarchical and refuses to have two different
# compiler/GCC versions loaded at once (confirmed on Ramses: loading fastp
# right after STAR fails with "Module cannot be loaded due to a conflict" on
# compiler/GCC/13.3.0). Loading both once at the top -- as
# chitarrini2020_star_featurecounts_ramses.sh does -- only worked there
# because that script's STAR/fastp/Subread modules all happened to share one
# GCC version; Arabidopsis's module set doesn't. Fix: never hold both loaded
# together -- unload compiler/GCC (which Lmod cascades into unloading every
# dependent, e.g. the previously-loaded tool) immediately before loading
# whichever tool is needed next, every single time it's needed, including
# once per sample inside the loop below (fastp and STAR alternate there).
load_fastp() { module unload compiler/GCC 2>/dev/null || true; module add bio/fastp/1.0.1-GCC-13.3.0; }  # CHECK EXACT MODULE NAME/VERSION: `module spider fastp` first.
load_star()  { module unload compiler/GCC 2>/dev/null || true; module add bio/STAR/2.7.10b-GCC-11.3.0; }

############################ 1. STAR GENOME INDEX (build once) #################
load_star
#if [ ! -f "${STAR_INDEX}/SAindex" ]; then
#  mkdir -p "$STAR_INDEX"
#  STAR --runThreadN 24 \
#       --runMode genomeGenerate \
#       --genomeDir "$STAR_INDEX" \
#       --genomeFastaFiles "$GENOME_FA" \
#       --sjdbGTFfile "$ANNOT_GFF3" \
#       --sjdbGTFtagExonParentTranscript Parent \
#       --sjdbGTFfeatureExon exon \
#       --sjdbOverhang 100 \
#       --genomeSAindexNbases 12
       # sjdbOverhang 100 = (read length 101) - 1, per ENA-verified 2x101bp reads
       # (base_count/read_count/2 = ~101 across all 17 runs, checked 2026-08-17).
       # genomeSAindexNbases 12, per the formula check printed above for the
       # ~120 Mbp TAIR10 genome (Vitis's PN40024 v4 used 13 for its ~500 Mbp
       # genome; Arabidopsis is smaller so a lower value is expected and correct,
       # not a copy-paste oversight).
# fi

############################ 2. PER-SAMPLE: Q30 filter -> STAR align ###########
tail -n +2 "$RUN_LIST" | while IFS=$'\t' read -r run sample exp mix url1 url2 md5_1 md5_2; do
  sample_dir="${WORK_DIR}/${run}_${sample}"
  mkdir -p "$sample_dir"
  cd "$sample_dir"

  [ -f "${FASTQ_DIR}/${run}_1.fastq.gz" ] || { echo "Missing ${FASTQ_DIR}/${run}_1.fastq.gz -- run froussios2019_download_ramses.sh first" >&2; exit 1; }
  [ -f "${FASTQ_DIR}/${run}_2.fastq.gz" ] || { echo "Missing ${FASTQ_DIR}/${run}_2.fastq.gz -- run froussios2019_download_ramses.sh first" >&2; exit 1; }
  [ -f "${run}_1.fastq.gz" ] || ln -s "${FASTQ_DIR}/${run}_1.fastq.gz" "${run}_1.fastq.gz"
  [ -f "${run}_2.fastq.gz" ] || ln -s "${FASTQ_DIR}/${run}_2.fastq.gz" "${run}_2.fastq.gz"

  if [ ! -f "${run}_Aligned.sortedByCoord.out.bam" ]; then
    # ONE quality filter, nothing else: average read quality >= Q30 -- identical
    # convention to chitarrini2020_star_featurecounts_ramses.sh.
    load_fastp
    fastp \
      -i "${run}_1.fastq.gz" -I "${run}_2.fastq.gz" \
      -o "${run}_1.q30.fastq.gz" -O "${run}_2.q30.fastq.gz" \
      --average_qual 30 \
      --disable_adapter_trimming \
      --disable_length_filtering \
      --disable_quality_filtering \
      -w 8 \
      -j "${run}_fastp.json" -h "${run}_fastp.html"

    load_star
    STAR --runThreadN 24 \
         --genomeDir "$STAR_INDEX" \
         --readFilesIn "${run}_1.q30.fastq.gz" "${run}_2.q30.fastq.gz" \
         --readFilesCommand zcat \
         --outFileNamePrefix "${run}_" \
         --outSAMtype BAM SortedByCoordinate \
         --outBAMsortingBinsN 20 \
         --limitBAMsortRAM 10000000000 \
         --outFilterMultimapNmax 1 \
         --twopassMode Basic
  fi
  cd "$WORK_DIR"
done

############################ 3. SINGLE featureCounts CALL, ALL SAMPLES ##########
module unload compiler/GCC 2>/dev/null || true
module add bio/Subread/2.1.1-GCC-13.2.0

BAMS=$(tail -n +2 "$RUN_LIST" | while IFS=$'\t' read -r run sample exp mix url1 url2 md5_1 md5_2; do
  echo -n "${WORK_DIR}/${run}_${sample}/${run}_Aligned.sortedByCoord.out.bam "
done)

featureCounts \
  -T 24 \
  -p -B -C \
  -t exon \
  -g exon_id \
  -a "$ANNOT_GFF3" \
  -o "${WORK_DIR}/froussios2019_all_samples.counts.tsv" \
  $BAMS

# NEGATIVE-LOGIC CHECK: assert the output actually has all 17 sample columns
# -- don't just assume the loop above produced every BAM. This exact failure
# mode (fewer BAMs than expected reaching featureCounts silently) is what
# produced the Shi2024 74/83 and Chitarrini2020 29/33 gaps caught elsewhere in
# this project; catch it here immediately instead of downstream in R.
N_SAMPLE_COLS=$(($(head -2 "${WORK_DIR}/froussios2019_all_samples.counts.tsv" | tail -1 | awk -F'\t' '{print NF}') - 6))
echo "Sample columns in count matrix: ${N_SAMPLE_COLS} (expected 17)"
if [ "$N_SAMPLE_COLS" -ne 17 ]; then
  echo "WARNING: expected 17 BAMs/columns in the count matrix, found ${N_SAMPLE_COLS}. Check which run(s) failed to align before trusting this matrix." >&2
fi

echo "Done. Combined count matrix: ${WORK_DIR}/froussios2019_all_samples.counts.tsv"
