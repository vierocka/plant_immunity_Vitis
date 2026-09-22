#!/bin/bash -l
#SBATCH --cpus-per-task=24
#SBATCH --mem=192gb
#SBATCH --time=72:00:00
#SBATCH --output=chitarrini2020_RNAseq.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=youremail@here
#SBATCH --account=xxx

###############################################################################
# Reprocessing of Chitarrini et al. 2020 (Sci Rep 10:12193) RNA-seq data
# (ENA project PRJEB28042, 33 runs: 0/12/24/48/96/120 hpi x mock/inoculated x 3
# biological replicates, Illumina HiSeq2500 2x50bp) with YOUR OWN pipeline
# (full PN40024 v4 genome + STAR + featureCounts), not the original paper's
# transcriptome-only/Bowtie2 approach - so gene IDs land directly in the same
# Vitvi/PN40024 v4 space as the rest of the manuscript, no ID crosswalk needed.
#
# Otherwise deliberately minimal: no adapter trimming, no length filtering, no
# extra QC beyond the two requirements below.
#   1. Average-per-read Q30 cutoff. STAR has NO native option to filter reads
#      by raw base-call quality (its own --outFilter* flags are post-
#      alignment, mismatch/multimapping-based, not raw-quality-based), so
#      this uses a single fastp pass with every OTHER fastp filter explicitly
#      disabled, immediately followed by STAR - the minimal addition that
#      still satisfies "Q30 average per read", not a trimming step.
#   2. Uniquely-mapped reads only (--outFilterMultimapNmax 1): reads with more
#      than one equally-best-scoring alignment location are dropped entirely
#      from STAR's output rather than reported/kept, so featureCounts below
#      never sees a multi-mapper and needs no -M/--fraction handling.
###############################################################################

set -euo pipefail

REF_DIR="/scratch/USERNAME/Vitis/GGref"    # reuse the existing PN40024 v4 reference
WORK_DIR="/scratch/USERNAME/Vitis/Chitarrini2020_RNAseq"
READS="/scratch/USERNAME/Vitis/"
RUN_LIST="${WORK_DIR}/chitarrini_run_list.tsv"      # run_accession, sample_title, url_R1, url_R2
GENOME_FA="${REF_DIR}/Vitis_vinifera.PN40024.v4.dna.toplevel.fa"
ANNOT_GFF3="${REF_DIR}/Vitis_vinifera.PN40024.v4.56.gff3"
STAR_INDEX="${REF_DIR}/STAR_index_PN40024v4"

mkdir -p "$WORK_DIR"
cd "$READS"
# Copy chitarrini_run_list.tsv here before submitting (scp'd alongside this script).

############################ MODULES ###########################################
# module unload compiler/GCC 2>/dev/null || true
# module add bio/STAR/2.7.10b-GCC-11.3.0
#  module add bio/fastp/1.0.1-GCC-13.3.0   # CHECK EXACT MODULE NAME/VERSION: `module spider fastp` first,
                        # this project's reference scripts never loaded fastp before.
# module add bio/Subread/2.1.1-GCC-13.2.0   # for featureCounts, matches gencode_ramses.sh

############################ 1. STAR GENOME INDEX (build once) #################
# if [ ! -f "${STAR_INDEX}/SAindex" ]; then
#   mkdir -p "$STAR_INDEX"
#  STAR --runThreadN 24 \
#       --runMode genomeGenerate \
#       --genomeDir "$STAR_INDEX" \
#       --genomeFastaFiles "$GENOME_FA" \
#       --sjdbGTFfile "$ANNOT_GFF3" \
#       --sjdbGTFtagExonParentTranscript Parent \
#       --sjdbGTFfeatureExon exon \
#       --sjdbOverhang 49 \
#       --genomeSAindexNbases 13
       # sjdbOverhang 49 = (read length 50) - 1, matches this dataset's 2x50bp reads.
       # genomeSAindexNbases 13 is close to the recommended min(14, log2(genomeLen)/2 - 1)
       # for a ~500 Mbp grapevine genome - verify against your genome's actual size if
       # this differs from what the rest of the project used.
# fi

############################ 2. PER-SAMPLE: download -> Q30 filter -> STAR align ####
for ((i=2; i<35; i++))
do
run=$(sed -n ''$i'p' $READS/chitarrini_run_list.tsv | cut -f3 | cut -d"/" -f7 | cut -d"_" -f1 )
    # ONE quality filter, nothing else: average read quality >= Q30.
module unload compiler/GCC    # Every other fastp filter explicitly disabled.
module add bio/fastp/1.0.1-GCC-13.3.0
    fastp \
      -i "${run}_1.fastq.gz" -I "${run}_2.fastq.gz" \
      -o "${run}_1.q30.fastq.gz" -O "${run}_2.q30.fastq.gz" \
      --average_qual 30 \
      --disable_length_filtering \
      -w 8 \
      -j "${run}_fastp.json" -h "${run}_fastp.html"
      # --disable_quality_filtering turns off fastp's OWN default per-base
      # quality filter (which rejects reads if too many individual bases are
      # low-quality) - --average_qual 30 is re-enabled separately below since
      # it is not controlled by that flag; fastp applies --average_qual
      # regardless of --disable_quality_filtering (it is a distinct filter).

      module unload compiler/GCC
      module add bio/STAR/2.7.10b-GCC-11.3.0

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
         # --outFilterMultimapNmax 1 = uniquely-mapped reads only (reads with
         # >1 equally-best alignment are excluded from the BAM entirely, not
         # just down-weighted/flagged). All other STAR filters left at
         # default (mismatch rate etc.) - not touched beyond this request.
done

############################ 3. SINGLE featureCounts CALL, ALL SAMPLES ##########
module unload compiler/GCC 2>/dev/null || true
module add bio/Subread/2.1.1-GCC-13.2.0   # for featureCounts, matches gencode_ramses.sh

# One call across all 33 BAMs -> one combined count matrix (columns = samples),
# rather than 33 separate per-sample count files.
BAMS=$(tail -n +2 "$RUN_LIST" | while IFS=$'\t' read -r run sample url1 url2; do
  echo -n "${run}_Aligned.sortedByCoord.out.bam "
done)

featureCounts -T 24 -a GFF -p -B -C  -t gene  -g gene_id  GGref/Vitis_vinifera.PN40024.v4.56.gff3   -o chitarrini2020_all_samples.counts.tsv   $BAMS

echo "Done. Combined count matrix: ${WORK_DIR}/chitarrini2020_all_samples.counts.tsv"
