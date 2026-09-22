#!/bin/bash -l
#SBATCH --cpus-per-task=24
#SBATCH --mem=128gb
#SBATCH --time=72:00:00
#SBATCH --output=froussios2019_mapping_%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=youremail@here
#SBATCH --account=xxx

###############################################################################
# STAR alignment for Froussios et al. 2019 (ArrayExpress E-MTAB-5446 / ENA
# ERP021226, 17 runs) against the Ensembl Plants TAIR10 reference, same
# pipeline conventions used throughout this project. Uniquely-mapped reads
# only (--outFilterMultimapNmax 1), two-pass mode. Run after
# froussios2019_filtering_ramses.sh.
###############################################################################

set -euo pipefail

REF_DIR="/scratch/USERNAME/Athal/TAIR10ref"
WORK_DIR="/scratch/USERNAME/Athal"
RUN_LIST="${WORK_DIR}/froussios2019_run_list.tsv"
GENOME_FA="${REF_DIR}/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
ANNOT_GFF3="${REF_DIR}/Arabidopsis_thaliana.TAIR10.63.gff3"
STAR_INDEX="${REF_DIR}/STAR_index_TAIR10"

cd "$WORK_DIR"

module unload compiler/GCC 2>/dev/null || true
module add bio/STAR/2.7.10b-GCC-11.3.0

############################ 1. STAR GENOME INDEX (build once) #################
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

############################ 2. PER-SAMPLE ALIGNMENT ############################
tail -n +2 "$RUN_LIST" | while IFS=$'\t' read -r run sample exp mix url1 url2 md5_1 md5_2; do
  sample_dir="${WORK_DIR}/${run}_${sample}"
  cd "$sample_dir"

  if [ ! -f "${run}_Aligned.sortedByCoord.out.bam" ]; then
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

echo "Done. Sorted BAMs in each \${WORK_DIR}/<run>_<sample>/ subfolder."
