#!/bin/bash -l
#SBATCH --cpus-per-task=24
#SBATCH --mem=192gb
#SBATCH --time=72:00:00
#SBATCH --output=chitarrini2020_mapping_%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=youremail@here
#SBATCH --account=xxx

###############################################################################
# STAR alignment for Chitarrini et al. 2020 (ENA PRJEB28042, 33 runs)
# against the full PN40024 v4 genome (not the original paper's
# transcriptome-only/Bowtie2 approach), so gene IDs land directly in the
# same Vitvi/PN40024 v4 space as the rest of the manuscript. Uniquely-
# mapped reads only (--outFilterMultimapNmax 1), two-pass mode. Run after
# chitarrini2020_filtering_ramses.sh.
###############################################################################

set -euo pipefail

REF_DIR="/scratch/USERNAME/Vitis/GGref"    # reuse the existing PN40024 v4 reference
WORK_DIR="/scratch/USERNAME/Vitis/Chitarrini2020_RNAseq"
READS="/scratch/USERNAME/Vitis/"
RUN_LIST="${WORK_DIR}/chitarrini_run_list.tsv"
GENOME_FA="${REF_DIR}/Vitis_vinifera.PN40024.v4.dna.toplevel.fa"
ANNOT_GFF3="${REF_DIR}/Vitis_vinifera.PN40024.v4.56.gff3"
STAR_INDEX="${REF_DIR}/STAR_index_PN40024v4"

mkdir -p "$WORK_DIR"
cd "$READS"

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
       --sjdbOverhang 49 \
       --genomeSAindexNbases 13
fi

############################ 2. PER-SAMPLE ALIGNMENT ############################
for ((i=2; i<35; i++)); do
  run=$(sed -n ''$i'p' "$RUN_LIST" | cut -f3 | cut -d"/" -f7 | cut -d"_" -f1)

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
done

echo "Done. Sorted BAMs in ${READS}."
