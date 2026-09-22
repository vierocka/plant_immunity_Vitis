#!/bin/bash -l
#SBATCH --cpus-per-task=32
#SBATCH --mem=64gb
#SBATCH --time=6:00:00
#SBATCH --output=froussios2019_STAR_featureCounts.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=youremail@here
#SBATCH --account=xxx

REF_DIR="/scratch/USERNAME/Athal/TAIR10ref"
WORK_DIR="/scratch/USERNAME/Athal"
FASTQ_DIR="${WORK_DIR}/fastq"                        # from froussios2019_download_ramses.sh
RUN_LIST="${WORK_DIR}/froussios2019_run_list.tsv"     # run_accession, sample_title, experiment_batch, ercc_mix, url_R1, url_R2, md5_R1, md5_R2
GENOME_FA="${REF_DIR}/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
ANNOT_GFF3="${REF_DIR}/Arabidopsis_thaliana.TAIR10.63.gff3"
STAR_INDEX="${REF_DIR}/STAR_index_TAIR10"

cd "$WORK_DIR"

module unload compiler/GCC
module add bio/Subread/2.1.1-GCC-13.2.0

BAMS=$(tail -n +2 "$RUN_LIST" | while IFS=$'\t' read -r run sample exp mix url1 url2 md5_1 md5_2; do
  echo -n "${WORK_DIR}/${run}_${sample}/${run}_Aligned.sortedByCoord.out.bam "
done | cut -d" " -f1-14)

featureCounts \
  -T 32 \
  -p -B -C \
  -t gene \
  -g gene_id \
  -a "$ANNOT_GFF3" \
  -o "${WORK_DIR}/froussios2019_samples_from1_to14_perGene.counts.tsv" \
  $BAMS


exit 0
