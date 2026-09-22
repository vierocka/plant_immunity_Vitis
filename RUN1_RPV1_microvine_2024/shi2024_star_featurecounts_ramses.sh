#!/bin/bash -l
#SBATCH --cpus-per-task=32
#SBATCH --mem=256gb
#SBATCH --time=120:00:00
#SBATCH --output=shi2024_RNAseq_%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=youremail@here
#SBATCH --account=xxx

###############################################################################
# Reprocessing of the RUN1/RPV1 microvine dataset (Plants 2024, PMC11314213;
# 83 runs across 3 BioProjects -- Syrah PRJNA862686, G5 PRJNA1118503, MV102/
# MV032 PRJNA1121615) with this project's own pipeline (PN40024 v4 genome
# only + STAR + featureCounts), same rationale as the Chitarrini2020
# reprocessing this is adapted from. Vitis-vinifera-only mapping is
# deliberate: reads from genuinely novel introgressed sequence simply won't
# align -- this reprocessing is host-genome-background-focused by design.
#
# Differences from chitarrini2020_star_featurecounts_ramses.sh:
#   - fastq already downloaded separately (sra_download_ramses.sh) into
#     ${WORK_DIR}/fastq/ as SRR*_1.fastq.gz / SRR*_2.fastq.gz -- no download
#     step here.
#   - Builds a SEPARATE STAR index with sjdbOverhang=149, matching this
#     dataset's actual 2x150bp reads (Illumina HiSeq3000, per the paper's
#     M&M) -- to omit the low-mappability issue a mismatched overhang would
#     cause. Saved to a new directory, NOT overwriting the existing
#     STAR_index_PN40024v4 (sjdbOverhang=49, built for Chitarrini's 2x50bp
#     reads and shared with other reprocessing in this project). Built once,
#     skipped on resubmission if already present.
#   - Per-tool module unload/load kept INSIDE the per-sample loop, same as
#     the original -- confirmed working on Ramses; each tool needs its own
#     compiler swapped in immediately before it runs, not loaded once
#     up front for a whole pass.
#   - Added skip-if-already-done and skip-if-fastq-not-yet-downloaded checks,
#     since fastq are still arriving from the download job as of this
#     writing -- safe to resubmit this whole script as more accessions land,
#     it will only process what's new.
#   - featureCounts runs across whatever BAMs exist at run time (not
#     necessarily all 83) and reports how many are missing, rather than
#     failing outright.
###############################################################################

set -euo pipefail

REF_DIR="/scratch/USERNAME/Vitis/GGref"
WORK_DIR="/scratch/USERNAME/Vitis/Shi2024"
FASTQ_DIR="${WORK_DIR}/fastq"
ACC_LIST="${WORK_DIR}/sra_accession_list.txt"
GENOME_FA="${REF_DIR}/Vitis_vinifera.PN40024.v4.dna.toplevel.fa"
ANNOT_GFF3="${REF_DIR}/Vitis_vinifera.PN40024.v4.56.gff3"
STAR_INDEX="${REF_DIR}/STAR_index_PN40024v4_ohd149"    # 150bp-specific index, NOT the shared 50bp one

mkdir -p "$WORK_DIR"
cd "$WORK_DIR"

############################ 0. STAR GENOME INDEX, sjdbOverhang=149 (build once) #
if [ ! -f "${STAR_INDEX}/SAindex" ]; then
    echo "=== building STAR index (sjdbOverhang=149, for 2x150bp reads) $(date) ==="
    mkdir -p "$STAR_INDEX"
    module unload compiler/GCC 2>/dev/null || true
    module add bio/STAR/2.7.10b-GCC-11.3.0

    STAR --runThreadN 32 \
         --runMode genomeGenerate \
         --genomeDir "$STAR_INDEX" \
         --genomeFastaFiles "$GENOME_FA" \
         --sjdbGTFfile "$ANNOT_GFF3" \
         --sjdbGTFtagExonParentTranscript Parent \
         --sjdbGTFfeatureExon exon \
         --sjdbOverhang 149 \
         --genomeSAindexNbases 13
         # sjdbOverhang 149 = (read length 150) - 1, matches this dataset's
         # actual 2x150bp reads. genomeSAindexNbases 13 unchanged from the
         # existing 50bp index (same genome, same recommended value for a
         # ~500 Mbp grapevine genome).
else
    echo "SKIP index build: ${STAR_INDEX}/SAindex already exists"
fi

############################ PER-SAMPLE: Q30 filter -> STAR align ##############
while read -r RUN; do
    R1="${FASTQ_DIR}/${RUN}_1.fastq.gz"
    R2="${FASTQ_DIR}/${RUN}_2.fastq.gz"
    OUT1="${WORK_DIR}/${RUN}_1.q30.fastq.gz"
    OUT2="${WORK_DIR}/${RUN}_2.q30.fastq.gz"
    BAM="${WORK_DIR}/${RUN}_Aligned.sortedByCoord.out.bam"

    if [[ -f "$BAM" ]]; then
        echo "SKIP ${RUN}: BAM already exists"
        continue
    fi
    if [[ ! -f "$R1" || ! -f "$R2" ]]; then
        echo "SKIP ${RUN}: fastq not yet downloaded"
        continue
    fi

    if [[ ! -f "$OUT1" || ! -f "$OUT2" ]]; then
        # ONE quality filter, nothing else: average read quality >= Q30.
        # Every other fastp filter explicitly disabled -- matches
        # Chitarrini2020 reprocessing exactly, for methods consistency
        # across every external dataset this project reprocesses.
        echo "=== fastp ${RUN} $(date) ==="
        module unload compiler/GCC 2>/dev/null || true
        module add bio/fastp/1.0.1-GCC-13.3.0

        fastp \
            -i "$R1" -I "$R2" \
            -o "$OUT1" -O "$OUT2" \
            --average_qual 30 \
            --disable_length_filtering \
            -w 8 \
            -j "${WORK_DIR}/${RUN}_fastp.json" -h "${WORK_DIR}/${RUN}_fastp.html"
    fi

    echo "=== STAR ${RUN} $(date) ==="
    module unload compiler/GCC 2>/dev/null || true
    module add bio/STAR/2.7.10b-GCC-11.3.0

    STAR --runThreadN 32 \
         --genomeDir "$STAR_INDEX" \
         --readFilesIn "$OUT1" "$OUT2" \
         --readFilesCommand zcat \
         --outFileNamePrefix "${WORK_DIR}/${RUN}_" \
         --outSAMtype BAM SortedByCoordinate \
         --outBAMsortingBinsN 20 \
         --limitBAMsortRAM 10000000000 \
         --outFilterMultimapNmax 1 \
         --twopassMode Basic
         # --outFilterMultimapNmax 1 = uniquely-mapped reads only (excluded
         # from the BAM entirely, not down-weighted). All other STAR filters
         # left at default, matching Chitarrini2020 reprocessing exactly.
done < "$ACC_LIST"

############################ SINGLE featureCounts CALL, all available BAMs #####
module unload compiler/GCC 2>/dev/null || true
module add bio/Subread/2.1.1-GCC-13.2.0

BAMS=()
MISSING=0
while read -r RUN; do
    BAM="${WORK_DIR}/${RUN}_Aligned.sortedByCoord.out.bam"
    if [[ -f "$BAM" ]]; then
        BAMS+=("$BAM")
    else
        MISSING=$((MISSING + 1))
    fi
done < "$ACC_LIST"

echo "featureCounts input: ${#BAMS[@]} BAMs present, ${MISSING} missing/not-yet-aligned"

featureCounts -T 32 -a "$ANNOT_GFF3" -p -B -C -t gene -g gene_id \
    -o "${WORK_DIR}/shi2024_all_samples.counts.tsv" \
    "${BAMS[@]}"

echo "Done. Combined count matrix: ${WORK_DIR}/shi2024_all_samples.counts.tsv"
