### FastQC v0.12.1 - default settings
## to check raw reads quality and set trimming
## run in a loop for read fq files (36 x2 files, followed by visual inspection)
# code example:
fastqc -q -t 8 -o $outDir $myDir/$ID"_1P.fq.gz"

### Trimmomatic v0.39
## http://www.usadellab.org/cms/?page=trimmomatic
## to trim and filter raw reads
## run in a loop for all read pairs
# code example for two paired files:
trimmomatic PE -threads 8 -trimlog $outDir/$ID".log" -summary $outDir/$ID"_summary.txt" -quiet -validatePairs $inDir/$ID"_1.fq.gz" $inDir/$ID"_2.fq.gz" $outDir/$ID"_1P.fq.gz" $outDir/$ID"_1U.fq.gz" $outDir/$ID"_2P.fq.gz" $outDir/$ID"_2U.fq.gz" ILLUMINACLIP:$myPath/SW/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa:2:28:8 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36

### STAR aligner version 2.7.4a
## to map trimmed reads onto the reference (Vitis vinifera PN40024.v4; genome-date 2021-05) 
## https://grapedia.org/files-download/
## https://academic.oup.com/g3journal/article/13/5/jkad067/7086178
# code to build the reference database:
cd $myPath/SW/STAR-2.7.10b/bin/Linux_x86_64
./STAR --runMode genomeGenerate --runThreadN 8 --genomeDir $myPath/reference/STARref --genomeSAindexNbases 12 --genomeFastaFiles $myPath/reference/Vitis_vinifera.fa
## run in a loop for all read pairs
# code example to map filtered reads, with --BySJout option to stricter filter splice junctions based on splice junction quality and annotation support:
./STAR --runThreadN 8 --genomeDir $myPath/reference/STARref --outFilterType BySJout --outFilterMultimapNmax 10 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 12 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --outFileNamePrefix $starOut/$ID"_" --outReadsUnmapped Fastx --genomeLoad LoadAndRemove --limitBAMsortRAM 12000000000 --outSAMtype BAM SortedByCoordinate --readFilesCommand gunzip -c --readFilesIn $outDir/$ID"_1P.fq.gz" $outDir/$ID"_2P.fq.gz"

### subread version 2.0.3-Linux-x86_64; function featureCounts
## to count mapped reads 
## run in a loop for all read pairs
# code example:
./featureCounts -T 4 -p -F SAF -a $myPath/reference/Vv_genes.saf -o $fCout/$ID".count" $starOut/$ID"_Aligned.sortedByCoord.out.bam"
