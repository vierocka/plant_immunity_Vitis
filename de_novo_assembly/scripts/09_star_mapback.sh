cd ~/Desktop/SW/STAR-2.7.10b/bin/Linux_x86_64
# version 2.7.4a

#./STAR --runMode genomeGenerate --runThreadN 4 --genomeDir ~/Dropbox/MendelUni_Vinselect/spades/unmapped/STARref --genomeSAindexNbases 11 --genomeFastaFiles ~/Dropbox/MendelUni_Vinselect/spades/unmapped/hardF_myFilter_all_cap.fa

mkdir ~/Dropbox/MendelUni_Vinselect/spades/unmapped/STARmap3_n100
readsInF="~/Dropbox/MendelUni_Vinselect/STARmap2"
refInF="~/Dropbox/MendelUni_Vinselect/spades/unmapped/STARref"
outInF="~/Dropbox/MendelUni_Vinselect/spades/unmapped/STARmap3_n100"

count=$(ls -1 $readsInF/ | grep -c "Unmapped.out.mate1$")
for ((i=1; i<$(($count+1)); i++))
do
ID=$(ls -1 $readsInF/ | grep "Unmapped.out.mate1$" |cut -d"_" -f1 | sed -n ''$i'p')
./STAR --runThreadN 8 --genomeDir $refInF/ --outFilterType BySJout --outFilterMultimapNmax 100 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 12 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --outFileNamePrefix $outInF/$ID"_" --outReadsUnmapped Fastx --genomeLoad LoadAndRemove --limitBAMsortRAM 12000000000 --outSAMtype BAM SortedByCoordinate --readFilesIn $readsInF/$ID"_Unmapped.out.mate1" $readsInF/$ID"_Unmapped.out.mate2"
done

cd ~/Desktop/SW/subread-2.0.3-Linux-x86_64/bin
inDir="~/Dropbox/MendelUni_Vinselect/spades/unmapped/STARmap3_n100"
mkdir ~/Dropbox/MendelUni_Vinselect/spades/unmapped/STARmap3_n100/fC

for ((i=1; i<37; i++))
do
ID=$(ls -1 $inDir | grep "_Aligned.sortedByCoord.out.bam" | sed -n ''$i'p' | cut -d"_" -f1)
./featureCounts -T 4 -p -F SAF -a ~/Dropbox/MendelUni_Vinselect/spades/unmapped/hardF_myFilter_all_cap.saf -o $inDir/fC/$ID".count" $inDir/$ID"_Aligned.sortedByCoord.out.bam"
done

exit 0
