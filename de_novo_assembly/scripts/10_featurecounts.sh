cd ~/Desktop/SW/subread-2.0.3-Linux-x86_64/bin
inDir="~/Dropbox/MendelUni_Vinselect/spades/unmapped/STARmap2"
mkdir ~/Dropbox/MendelUni_Vinselect/spades/unmapped/STARmap2/fC

for ((i=1; i<37; i++))
do
ID=$(ls -1 $inDir | grep "_Aligned.sortedByCoord.out.bam" | sed -n ''$i'p' | cut -d"_" -f1)
./featureCounts -T 4 -p -F SAF -a ~/Dropbox/MendelUni_Vinselect/spades/unmapped/hardF_myFilter_all_cap.saf -o $inDir/fC/$ID".count" $inDir/$ID"_Aligned.sortedByCoord.out.bam"
done
