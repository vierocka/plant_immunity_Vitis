cd ~/Dropbox/MendelUni_Vinselect/STARmap1
for ((i=1; i<37; i++))
do
# raw outputs from STAR aligner
fileID=$(ls -1 | grep "_SJ.out.tab" | sed -n ''$i'p')
ID=$(echo $fileID | cut -d"_" -f1)
myTotal=$(wc -l $fileID | cut -d" " -f1)
# column 7 - Unique read count - Number of uniquely mapped reads supporting this junction
# column 9 - Maximum spliced alignment overhang - Longest length of perfectly matched bases on either side of junction (important for confidence)
myCo=$(awk ' { if ( $7 > 10 && $9 > 20 ) { print $0 }}' $fileID | wc -l)
allUniqMapR=$(grep "Uniquely mapped reads number" $ID"_Log.final.out" | cut -d"|" -f2 | tr -d ' ')
# overview: sample, total count of filtred SJs, total count of raw SJs, count of uniquely mapped reads (library size) 
echo -e $ID"\t"$myCo"\t"$myTotal"\t"$allUniqMapR  >> SJ_ID_filteredCount_allDetected_allUniqMapInSTARalign.csv
done
