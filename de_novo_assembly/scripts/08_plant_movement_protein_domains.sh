cd ~/Dropbox/MendelUni_Vinselect/spades/unmapped

counts=$(grep -P "Plant viral movement protein" hardF_myFilter_all_cap_augustus_7658proteins.hmmscan.tblout  | awk ' { print $3 } ' | sort | uniq | wc -l | cut -d" " -f1)
for ((i=1; i<$(($counts+1)); i++))
do
ID=$(grep -P "Plant viral" hardF_myFilter_all_cap_augustus_7658proteins.hmmscan.tblout  | awk ' { print $3 } ' | sort | uniq | sed -n ''$i'p')
#OGn=$(grep -Pn "$ID\b|$ID\t|$ID$" hardF_myFilter_all_cap_augustus_proteins.bdBlP.I2mcl | cut -d":" -f1)
#grep "OG_$OGn," pattern_boolean.csv >> PlantViralMovementProtein_boolean.csv
#grep "OG_$OGn," pattern_counts.csv >> PlantViralMovementProtein_count.csv
contig=$(grep -P "Protein=$ID$" hardF_myFilter_all_cap_augustus_9columns.gff3 | cut -f1) 
echo -e $ID"\t"$contig >>  PlantViralMovementProtein_protein_contigIDs.csv
done
#sort PlantViralMovementProtein_boolean.csv | uniq
cut -f2 PlantViralMovementProtein_protein_contigIDs.csv | sort | uniq > PlantViralMovementProtein_contigIDs.csv


cPMPr=$(wc -l PlantViralMovementProtein_contigIDs.csv | cut -d" " -f1)
echo $cPMPr
head -n 1 STARmap3_n100/fC/Rlogs_deNovo_4036transcripts_more15reads_combat.csv >  STARmap3_n100/fC/PlantViralMovementProtein_rlogCombat.csv
for ((i=1; i<$(($cPMPr+1)); i++))
do
ContigID=$(sed -n ''$i'p' PlantViralMovementProtein_contigIDs.csv)
# awk ' BEGIN{FS="\t"}; /"'$ContigID'"/ { if ($5 <= 0.05 ) { print FILENAME"\t"$0 }} ' STARmap3_n100/fC/Denovo_*.csv
awk ' BEGIN{FS="\t"}; /"'$ContigID'"/ { if ($5 <= 0.05 ) { print FILENAME"\t"$0 }} ' STARmap3_n100/fC/Rlogs_deNovo_4036transcripts_more15reads_combat.csv >>  STARmap3_n100/fC/PlantViralMovementProtein_rlogCombat.csv
done
