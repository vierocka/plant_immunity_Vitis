cd /home/veve/Dropbox/MendelUni_Vinselect/reference/stringDB
count=$(wc -l COG_withVitisVinifera_29760.list | cut -d" " -f1)
for ((i=1; i<$(($count+1)); i++))
do
ID_OG=$(sed -n ''$i'p' COG_withVitisVinifera_29760.list)
grep -P "\t$ID_OG\t" COG.mappings.v12.0.txt >> 29760all.COG.mappings.v12.0.txt
done
awk ' /^3702./ || /^29760./ { print $0 } ' 29760all.COG.mappings.v12.0.txt  > 29760_3702.COG.mappings.v12.0.txt

echo -e "OG\tAthalProteins\tVvinProteins\tAthalAnnotation\tVvinAnnotation" > OGs_29760_3702_stringDB.tsv
for ((i=1; i<$(($count+1)); i++))
do
ID_OG=$(sed -n ''$i'p' COG_withVitisVinifera_29760.list)
Athal=$(grep -P "\t$ID_OG\t" 29760_3702.COG.mappings.v12.0.txt | grep "^3702." | cut -f1 | tr -s '\n' ',')
Nathal=$(grep -P "\t$ID_OG\t" 29760_3702.COG.mappings.v12.0.txt | grep -c "^3702." )
if [[ $Nathal > 0 ]]
then
Vvin=$(grep -P "\t$ID_OG\t" 29760_3702.COG.mappings.v12.0.txt | grep "^29760." | cut -f1 | tr -s '\n' ',')
AthalAnnot=$(grep -P "\t$ID_OG\t" 29760_3702.COG.mappings.v12.0.txt | grep "^3702." | cut -f5 | tr -s '\n' ';')
VvinAnnot=$(grep -P "\t$ID_OG\t" 29760_3702.COG.mappings.v12.0.txt | grep "^29760." | cut -f5 | tr -s '\n' ';')
echo -e $ID_OG"\t"$Athal"\t"$Vvin"\t"$AthalAnnot"\t"$VvinAnnot | sed s/,\\t/\\t/ >> OGs_29760_3702_stringDB.tsv
fi
done

