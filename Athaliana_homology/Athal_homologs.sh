cd ~/Dropbox/MendelUni_Vinselect/reference/Vitis_vinifera
for (( i=1; i<26170; i++ ))
do
IDvv=$( sed -n ''$i'p' 26169genes_with_rlogs_annotation.csv | cut -f1 )
IDathal=$( grep $IDvv qPN40024.v4.pep_refTAIR10_blp.out | head -n 1 | cut -f2)
Athalhomol=$(grep $IDvv qPN40024.v4.pep_refTAIR10_blp.out | head -n 1 | cut -f3)
fullINFOAthal=$( grep -P "$IDathal" TAIR10/GCF_000001735.3_TAIR10_protein.faa | head -n 1 | tr -d '>' )
echo -e $IDvv"\t"$fullINFOAthal"\t"$Athalhomol >> Vvin_AthalTAIR10_homolPerc.csv
done
