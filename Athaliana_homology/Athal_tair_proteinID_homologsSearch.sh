echo "Athal protein ID search script is running." 
for ((i=2; i<26171; i++))
do 
ID=$(sed -n ''$i'p' DGE_Rpv1_vs_WtCultivar_ParamCombat_rlogs_26169genes_sumReadsOver15.csv | cut -f8 | tr -d "\"")
 if [[ $ID != "" ]]
 then 
 Athal=$(grep "29760.$ID," ~/Dropbox/MendelUni_Vinselect/reference/stringDB/OGs_29760_3702_stringDB.tsv | cut -f2 | cut -d"," -f1)
 
 if [[ Athal != "" ]] 
 then
 echo $Athal >> rpv1wt_tair_id.list
 else
 echo "" >> rpv1wt_tair_id.list
 fi
 
 else 
 echo ""  >> rpv1wt_tair_id.list
 fi
 done
 
 echo "DONE!"
