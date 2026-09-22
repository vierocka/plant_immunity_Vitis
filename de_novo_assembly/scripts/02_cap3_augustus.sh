# AUGUSTUS (3.4.0) is a gene prediction tool.
# Sources and documentation at https://github.com/Gaius-Augustus/Augustus

cd ~/Dropbox/MendelUni_Vinselect/spades/unmapped/
augustus --strand=both --genemodel=complete --alternatives-from-sampling=true --species=cacao  --gff3=on hardF_myFilter_all_cap.fa  > hardF_myFilter_all_cap_augustus.gff3

cat hardF_myFilter_all_cap_augustus.gff3 | tr -d '\n' | grep -oP "protein sequence = [[A-Z\ #]{1,6000}]" | tr -d '#' | tr -d ' ' | tr -d ']' | tr -d '[' | awk ' BEGIN {RS="proteinsequence="}; NR>1 { print ">pr"NR-1"\n"$0 } ' > hardF_myFilter_all_cap_augustus_proteins.fa
grep -c ">" hardF_myFilter_all_cap_augustus_proteins.fa
