cd ~/Dropbox/MendelUni_Vinselect/spades/unmapped/Rpv1/rna
cat Rpv1_hardF_myFilter_augustus.gff3 | tr -d '\n' | grep -oP "protein sequence = [[A-Z\ #]{1,6000}]" | tr -d '#' | tr -d ' ' | tr -d ']' | tr -d '[' | awk ' BEGIN{RS="proteinsequence="}; NR>1 { print ">Rpv1_pr"NR-1"\n"$0}' > Rpv1-deNovo-extracted_proteins.fa
cd ~/Dropbox/MendelUni_Vinselect/spades/unmapped/Rpv1_12/rna
cat Rpv1-12_hardF_myFilter_augustus.gff3 | tr -d '\n' | grep -oP "protein sequence = [[A-Z\ #]{1,6000}]" | tr -d '#' | tr -d ' ' | tr -d ']' | tr -d '[' | awk ' BEGIN{RS="proteinsequence="}; NR>1 { print ">Rpv1-12_pr"NR-1"\n"$0}' > Rpv1-12-deNovo-extracted_proteins.fa
cd ~/Dropbox/MendelUni_Vinselect/spades/unmapped/Rpv1_12_3/rna
cat Rpv1-12-3_hardF_myFilter_augustus.gff3 | tr -d '\n' | grep -oP "protein sequence = [[A-Z\ #]{1,6000}]" | tr -d '#' | tr -d ' ' | tr -d ']' | tr -d '[' | awk ' BEGIN{RS="proteinsequence="}; NR>1 { print ">Rpv1-12-3_pr"NR-1"\n"$0}' > Rpv1-12-3-deNovo-extracted_proteins.fa
cd ~/Dropbox/MendelUni_Vinselect/spades/unmapped/Suscp/rna
cat Suscpt_hardF_myFilter_augustus.gff3 | tr -d '\n' | grep -oP "protein sequence = [[A-Z\ #]{1,6000}]" | tr -d '#' | tr -d ' ' | tr -d ']' | tr -d '[' | awk ' BEGIN{RS="proteinsequence="}; NR>1 { print ">Suscpt_pr"NR-1"\n"$0}' > Suscpt-deNovo-extracted_proteins.fa


grep -P '\ttranscript\t'  hardF_myFilter_all_cap_augustus.gff3 | awk ' { print $0";Protein=pr"NR} ' >> hardF_myFilter_all_cap_augustus_9columns.gff3

