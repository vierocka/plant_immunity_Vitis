cd ~/Dropbox/MendelUni_Vinselect/spades/unmapped
makeblastdb -dbtype prot -in hardF_myFilter_all_cap_augustus_proteins.fa -out hardF_myFilter_all_cap_augustus_proteins.db
blastp -db hardF_myFilter_all_cap_augustus_proteins.db -query hardF_myFilter_all_cap_augustus_proteins.fa -out hardF_myFilter_all_cap_augustus_proteins.bdBlP.out -evalue 0.000001 -num_threads 4 -outfmt 6
mcl hardF_myFilter_all_cap_augustus_proteins.bdBlP.out --abc -I 2 -o hardF_myFilter_all_cap_augustus_proteins.bdBlP.I2mcl
mcl hardF_myFilter_all_cap_augustus_proteins.bdBlP.out --abc -I 1 -o hardF_myFilter_all_cap_augustus_proteins.bdBlP.I1mcl
mcl hardF_myFilter_all_cap_augustus_proteins.bdBlP.out --abc -I 1.6 -o hardF_myFilter_all_cap_augustus_proteins.bdBlP.I1.6mcl
mcl hardF_myFilter_all_cap_augustus_proteins.bdBlP.out --abc -I 1.3 -o hardF_myFilter_all_cap_augustus_proteins.bdBlP.I1.3mcl
mcl hardF_myFilter_all_cap_augustus_proteins.bdBlP.out --abc -I 2.5 -o hardF_myFilter_all_cap_augustus_proteins.bdBlP.I2.5mcl
