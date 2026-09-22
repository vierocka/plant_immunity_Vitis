cd ~/Dropbox/MendelUni_Vinselect/spades/unmapped
hmmpress ~/anaconda3/bin/Pfam-A.hmm
hmmscan  -o hardF_myFilter_all_cap_augustus_proteins.hmmscan.out --tblout hardF_myFilter_all_cap_augustus_proteins.hmmscan.tblout -E 0.001 --domE 0.001 --max --cpu 8 ~/anaconda3/bin/Pfam-A.hmm hardF_myFilter_all_cap_augustus_proteins.fa

# the pipeline has been interrupted; continue:
hmmscan  -o hardF_myFilter_all_cap_augustus_proteins_p2.hmmscan.out --tblout hardF_myFilter_all_cap_augustus_proteins_p2.hmmscan.tblout -E 0.001 --domE 0.001 --max --cpu 8 ~/anaconda3/bin/Pfam-A.hmm hardF_myFilter_all_cap_augustus_proteins_p2.fa

rm hardF_myFilter_all_cap_augustus_proteins_p2.fa
cat  hardF_myFilter_all_cap_augustus_proteins.hmmscan.tblout hardF_myFilter_all_cap_augustus_proteins_p2.hmmscan.tblout > hardF_myFilter_all_cap_augustus_7658proteins.hmmscan.tblout
rm hardF_myFilter_all_cap_augustus_proteins.hmmscan.tblout
rm hardF_myFilter_all_cap_augustus_proteins_p2.hmmscan.tblout
cat hardF_myFilter_all_cap_augustus_proteins.hmmscan.out hardF_myFilter_all_cap_augustus_proteins_p2.hmmscan.out > hardF_myFilter_all_cap_augustus_7658proteins.hmmscan.out
rm hardF_myFilter_all_cap_augustus_proteins.hmmscan.out
rm hardF_myFilter_all_cap_augustus_proteins_p2.hmmscan.out

# I did not run this one:
# hmmscan  -o All_deNovo-extracted_proteins_perGenotype.hmmscan.out --tblout All_deNovo-extracted_proteins_perGenotype.hmmscan.tblout -E 0.001 --domE 0.001 --cpu 8 ~/anaconda3/bin/Pfam-A.hmm All_deNovo-extracted_proteins.fa

# ~/Dropbox/MendelUni_Vinselect/reference/Vitis_vinifera/GG/Vitis_vinifera.PN40024.v4.pep.all.fa
hmmscan  -o ~/Dropbox/MendelUni_Vinselect/reference/Vitis_vinifera/GG/Vitis_vinifera.PN40024.v4.pep.all.hmmscan.out --tblout ~/Dropbox/MendelUni_Vinselect/reference/Vitis_vinifera/GG/Vitis_vinifera.PN40024.v4.pep.all.hmmscan.tblout -E 0.001 --domE 0.001 --cpu 8 ~/anaconda3/bin/Pfam-A.hmm ~/Dropbox/MendelUni_Vinselect/reference/Vitis_vinifera/GG/Vitis_vinifera.PN40024.v4.pep.all.fa
