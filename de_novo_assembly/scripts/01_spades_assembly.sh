#cd ~/Dropbox/MendelUni_Vinselect/spades
#mkdir unmapped
#mkdir unmapped/Rpv1
#mkdir unmapped/Rpv1/rna
#mkdir unmapped/Rpv1/meta

#mkdir unmapped/Rpv1_12
#mkdir unmapped/Rpv1_12/rna
#mkdir unmapped/Rpv1_12/meta

#mkdir unmapped/Rpv1_12_3
#mkdir unmapped/Rpv1_12_3/rna
#mkdir unmapped/Rpv1_12_3/meta

#mkdir unmapped/Suscp
#mkdir unmapped/Suscp/rna
#mkdir unmapped/Suscp/meta

#cd ~/Dropbox/MendelUni_Vinselect/
#cat STARmap1/?1?_Unmapped.out.mate1 > spades/unmapped/Rpv1_1.fq
#cat STARmap1/?1?_Unmapped.out.mate2 > spades/unmapped/Rpv1_2.fq

#cat STARmap1/?2?_Unmapped.out.mate1 > spades/unmapped/Rpv1_12_1.fq
#cat STARmap1/?2?_Unmapped.out.mate2 > spades/unmapped/Rpv1_12_2.fq

#cat STARmap1/?3?_Unmapped.out.mate1 > spades/unmapped/Rpv1_12_3_1.fq
#cat STARmap1/?3?_Unmapped.out.mate2 > spades/unmapped/Rpv1_12_3_2.fq

#cat STARmap1/?4?_Unmapped.out.mate1 > spades/unmapped/Suscp_1.fq
#cat STARmap1/?4?_Unmapped.out.mate2 > spades/unmapped/Suscp_2.fq

cd ~/Dropbox/MendelUni_Vinselect/spades
#spades -o unmapped/Rpv1/rna --only-assembler --rna -t 8 -m 24 -1 unmapped/Rpv1_1.fq -2 unmapped/Rpv1_2.fq
#spades -o unmapped/Rpv1/meta --only-assembler --meta -t 8 -m 24 -1 unmapped/Rpv1_1.fq -2 unmapped/Rpv1_2.fq

#spades -o unmapped/Rpv1_12/rna --only-assembler --rna -t 8 -m 24 -1 unmapped/Rpv1_12_1.fq -2 unmapped/Rpv1_12_2.fq
#spades -o unmapped/Rpv1_12/meta --only-assembler --meta -t 8 -m 24 -1 unmapped/Rpv1_12_1.fq -2 unmapped/Rpv1_12_2.fq

#spades -o unmapped/Rpv1_12_3/rna --only-assembler --rna -t 8 -m 24 -1 unmapped/Rpv1_12_3_1.fq -2 unmapped/Rpv1_12_3_2.fq
#spades -o unmapped/Rpv1_12_3/meta --only-assembler --meta -t 8 -m 24 -1 unmapped/Rpv1_12_3_1.fq -2 unmapped/Rpv1_12_3_2.fq

spades -o unmapped/Suscp/rna --only-assembler --rna -t 8 -m 24 -1 unmapped/Suscp_1.fq -2 unmapped/Suscp_2.fq
spades -o unmapped/Suscp/meta --only-assembler --meta -t 8 -m 24 -1 unmapped/Suscp_1.fq -2 unmapped/Suscp_2.fq

exit 0
