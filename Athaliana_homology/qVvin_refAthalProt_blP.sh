 # TAIR10 dataset downloaded 11. April 2025 
 # wget https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCF_000001735.3/download?include_annotation_type=GENOME_FASTA
makeblastdb -dbtype prot -out TAIR10.proteins -in GCF_000001735.3_TAIR10_protein.faa
blastp -query GG/Vitis_vinifera.PN40024.v4.pep.all.fa -db TAIR10/TAIR10.proteins -outfmt 6 -max_target_seqs 10 -num_threads 4 -evalue 0.000001 -out qPN40024.v4.pep_refTAIR10_blp.out
