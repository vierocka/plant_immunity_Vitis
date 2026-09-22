#!/bin/bash
# AIM: Pfam-A domain annotation (HMMER hmmscan) for every gene appearing in
# any high-confidence AF2-Multimer interaction -- the domain architecture
# reported in Supplementary_Table_8.xlsx (Domain_architecture, domain_summary).
#
# Input: hmmscan_AF2PPIs_seqs.fa (protein FASTA of all genes in the
# high-confidence interaction set, one entry per gene).

hmmpress ~/Desktop/SW/hmmdb/Pfam-A.hmm

hmmscan \
  --tblout hmmscan_AF2PPIs_seqs_tbl.out \
  --domtblout hmmscan_AF2PPIs_seqs_domtbl.out \
  --noali -E 0.001 --domE 0.001 --seed 93092 --cpu 8 \
  ~/Desktop/SW/hmmdb/Pfam-A.hmm \
  hmmscan_AF2PPIs_seqs.fa
