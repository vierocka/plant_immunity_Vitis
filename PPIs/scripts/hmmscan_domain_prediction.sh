#!/bin/bash
# AIM: Pfam-A domain annotation (HMMER hmmscan) for every gene appearing in
# any high-confidence AF2-Multimer interaction -- the domain architecture
# reported in Supplementary_Table_8.xlsx (Domain_architecture, domain_summary).
#
# Input: hmmscan_AF2PPIs_seqs.fa (protein FASTA of all genes in the
# high-confidence interaction set, one entry per gene).
#
# NOTE 2026-10-15: this run used a newer Pfam-A (Jan 2026; e.g. LRR_14). ST8 sheets 2-3 are now built from the
# genome-wide HMMER 3.3.2 run (reference/Vitis_vinifera/GG/Vitis_vinifera.PN40024.v4.pep.all.hmmscan.*, May 2023)
# by domain_filter_ST8.py, matching the Methods.

hmmpress ~/Desktop/SW/hmmdb/Pfam-A.hmm

hmmscan \
  --tblout hmmscan_AF2PPIs_seqs_tbl.out \
  --domtblout hmmscan_AF2PPIs_seqs_domtbl.out \
  --noali -E 0.001 --domE 0.001 --seed 93092 --cpu 8 \
  ~/Desktop/SW/hmmdb/Pfam-A.hmm \
  hmmscan_AF2PPIs_seqs.fa
