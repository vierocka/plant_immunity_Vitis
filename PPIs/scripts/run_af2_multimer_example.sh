#!/bin/bash
# AIM: run one AlphaFold2-Multimer (ColabFold) pairwise interaction
# prediction, end to end. Generic template -- swap in any two protein
# FASTA files for a different pair.
#
# Worked example here: SOBIR1 (Vitvi17g00964) x LRIP1 (Vitvi09g04475),
# the manuscript's core positive finding -- sequences in ../full_example/.
#
# Every real screen (within-module, cross-module, positive/negative
# controls, stickiness panel) in this project is the same command run
# once per candidate pair; only the input FASTA and output folder change.

set -euo pipefail
source ~/localcolabfold/conda/etc/profile.d/conda.sh
conda activate ~/localcolabfold/colabfold-conda

SEQ_A="../full_example/SOBIR1_Vitvi17g00964.fa"
SEQ_B="../full_example/LRIP1_Vitvi09g04475.fa"
PAIR_NAME="Vitvi17g00964_x_Vitvi09g04475"
OUT_DIR="../full_example/${PAIR_NAME}_output"

mkdir -p "$OUT_DIR"

# ColabFold multimer input: one FASTA with both chains, ':'-joined.
PAIR_FASTA="${OUT_DIR}/${PAIR_NAME}.fasta"
{
  echo ">${PAIR_NAME}"
  seqA=$(grep -v "^>" "$SEQ_A" | tr -d '\n')
  seqB=$(grep -v "^>" "$SEQ_B" | tr -d '\n')
  echo "${seqA}:${seqB}"
} > "$PAIR_FASTA"

# 5 independent models, 3 recycles each -- the manuscript's standard
# per-pair setting throughout (see ../README.md for the criteria this feeds).
colabfold_batch \
  --model-type alphafold2_multimer_v3 \
  --num-models 5 \
  --num-recycle 3 \
  "$PAIR_FASTA" "$OUT_DIR"

echo "Done. Per-model JSON (ipTM/pTM/pLDDT/PAE) and PDB structures in: $OUT_DIR"
echo "Next: extract ipTM/contacts/PAE from the JSON files and apply the"
echo "5-model strict criterion described in ../README.md."
