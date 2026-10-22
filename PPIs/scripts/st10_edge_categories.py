#!/usr/bin/env python3
# Fills the empty declared_category(ies) of four candidate pairs in Table S10, sheet all_positive_hits (2026-10-22).
# These pairs come from the literature-based receptor-family screen
# (reference/Vitis_vinifera/GG/colabfold_pairs_BRI1_BAK1_CLV1_BAM1_signaling.fa, built by make_colabfold_pairs.py:
# BRI1/BAK1/CLV1/BAM1 homologs + LRR genes of the SOBIR1 and LRIP1 modules, all combinations).
# Usage (repo root): python3 PPIs/scripts/st10_edge_categories.py
import shutil, openpyxl
SRC = "PPIs/results/Supplementary_Table_8.xlsx"
BMC = "../draft/Sections_by_VK/BMC_PlantBiology/Supplementary_tables/Additional_file_18_TableS10.xlsx"
LABEL = "receptor-family screen"
PAIRS = {frozenset(p) for p in [("Vitvi01g04416", "Vitvi12g01277"), ("Vitvi09g04548", "Vitvi12g01277"),
                                ("Vitvi01g04416", "Vitvi09g04536"), ("Vitvi01g04416", "Vitvi09g04548")]}
wb = openpyxl.load_workbook(SRC); ws = wb["all_positive_hits"]
hdr = [c.value for c in ws[5]]
ia, ib, ic = hdr.index("geneA_gene_id") + 1, hdr.index("geneB_gene_id") + 1, hdr.index("declared_category(ies)") + 1
n = 0
for r in range(6, ws.max_row + 1):
    a, b = ws.cell(r, ia).value, ws.cell(r, ib).value
    if a and b and frozenset((a, b)) in PAIRS and not ws.cell(r, ic).value:
        ws.cell(r, ic).value = LABEL; n += 1
assert n == 20, n            # 4 pairs x 5 models
wb.save(SRC); shutil.copyfile(SRC, BMC)
print(n, "rows labelled; copied to", BMC)
