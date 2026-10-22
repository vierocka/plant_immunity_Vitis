#!/usr/bin/env python3
# Adds sheet "leave_one_out" to Additional file 10: Table S6 (2026-10-22): DEG counts per genotype x time point across
# the 36 leave-one-library-out DESeq2 refits (02_Normalization_and_DGEA/additional_tests/scripts/DGEA_leave_one_out_sensitivity.R;
# same design ~ batch + condition and DE criterion FDR < 0.05, |apeglm log2FC| > 1).
# Usage (repo root): python3 05_transcriptional_dynamics/canonical_DESeq2/scripts/tableS6_leave_one_out_sheet.py
import csv, collections, statistics, openpyxl
LOO = "02_Normalization_and_DGEA/additional_tests/results/LOO_sensitivity/LOO_DESeq2_primary_DEG_count_deltas.csv"
XLSX = "../draft/Sections_by_VK/BMC_PlantBiology/Supplementary_tables/Additional_file_10_TableS6.xlsx"
tot, full, own = collections.Counter(), collections.Counter(), collections.defaultdict(set)
rows = list(csv.DictReader(open(LOO)))
first = rows[0]["dropped_sample"]                         # n_full is the same in every refit; take it from one
for r in rows:
    k = (r["genotype"], int(r["timing"]), r["dropped_sample"])
    tot[k] += int(r["n_LOO"])
    if r["dropped_sample"] == first: full[k[:2]] += int(r["n_full"])
    if r["dropped_role"] != "unrelated": own[k[:2]].add(r["dropped_sample"])
genos, times = ["Rpv12", "Rpv12+1", "Rpv12+1+3"], [0, 6, 24]
folds = sorted({k[2] for k in tot}); assert len(folds) == 36
wb = openpyxl.load_workbook(XLSX)
if "leave_one_out" in wb.sheetnames: del wb["leave_one_out"]
ws = wb.create_sheet("leave_one_out")
ws.append(["Leave-one-library-out DESeq2 refits (36 refits, each omitting one of the 36 libraries; design ~ batch + condition; "
           "DE = FDR < 0.05 and |apeglm log2FC| > 1). Counts = up- plus downregulated genes vs the susceptible control."])
ws.append([])
ws.append(["genotype", "time_hpi", "DEGs_full_data", "min_6_refits_without_own_library", "max_6_refits_without_own_library",
           "min_all_36_refits", "max_all_36_refits"])
for t in times:
    for g in genos:
        a = [tot[(g, t, f)] for f in folds]; o = [tot[(g, t, f)] for f in own[(g, t)]]
        ws.append([g, t, full[(g, t)], min(o), max(o), min(a), max(a)])
ws.append([])
ws.append(["time_hpi", "genotype_with_most_DEGs", "refits_where_it_has_most (of 36)", "median_ratio_most_vs_second_ranked_genotype"])
for t in times:
    tops = collections.Counter(); ratios = []
    for f in folds:
        v = sorted(((tot[(g, t, f)], g) for g in genos), reverse=True); tops[v[0][1]] += 1; ratios.append(v[0][0] / v[1][0])
    g, n = tops.most_common(1)[0]
    ws.append([t, g, n, round(statistics.median(ratios), 2)])
wb.save(XLSX); print("sheet leave_one_out written")
