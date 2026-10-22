#!/usr/bin/env python3
# Rewrites sheet "Module_level_summary" of Additional file 16: Table S9 (2026-10-22): all candidate modules at the three
# correlation thresholds (script 15; no trait-association filter) and STRING edge support (script 16).
# Usage (from 06_PCNWA/): python3 network_robustness/tableS9_module_sheet.py
import csv, statistics, collections, openpyxl
XLSX = "../../draft/Sections_by_VK/BMC_PlantBiology/Supplementary_tables/Additional_file_16_TableS9.xlsx"
D = "network_robustness/15_candidate_modules/"
sizes = collections.defaultdict(list)
for r in csv.DictReader(open(D + "module_sizes_by_cutoff.csv")): sizes[r["r"]].append(int(r["n_partners"]))
sup = list(csv.DictReader(open(D + "string_edge_support_by_threshold.csv")))
pve = [float(r["pve_pc1"]) for r in csv.DictReader(open(D + "module_pc1_variance_explained.csv")) if r["method"].startswith("correlation")]
wb = openpyxl.load_workbook(XLSX)
old = wb["Module_level_summary"]; pos = wb.sheetnames.index("Module_level_summary"); wb.remove(old)
ws = wb.create_sheet("Module_level_summary", pos)
ws.append(["Module = an anchor gene (one of the 9,459 canonical DESeq2 differentially expressed genes) plus every gene of the 26,169-gene "
           "matrix positively correlated with it at r >= threshold (rlog, condition-protected ComBat); anchors with >= 3 partners. "
           "Module sizes count partner genes, excluding the anchor. No module-trait significance filter is applied."])
ws.append([])
ws.append(["Metric", "r >= 0.757 (top 1%)", "r >= 0.8077 (top 0.5%, used throughout)", "r >= 0.8839 (top 0.1%)"])
ks = ["0.757", "0.8077", "0.8839"]
ws.append(["Anchor genes tested"] + [9459] * 3)
ws.append(["Modules (anchors with >= 3 partners)"] + [len(sizes[k]) for k in ks])
ws.append(["Partners per module, median"] + [statistics.median(sizes[k]) for k in ks])
ws.append(["Partners per module, max"] + [max(sizes[k]) for k in ks])
for ev, lab in (("all_evidence", "all evidence"), ("experimental", "experimental only"), ("no_literature", "without text-mining/database")):
    row = {r["r_cutoff"]: r for r in sup if r["evidence"] == ev}
    ws.append([f"Edges with STRING support, {lab} (score >= 400), % (fold vs random pairs)"] +
              [f'{float(row[k]["pct_supported"]):.2f}% ({float(row[k]["fold_enrichment"]):.1f}x)' for k in ks])
ws.append(["Module PC1 variance explained, mean", "", round(statistics.mean(pve), 3), ""])
# sheet NW_density_thresholds (2026-10-22): density of the per-genotype and per-time-point networks at three thresholds
# (combat_protected/scripts/network_density_recompute_r0.8077.R, last section)
if "NW_density_thresholds" in wb.sheetnames: del wb["NW_density_thresholds"]
w3 = wb.create_sheet("NW_density_thresholds")
w3.append(["Co-transcriptional networks of the 9,459 DE genes (rlog, condition-protected ComBat) per genotype (9 libraries) and per time point "
           "(12 libraries) at r >= 0.757 / 0.8077 / 0.8839 (top 1% / 0.5% / 0.1% of all pairwise correlations); positive correlations only."])
w3.append([])
w3.append(["r_cutoff", "network", "type", "n_nodes", "n_edges", "density", "mean_degree", "clustering"])
for r in csv.DictReader(open("combat_protected/results/network_density_three_thresholds_protected_9459panel.csv")):
    w3.append([float(r["r_cutoff"]), r["network"], r["type"], int(r["n_nodes"]), int(r["n_edges"]), round(float(r["density"]), 4),
               round(float(r["mean_degree"]), 1), round(float(r["clustering"]), 3)])
wb.save(XLSX); print("Module_level_summary + NW_density_thresholds written")
