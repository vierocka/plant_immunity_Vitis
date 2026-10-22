#!/usr/bin/env python3
# Table S10 sheets 13-16 (gene-model comparison and repeated AF2 predictions). Inputs/outputs outside the repository
# (reference annotations, AF2 run, audit results) are read from VITIS_BASE.
# Adds sheet "gene_models_network" to Table S10 (repo PPIs/results/Supplementary_Table_8.xlsx and the BMC copy
# Additional_file_18_TableS10.xlsx), 2026-10-22: PN40024 v4 model vs best available complete model for the nine
# Figure 4A genes (cross-annotation check, SignalP 6.0). Usage (from the folder containing plant_immunity_Vitis):
#   python3 plant_immunity_Vitis/PPIs/gene_models/st10_gene_model_sheet.py
import os, csv, shutil, re, collections, openpyxl
from openpyxl.styles import Font
B = os.environ.get("VITIS_BASE", os.path.abspath(os.path.join(os.path.dirname(__file__), "../../..")) + "/")   # folder that contains plant_immunity_Vitis/
REPO = B + "plant_immunity_Vitis/PPIs/results/Supplementary_Table_8.xlsx"
BMC = B + "draft/Sections_by_VK/BMC_PlantBiology/Supplementary_tables/Additional_file_18_TableS10.xlsx"
best = {r["gene"]: r for r in csv.DictReader(open(B + "gene_model_audit/corrected_models/network_proteins_best_models.tsv"), delimiter="\t")}
aud = {r["gene"]: r for r in csv.DictReader(open(B + "gene_model_audit/results/gene_model_cross_annotation_43genes.tsv"), delimiter="\t")}
stage1 = {r["gene"]: r for r in csv.DictReader(open(B + "plant_immunity_Vitis/PPIs/results/gene_model_audit_43genes.tsv"), delimiter="\t")}
sp = {}
for l in open(B + "gene_model_audit/results/signalp6/prediction_results.txt"):
    if l.startswith("#"): continue
    f = l.rstrip("\n").split("\t"); cs = re.search(r"CS pos: (\d+)-", f[4]) if len(f) > 4 else None
    sp[f[0]] = f"yes, cleavage after {cs.group(1)}" if f[1] == "SP" else "no"
names = {"Vitvi17g00964": "SOBIR1", "Vitvi09g04475": "LRIP1", "Vitvi01g04416": "LRIP1-like", "Vitvi09g04548": "LRIP1-like",
         "Vitvi01g04419": "RLP34", "Vitvi09g01951": "RGI1", "Vitvi09g04536": "GSO2", "Vitvi12g01277": "BAK1-like", "Vitvi19g00098": "MRH1"}
order = ["Vitvi09g04475", "Vitvi01g04416", "Vitvi01g04419", "Vitvi12g01277", "Vitvi09g04548", "Vitvi09g01951", "Vitvi09g04536",
         "Vitvi17g00964", "Vitvi19g00098"]
srcname = {"v4": "PN40024 v4 (unchanged)", "RefSeq": "NCBI RefSeq, T2T assembly GCF_030704535.1",
           "T2T_GenBank": "GenBank, T2T assembly GCA_030704535.1"}
for path in (REPO,):
    wb = openpyxl.load_workbook(path)
    if "gene_models_network" in wb.sheetnames: del wb["gene_models_network"]
    ws = wb.create_sheet("gene_models_network")
    ws.append(["Gene models of the nine Figure 4A genes in four annotations of the PN40024 genome. Missing N / C = residues of the closest "
               "Arabidopsis homolog not covered at the N or C end of the best BLASTP alignment (negative = model longer). Complete = <= 25 "
               "residues missing at both ends and length <= 1.25 x homolog (longer, fused models are not counted). Signal peptide: SignalP 6.0. "
               "The genome sequence at these loci is nearly identical between the v4 and telomere-to-telomere assemblies, and the annotations "
               "differ mainly in exon structure (sheet gene_models_loci). All models are computational predictions; a complete model is not "
               "thereby validated."])
    ws["A1"].font = Font(italic=True)
    ws.append([])
    ws.append(["Annotations compared", "assembly", "version / accession", "obtained"])
    for row in (["PN40024 v4 (Ensembl Plants; used in this study)", "PN40024 v4", "releases 56 and 57 (identical at these loci)", "2023"],
                ["GenBank annotation of the telomere-to-telomere assembly", "PN40024 T2T", "GCA_030704535.1", "Oct 2023"],
                ["NCBI RefSeq annotation of the telomere-to-telomere assembly", "PN40024 T2T", "GCF_030704535.1", "Oct 2023"],
                ["Phytozome", "PN40024 12X", "Vvinifera v2.1", "Dec 2024"]):
        ws.append(row)
    for c in ws[3]: c.font = Font(bold=True)
    ws.append([])
    hdr = ["gene", "name", "Arabidopsis_homolog", "homolog_aa", "v4_aa", "v4_missing_N", "v4_missing_C", "v4_signal_peptide",
           "T2T_GenBank_model", "T2T_GenBank_aa", "T2T_GenBank_missing_N", "T2T_GenBank_missing_C",
           "RefSeq_model", "RefSeq_aa", "RefSeq_missing_N", "RefSeq_missing_C",
           "Phytozome_model", "Phytozome_aa", "Phytozome_missing_N", "Phytozome_missing_C",
           "model_used_for_AF2_rerun", "used_signal_peptide", "status"]
    ws.append(hdr); hr = ws.max_row
    for c in ws[hr]: c.font = Font(bold=True)
    def blk(a, k):
        m = a[f"{k}_model"]
        if m == "NA": return ["no model of plausible size at locus", None, None, None]
        m = re.sub(r"^lcl\|NC_\d+\.\d_prot_", "", m); m = re.sub(r"_\d+$", "", m) if k == "RefSeq" else m
        return [m, int(a[f"{k}_len"]), int(a[f"{k}_N_missing"]), int(a[f"{k}_C_missing"])]
    for g in order:
        b, a, s1 = best[g], aud[g], stage1[g]
        src = b["source"]; pid = re.sub(r"^lcl\|NC_\d+\.\d_prot_", "", b["model"]); pid = re.sub(r"_\d+$", "", pid) if src == "RefSeq" else pid
        status = "complete (v4 unchanged)" if src == "v4" else ("complete" if b["status"] == "full" else "longer than v4, still partial")
        if g in ("Vitvi09g04475", "Vitvi12g01277"): status += "; uncertain, see sheet gene_models_loci, Part 2"
        ws.append([g, names[g], s1["best_Athal_hit"], int(s1["Athal_len"]), int(a["v4_len"]), int(a["v4_N_missing"]), int(a["v4_C_missing"]),
                   sp["v4_" + g]] + blk(a, "T2T_GenBank") + blk(a, "RefSeq") + blk(a, "Phytozome_v2.1") +
                  [f"{srcname[src]}: {pid}", sp[f"best_{g}_{src}"], status])

    # ---- sheet gene_models_43: all genes of the structural interaction set ----
    if "gene_models_43" in wb.sheetnames: del wb["gene_models_43"]
    w2 = wb.create_sheet("gene_models_43")
    w2.append(["Gene models of the 43 genes of the structural interaction set in four annotations of the PN40024 genome (as in sheet "
               "gene_models_network). complete = <= 25 residues of the closest Arabidopsis homolog missing at both ends and length <= 1.25 x "
               "homolog; '-' = no model of plausible size at the locus. Vitvi14g04335 has no Arabidopsis homolog hit."])
    w2["A1"].font = Font(italic=True); w2.append([])
    w2.append(["gene", "name", "Figure_4", "Arabidopsis_homolog", "homolog_aa", "v4_aa", "v4_complete",
               "T2T_GenBank_aa", "T2T_GenBank_complete", "RefSeq_aa", "RefSeq_complete", "Phytozome_aa", "Phytozome_complete", "complete_in_any_annotation"])
    for c in w2[3]: c.font = Font(bold=True)
    figA = set(order); figC = {"Vitvi17g04216", "Vitvi05g00062", "Vitvi14g03030", "Vitvi14g03031", "Vitvi14g03032", "Vitvi14g04642",
                               "Vitvi14g03033", "Vitvi07g01908"}
    def full(r, k, hl):   # 2026-10-22: length rule (<= 1.25 x homolog) now applied to v4 too (was applied only to the other annotations)
        n_, c_ = r[f"{k}_N_missing"], r[f"{k}_C_missing"]
        return "-" if n_ in ("NA", "") else ("yes" if int(n_) <= 25 and int(c_) <= 25 and int(r[f"{k}_len"]) <= 1.25 * hl else "no")
    tot = collections.Counter(); n43 = 0
    for g, r in aud.items():
        s1 = stage1[g]; hom = s1["best_Athal_hit"] != "NA"
        f = {k: full(r, k, int(s1["Athal_len"]) if hom else 0) for k in ("v4", "T2T_GenBank", "RefSeq", "Phytozome_v2.1")}
        anyf = ("yes" if "yes" in f.values() else "no") if hom else "NA"
        if hom:
            n43 += 1
            for k, v in f.items(): tot[k] += v == "yes"
            tot["any"] += anyf == "yes"
        L = lambda k: (int(r[f"{k}_len"]) if r.get(f"{k}_len", "NA") not in ("NA", "") else "-")
        w2.append([g, r["name"] or "", "A" if g in figA else ("C" if g in figC else ""), s1["best_Athal_hit"],
                   int(s1["Athal_len"]) if hom else None, int(r["v4_len"]), f["v4"] if hom else "NA",
                   L("T2T_GenBank"), f["T2T_GenBank"] if hom else "NA", L("RefSeq"), f["RefSeq"] if hom else "NA",
                   L("Phytozome_v2.1"), f["Phytozome_v2.1"] if hom else "NA", anyf])
    w2.append([])
    w2.append([f"Complete models among the {n43} genes with an Arabidopsis homolog", "", "", "", "", "", tot["v4"], "", tot["T2T_GenBank"],
               "", tot["RefSeq"], "", tot["Phytozome_v2.1"], tot["any"]])
    w2.cell(w2.max_row, 1).font = Font(bold=True)
    ORDER = ["all_positive_hits", "Domain_architecture", "domain_summary", "roadmap_all_PPI_tested_pairs", "within_module_all_models",
             "cross_module_all_models", "negative_controls_all_models", "positive_control_EDS1_SAG101_PA", "promiscuity_stickiness_tests_al",
             "failures_zero_contacts", "failures_chainA_empty", "failures_chainB_empty", "gene_models_network", "gene_models_43",
             "AF2_rerun_gene_models", "gene_models_loci"]   # sheet order of the Table S10 caption
    wb._sheets.sort(key=lambda x: ORDER.index(x.title) if x.title in ORDER else 99)
    wb.save(path)
shutil.copyfile(REPO, BMC)
print("sheet gene_models_network written to Table S10 (repo + BMC copy)")
