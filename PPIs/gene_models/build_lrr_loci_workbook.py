#!/usr/bin/env python3
# Table S10 sheets 13-16 (gene-model comparison and repeated AF2 predictions). Inputs/outputs outside the repository
# (reference annotations, AF2 run, audit results) are read from VITIS_BASE.
# LRR receptor loci across PN40024 annotations (2026-10-22): one workbook with
#   1 locus_placement     - each PN40024 v4 gene locus aligned to the T2T assembly (minimap2 -x asm5, up to 3 hits)
#   2 longer_model_at_locus - the longer model used for each gene: share of its CDS inside the v4 gene, protein identity (BLASTP)
#   3 construct_boundaries - signal peptide (SignalP 6.0), TM helix (TMHMM 2.0), ectodomain constructs for an AF2 test
# Usage (from the folder containing plant_immunity_Vitis): python3 plant_immunity_Vitis/PPIs/gene_models/build_lrr_loci_workbook.py
import re, csv, subprocess, tempfile, os, openpyxl
from openpyxl.styles import Font, Alignment
B = os.environ.get("VITIS_BASE", os.path.abspath(os.path.join(os.path.dirname(__file__), "../../..")) + "/")   # folder that contains plant_immunity_Vitis/
R = B + "reference/Vitis_vinifera/"
T2T = R + "ASM3070453v1/GCA_030704535.1/"
OUT = B + "gene_model_audit/results/LRR_loci_PN40024_annotations.xlsx"
TMP = tempfile.mkdtemp()
names = {"Vitvi09g04475": "LRIP1", "Vitvi01g04416": "LRIP1-like", "Vitvi09g04548": "LRIP1-like", "Vitvi01g04419": "RLP34",
         "Vitvi12g01277": "BAK1-like", "Vitvi09g01951": "RGI1", "Vitvi09g04536": "GSO2", "Vitvi17g00964": "SOBIR1", "Vitvi19g00098": "MRH1"}
genes = list(names)

def fasta(path, keep=None):
    d, n = {}, None
    for l in open(path):
        if l.startswith(">"): n = l[1:].split()[0]; d[n] = [] if keep is None or keep(n) else None
        elif n and d.get(n) is not None: d[n].append(l.strip())
    return {k: "".join(v).rstrip("*") for k, v in d.items() if v is not None}

best = {r["gene"]: r for r in csv.DictReader(open(B + "gene_model_audit/corrected_models/network_proteins_best_models.tsv"), delimiter="\t")}
longseq = {b.split()[0]: "".join(b.split("\n")[1:]) for b in open(B + "gene_model_audit/corrected_models/network_proteins_best_models.faa").read().split(">")[1:]}
v4 = fasta(R + "GG/Vitis_vinifera.PN40024.v4.pep.all.fa", keep=lambda n: n.split("_")[0] in names and n.endswith("_P001"))
loc = {}
for l in open(R + "GG/Vitis_vinifera.PN40024.v4.56.gff3"):
    f = l.split("\t")
    if len(f) > 8 and f[2] == "gene":
        m = re.search(r"ID=gene:(\w+)", f[8])
        if m and m.group(1) in names: loc[m.group(1)] = (f[0], int(f[3]), int(f[4]), f[6])

# ---- 1 locus placement ----
with open(f"{TMP}/loci.fa", "w") as o:
    for g in genes:
        c, s, e, _ = loc[g]
        o.write(subprocess.run(["samtools", "faidx", R + "GG/reference.fa", f"{c}:{s}-{e}"], capture_output=True, text=True, check=True).stdout.replace(f">{c}:{s}-{e}", f">{g}"))
paf = subprocess.run(["minimap2", "-x", "asm5", "-N", "3", T2T + "GCA_030704535.1_ASM3070453v1_genomic.fna", f"{TMP}/loci.fa"],
                     capture_output=True, text=True, check=True).stdout.splitlines()
place, prim = [], {}
for l in paf:
    f = l.split("\t"); tp = [x for x in f if x.startswith("tp:A:")][0][-1]
    row = [f[0], names[f[0]], f"{loc[f[0]][0]}:{loc[f[0]][1]}-{loc[f[0]][2]} ({loc[f[0]][3]})", int(f[1]), f[5], int(f[7]) + 1, int(f[8]),
           int(f[9]), int(f[10]), round(100 * int(f[9]) / int(f[10]), 2), int(f[11]), "primary" if tp == "P" else "secondary"]
    place.append(row)
    if tp == "P" and f[0] not in prim: prim[f[0]] = (f[5], int(f[7]) + 1, int(f[8]))

# ---- 2 longer model at locus ----
segs = {}
for l in open(T2T + "GCF_030704535.1_ASM3070453v1_translated_cds.faa"):
    if l.startswith(">lcl|NC_0818"):
        m = re.match(r">lcl\|NC_0818(\d\d)\.1_prot_(\S+)", l); lo = re.search(r"\[location=([^\]]+)\]", l); gn = re.search(r"\[gene=([^\]]+)\]", l)
        segs[l[1:].split()[0]] = (f"CP1266{48 + int(m.group(1)) - 5}.1", [(int(a), int(b)) for a, b in re.findall(r"(\d+)\.\.(\d+)", lo.group(1))], gn.group(1) if gn else "")
for l in open(T2T + "genomic.gff"):
    f = l.split("\t")
    if len(f) > 8 and f[2] == "CDS":
        p = re.search(r"protein_id=([^;\n]+)", f[8]).group(1); lt = re.search(r"locus_tag=([^;\n]+)", f[8])
        c = segs.setdefault(p, (f[0], [], lt.group(1) if lt else "")); c[1].append((int(f[3]), int(f[4])))
def blast(q, s):
    open(f"{TMP}/q.faa", "w").write(f">q\n{q}\n"); open(f"{TMP}/s.faa", "w").write(f">s\n{s}\n")
    o = subprocess.run(["blastp", "-query", f"{TMP}/q.faa", "-subject", f"{TMP}/s.faa", "-outfmt", "6 pident length mismatch gaps qstart qend sstart send evalue"],
                       capture_output=True, text=True).stdout.split("\n")[0].split("\t")
    return o
model_rows = []
for g in genes:
    b = best[g]; src = b["source"]; pid = b["model"]; q = v4[g + "_P001"]; sq = longseq[g]
    clean = re.sub(r"^lcl\|NC_[0-9.]+_prot_", "", pid); clean = re.sub(r"_\d+$", "", clean) if src == "RefSeq" else clean
    if src == "v4":
        model_rows.append([g, names[g], len(q), "PN40024 v4 (kept)", clean, len(sq), "", "", "", "", 100.0, len(q), 0, 0, "1-%d" % len(q), "1-%d" % len(q), "same model"])
        continue
    ch, sg, gn = segs[pid]; c, s, e = prim[g]
    tot = sum(y - x + 1 for x, y in sg); ins = sum(max(0, min(y, e) - max(x, s) + 1) for x, y in sg)
    h = blast(q, sq)
    model_rows.append([g, names[g], len(q), {"RefSeq": "NCBI RefSeq (GCF_030704535.1)", "T2T_GenBank": "GenBank (GCA_030704535.1)"}[src], clean, len(sq), gn,
                       f"{min(x for x, _ in sg)}-{max(y for _, y in sg)}", tot, round(100 * ins / tot, 1),
                       float(h[0]), int(h[1]), int(h[2]), int(h[3]), f"{h[4]}-{h[5]}", f"{h[6]}-{h[7]}", ""])

# RNA-seq junction support of each longer model (2026-10-22): its introns vs STAR splice junctions of two libraries of this study
# mapped to the T2T assembly (A11, A12 = Rpv12 0 and 6 hpi; the other SJ files are on the HPC server). SJ.out.tab column 7 =
# uniquely mapped, column 8 = multi-mapped reads; chromosome PNx = CP1266(47+x).1.
import glob, collections
uq, mm = collections.Counter(), collections.Counter()
for fn in sorted(glob.glob(R + "PN40024.T2T/STARmap/*_SJ.out.tab")):
    for l in open(fn):
        f = l.split("\t"); k = (f"CP1266{47 + int(f[0][2:])}.1", int(f[1]), int(f[2])); uq[k] += int(f[6]); mm[k] += int(f[7])
notes = {"Vitvi09g04475": "1,053 aa = 1.16 x Arabidopsis LRIP1 (905 aa); extends 3.3 kb beyond the 5' end of the v4 gene; the locus lies in a "
                          "cluster of LRIP1-like genes and fragments (Vitvi09g04475, 04477, 04478, 04480 within 140 kb), where fused models can arise",
         "Vitvi12g01277": "11 coding exons with lengths as in Arabidopsis BAK1 (AT4G33430, 75% identity; CDS span 3.6 kb), but spread over 34.9 kb; "
                          "the last three introns are seen only in multi-mapped reads, and the last intron (17.5 kb) contains the v4 gene "
                          "Vitvi12g04322 (170 aa, unrelated); the 3' end of the model is uncertain",
         "Vitvi01g04416": "v4 region lies in a duplicated segment (second T2T hit, 84.5% identity, MAPQ 0; Part 1); the closest v4 paralog "
                          "(Vitvi01g04429) shares 92.5% coding-sequence identity; read counts used only uniquely mapped reads"}
for r in model_rows:
    g = r[0]
    if r[3] == "PN40024 v4 (kept)": r[16:16] = ["", "", ""]; continue
    pid = best[g]["model"]; ch, sg, _ = segs[pid]; sg = sorted(sg)
    intr = [(ch, sg[i][1] + 1, sg[i + 1][0] - 1) for i in range(len(sg) - 1)]
    r[16:16] = [len(intr), sum(uq[i] >= 3 for i in intr), sum(uq[i] < 3 and mm[i] >= 3 for i in intr)]
    r[-1] = notes.get(g, r[-1])

# ---- 3 construct boundaries (SignalP 6.0 + TMHMM 2.0; TMHMM reads signal peptides as a first helix) ----
sp = {}
for l in open(B + "gene_model_audit/results/signalp6/prediction_results.txt"):
    if l.startswith("best_"):
        f = l.rstrip("\n").split("\t"); m = re.search(r"CS pos: (\d+)-", f[4] if len(f) > 4 else "")
        sp[f[0].split("_")[1]] = int(m.group(1)) if f[1] == "SP" and m else 0
TM = os.environ.get("TMHMM_DIR", "tmhmm-2.0d") + "/"   # TMHMM 2.0 installation
with open(f"{TMP}/long.faa", "w") as o:
    for g in genes: o.write(f">{g}\n{longseq[g]}\n")
tm = {}
out = subprocess.run([TM + "bin/decodeanhmm.Linux_x86_64", "-f", TM + "lib/TMHMM2.0.options", "-modelfile", TM + "lib/TMHMM2.0.model"],
                     stdin=open(f"{TMP}/long.faa"), capture_output=True, text=True).stdout
cur = None
for l in out.splitlines():
    if l.startswith(">"): cur = l[1:].split()[0]
    elif l.startswith("%pred"): tm[cur] = [(int(a), int(b)) for a, b in re.findall(r"M (\d+) (\d+)", l)]
cons = []
for g in genes:
    L = len(longseq[g]); s = sp.get(g, 0)
    helices = [h for h in tm.get(g, []) if not (s and h[0] <= s)]          # drop the helix that is the signal peptide
    t = helices[0] if helices else None
    ecd = (s + 1, (t[0] - 1) if t else L); ecdtm = (s + 1, min(L, t[1] + 5)) if t else None
    removed = f"{t[1] + 1}-{L}" if t and t[1] < L else ""
    note = {"Vitvi09g04548": "N-terminus incomplete in every annotation", "Vitvi09g01951": "N-terminus and kinase incomplete",
            "Vitvi09g04536": "ectodomain fragment only (no TM in any model)"}.get(g, "")
    cons.append([g, names[g], best[g]["source"], L, f"1-{s}" if s else "none", f"{t[0]}-{t[1]}" if t else "none",
                 f"{ecd[0]}-{ecd[1]}", ecd[1] - ecd[0] + 1, f"{ecdtm[0]}-{ecdtm[1]}" if ecdtm else "", (ecdtm[1] - ecdtm[0] + 1) if ecdtm else "", removed, note])

# ---- write ----
wb = openpyxl.Workbook(); wb.remove(wb.active)
def sheet(title, note, hdr, rows):
    ws = wb.create_sheet(title); ws.append([note]); ws["A1"].font = Font(italic=True); ws.append([]); ws.append(hdr)
    for c in ws[3]: c.font = Font(bold=True)
    for r in rows: ws.append(r)
    for i, h in enumerate(hdr, 1): ws.column_dimensions[openpyxl.utils.get_column_letter(i)].width = max(12, min(40, len(str(h)) + 2))
sheet("README", "LRR receptor and co-receptor loci of a predicted grapevine SOBIR1/LRIP1 network across annotations of the same genome "
      "(Vitis vinifera PN40024). Sheets: 1_locus_placement, 2_longer_model_at_locus, 3_construct_boundaries. Sources: Ensembl PN40024 v4 "
      "(releases 56/57); GenBank and NCBI RefSeq annotations of the telomere-to-telomere assembly GCA/GCF_030704535.1; SignalP 6.0; TMHMM 2.0. "
      "Built by plant_immunity_Vitis/PPIs/gene_models/build_lrr_loci_workbook.py.", ["sheet", "content"],
      [["1_locus_placement", "each v4 gene region aligned to the T2T assembly (minimap2 -x asm5, primary and up to 3 secondary hits)"],
       ["2_longer_model_at_locus", "longer model of the same locus: share of its CDS inside the v4 gene and BLASTP identity with the v4 protein"],
       ["3_construct_boundaries", "signal peptide, first TM helix after it, ectodomain and ectodomain+TM constructs (TM end + 5 residues)"]])
sheet("1_locus_placement", "v4 gene region (genomic, gene start-end) aligned to the T2T assembly with minimap2 -x asm5 -N 3. identity = matches / alignment length.",
      ["gene", "name", "v4 locus", "v4 locus bp", "T2T chromosome", "T2T start", "T2T end", "matches", "alignment length", "identity %", "MAPQ", "hit"], place)
sheet("2_longer_model_at_locus", "Longer protein model at the same locus (T2T coordinates). CDS inside v4 gene = share of the model's coding bases within the "
      "primary T2T placement of the v4 gene. BLASTP = v4 protein (query) vs longer model (subject), best HSP. Junctions: STAR SJ.out.tab "
      "of two libraries of this study mapped to the T2T assembly.",
      ["gene", "name", "v4 aa", "annotation", "model", "model aa", "locus tag / gene", "model span (T2T)", "model CDS bp", "CDS inside v4 gene %",
       "BLASTP % identity", "aligned aa", "mismatches", "gaps", "v4 range", "model range", "introns",
       "introns with >= 3 uniquely mapped junction reads", "introns with multi-mapped junction reads only", "note"], model_rows)
sheet("3_construct_boundaries", "Model used = longest available model (sheet 2). Signal peptide: SignalP 6.0 (eukarya; cleavage after the given residue). "
      "TM: TMHMM 2.0, first helix after the signal peptide (TMHMM also calls signal peptides as helices; those are excluded). "
      "ECD = mature ectodomain (after the signal peptide, before the TM helix); ECD+TM = ECD plus TM helix and 5 juxtamembrane residues.",
      ["gene", "name", "model source", "model aa", "signal peptide", "TM helix", "ECD", "ECD aa", "ECD+TM", "ECD+TM aa", "removed (cytoplasmic)", "note"], cons)
wb.save(OUT); print("written", OUT)
for r in cons: print(r)

# ---- Table S10 sheet "gene_models_loci" (2026-10-22): sheets 1 and 2 above, for the BMC manuscript ----
import shutil
REPO = B + "plant_immunity_Vitis/PPIs/results/Supplementary_Table_8.xlsx"
BMC = B + "draft/Sections_by_VK/BMC_PlantBiology/Supplementary_tables/Additional_file_18_TableS10.xlsx"
w8 = openpyxl.load_workbook(REPO)
if "gene_models_loci" in w8.sheetnames: del w8["gene_models_loci"]
s = w8.create_sheet("gene_models_loci")
s.append(["Genome sequence and protein models of the Figure 4A loci in two PN40024 assemblies. Part 1: each PN40024 v4 gene region "
          "aligned to the telomere-to-telomere (T2T) assembly (minimap2 -x asm5); identity = matches / alignment length. Part 2: the longer "
          "model of the same locus used for the repeated predictions (sheet gene_models_network): share of its coding sequence inside the "
          "v4 gene, BLASTP comparison with the v4 protein (best alignment), and support of its introns by splice junctions in two RNA-seq libraries of "
          "this study (Rpv12, 0 and 6 hpi) mapped to the T2T assembly with STAR (>= 3 uniquely mapped reads; or seen only in multi-mapped reads, "
          "i.e. the sequence also occurs elsewhere in the genome). Both annotations are computational predictions; a longer model is not "
          "necessarily the correct one."]); s["A1"].font = Font(italic=True)
s.append([]); s.append(["Part 1. Locus placement"]); s.cell(s.max_row, 1).font = Font(bold=True)
s.append(["gene", "name", "v4 locus", "v4 locus bp", "T2T chromosome", "T2T start", "T2T end", "matches", "alignment length", "identity %", "MAPQ", "hit"])
for c in s[s.max_row]: c.font = Font(bold=True)
for r in place: s.append(r)
s.append([]); s.append(["Part 2. Longer model at the same locus"]); s.cell(s.max_row, 1).font = Font(bold=True)
s.append(["gene", "name", "v4 aa", "annotation", "model", "model aa", "locus tag / gene", "model span (T2T)", "model CDS bp", "CDS inside v4 gene %",
          "BLASTP % identity", "aligned aa", "mismatches", "gaps", "v4 range", "model range", "introns",
       "introns with >= 3 uniquely mapped junction reads", "introns with multi-mapped junction reads only", "note"])
for c in s[s.max_row]: c.font = Font(bold=True)
for r in model_rows: s.append(r)
ORDER = ["all_positive_hits", "Domain_architecture", "domain_summary", "roadmap_all_PPI_tested_pairs", "within_module_all_models",
         "cross_module_all_models", "negative_controls_all_models", "positive_control_EDS1_SAG101_PA", "promiscuity_stickiness_tests_al",
         "failures_zero_contacts", "failures_chainA_empty", "failures_chainB_empty", "gene_models_network", "gene_models_43",
         "AF2_rerun_gene_models", "gene_models_loci"]   # sheet order of the Table S10 caption
w8._sheets.sort(key=lambda x: ORDER.index(x.title) if x.title in ORDER else 99)
w8.save(REPO); shutil.copyfile(REPO, BMC); print("Table S10 sheet gene_models_loci written")
