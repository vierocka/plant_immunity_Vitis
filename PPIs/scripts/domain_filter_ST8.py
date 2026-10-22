"""Strict Pfam filter for Supplementary Table 8 (2026-10-15; reviewer: summary did not follow the stated rules).
hmmscan was run with -E/--domE 0.001 (hmmscan_domain_prediction.sh), so all 623 rows pass the conditional E-value.
Here a hit is kept if full-sequence E < 0.001 AND independent-domain E (i-Evalue) < 0.001; overlapping hits in one
protein are resolved by keeping the highest domain score (no overlap allowed, alignment coordinates).
Sheet 2 gets a 'kept_strict' column (all rows stay visible); sheet 3 is rebuilt as genes per category.
Source (user, 2026-10-15): genome-wide hmmscan of all PN40024 v4 proteins (HMMER 3.3.2, -E/--domE 0.001, May 2023;
reference/Vitis_vinifera/GG/Vitis_vinifera.PN40024.v4.pep.all.hmmscan.{out,tblout}), isoform _P001 - this replaces the
sheet's earlier rows, which came from a later run against a newer Pfam release (e.g. LRR_14).
Usage (repo root): python3 PPIs/scripts/domain_filter_ST8.py <Supplementary_Table_8.xlsx>"""
import sys, re, os, collections, openpyxl
f = sys.argv[1]; wb = openpyxl.load_workbook(f); ws = wb["Domain_architecture"]
GG = os.path.expanduser("~/Dropbox/MendelUni_Vinselect/reference/Vitis_vinifera/GG/Vitis_vinifera.PN40024.v4.pep.all.hmmscan")
# --- rebuild sheet rows from the genome-wide run (gene set and annotation columns taken from the existing sheet)
H0 = [c.value for c in ws[3]]; i0 = {k: i for i, k in enumerate(H0)}
meta, ppi0 = {}, {}
for r in ws.iter_rows(min_row=4, values_only=True):
    if r[i0["query"]]:
        g = r[i0["query"]].split("_")[0]; ppi0.setdefault(g, r[i0["is_PPI_forming_gene"]])
        if r[0] and str(r[0]).startswith("Vitvi"): meta[g] = r[1:4]
genes = sorted(set(ppi0) | {"Vitvi14g04643"})
full = {}
for line in open(GG + ".tblout"):
    if line.startswith("#"): continue
    x = line.split(None, 18)
    if x[2].endswith("_P001") and x[2].split("_")[0] in genes: full[(x[2], x[0])] = (x[1], x[4], x[5], x[6], x[18].strip())
dom = re.compile(r"^\s+(\d+)\s+[!?]\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+\S\S\s+(\d+)\s+(\d+)\s+\S\S\s+(\d+)\s+(\d+)\s+\S\S\s+(\S+)")
new, q, keepq, model, qlen = [], None, False, None, None
for line in open(GG + ".out"):
    if line.startswith("Query:"):
        q = line.split()[1]; qlen = int(line.split("L=")[1].rstrip("]\n")); keepq = q.endswith("_P001") and q.split("_")[0] in genes; model = None; continue
    if not keepq: continue
    if line.startswith(">>"): model = line.split()[1]; continue
    m = dom.match(line)
    if m and model: new.append((q, model, qlen) + m.groups())
nd = collections.Counter((x[0], x[1]) for x in new)
hdr = ["ENSMBL", "UniProt", "Athal_common_name", "Athal_homolog_description", "target", "accession", "query", "qlen",
       "full_seq_E-value", "full_seq_score", "full_seq_bias", "domain_num", "domain_of", "c-Evalue", "i-Evalue", "domain_score",
       "domain_bias", "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to", "acc", "description of target", "is_PPI_forming_gene"]
ws.delete_rows(3, ws.max_row); ws.append(hdr) if ws.max_row < 3 else None
for i, v in enumerate(hdr, 1): ws.cell(row=3, column=i, value=v)
r_i, last = 4, None
for (qq, mdl, ql, n, sc, bias, ce, ie, hf, ht, af, at, ef, et, acc) in new:
    g = qq.split("_")[0]; fa = full[(qq, mdl)]
    head = [g] + list(meta.get(g, (None, None, None))) if g != last else [None] * 4
    vals = head + [mdl, fa[0], qq, ql, float(fa[1]), float(fa[2]), float(fa[3]), int(n), nd[(qq, mdl)], float(ce), float(ie), float(sc),
                   float(bias), int(hf), int(ht), int(af), int(at), int(ef), int(et), float(acc), fa[4], "TRUE" in str(ppi0.get(g))]
    for c, v in enumerate(vals, 1): ws.cell(row=r_i, column=c, value=v)
    r_i += 1; last = g
print("genome-wide rows", len(new), "genes", len({x[0] for x in new}), "missing", sorted(set(genes) - {x[0].split('_')[0] for x in new}))
H = [c.value for c in ws[3]]; ix = {k: i for i, k in enumerate(H)}
rows = [r for r in ws.iter_rows(min_row=4) if r[ix["query"]].value]
gene = lambda r: r[ix["query"]].value.split("_")[0]
ppi = {}
for r in rows: ppi.setdefault(gene(r), bool(r[ix["is_PPI_forming_gene"]].value))
keep = collections.defaultdict(list); kept = set()
cand = [r for r in rows if float(r[ix["full_seq_E-value"]].value) < 1e-3 and float(r[ix["i-Evalue"]].value) < 1e-3]
for r in sorted(cand, key=lambda r: -float(r[ix["domain_score"]].value)):
    a, b = int(r[ix["ali_from"]].value), int(r[ix["ali_to"]].value)
    if all(min(b, d) - max(a, c) + 1 <= 0 for c, d, _ in keep[gene(r)]):
        keep[gene(r)].append((a, b, r[ix["target"]].value)); kept.add(id(r))
col = len(H) + 1
ws.cell(row=3, column=col, value="kept_strict")
for r in rows: ws.cell(row=r[0].row, column=col, value=id(r) in kept)
ws["A1"] = ("Domain_architecture: genome-wide hmmscan (HMMER 3.3.2, Pfam-A; isoform P001) hits for the 43 genes of the structural interaction set "
            "(41 with hits; Vitvi14g04335 and Vitvi14g04643 none). kept_strict = TRUE: full-sequence and independent-domain E-value < 0.001, "
            "best-scoring non-overlapping hit per region. Rows with kept_strict = FALSE are shown for completeness.")
CAT = [("Leucine-rich repeat (LRR)", lambda t: t.startswith("LRR")), ("Protein kinase", lambda t: t in ("Pkinase", "PK_Tyr_Ser-Thr")),
       ("EDS1/PAD4/SAG101 family (EDS1_EP, Lipase_3)", lambda t: t in ("EDS1_EP", "Lipase_3")), ("Other", lambda t: True)]
def cat(t): return next(n for n, fn in CAT if fn(t))
types = {g: sorted(set(x[2] for x in v)) for g, v in keep.items()}
if "domain_summary" in wb.sheetnames: del wb["domain_summary"]
ss = wb.create_sheet("domain_summary", wb.sheetnames.index("Domain_architecture") + 1)
nG, nP = len(types), sum(ppi.get(g, False) for g in set(ppi) | {"Vitvi14g04643"}) + 1  # +1: Vitvi14g04643 (no hit) is PPI-forming
ss.append([f"domain_summary: genes per domain category after the strict filter (43 genes; {nG} with a retained hit; 32 PPI-forming = in a pair passing the strict five-model criterion). Counts are genes, not domain types."])
ss.append([])
ss.append(["Category", "Pfam models kept", f"Genes (of {nG} with hits)", "PPI-forming genes (of 32)", "Genes"])
for n, _ in CAT:
    gs = sorted(g for g, t in types.items() if any(cat(x) == n for x in t))
    ss.append([n, ", ".join(sorted(set(x for g in gs for x in types[g] if cat(x) == n))), len(gs), sum(ppi[g] for g in gs), ", ".join(gs)])
ss.append([]); ss.append(["Per gene (kept Pfam models)"]); ss.append(["gene", "PPI_forming", "kept_models"])
for g in sorted(types): ss.append([g, ppi[g], ", ".join(types[g])])
wb.save(f)
print(len(rows), "rows;", len(kept), "kept;", len(types), "genes;", len(set(x for t in types.values() for x in t)), "Pfam models")
for r in list(ss.iter_rows(min_row=3, max_row=7, values_only=True)): print(r[:4])
