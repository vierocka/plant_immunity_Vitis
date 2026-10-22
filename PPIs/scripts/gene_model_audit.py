#!/usr/bin/env python3
# Gene-model audit of the 43-gene structural interaction set (2026-10-19).
# PN40024 v4 models are computational predictions; this flags models whose length disagrees with the closest
# Arabidopsis homolog (> 25 aa of the homolog missing at either end of the best HSP) and checks whether the
# independent NCBI RefSeq model of the same locus (GCF_030704535.1, T2T PN40024) is longer at that end.
# Usage: python3 PPIs/scripts/gene_model_audit.py > PPIs/results/gene_model_audit_43genes.tsv
import subprocess, tempfile, os, re, sys, openpyxl
R = os.path.expanduser("~/Dropbox/MendelUni_Vinselect/reference/Vitis_vinifera/")
V4 = R + "GG/Vitis_vinifera.PN40024.v4.pep.all.fa"
TAIR_DB = R + "TAIR10/TAIR10.proteins"
REFSEQ_DB = R + "ASM3070453v1/GCA_030704535.1/GCF_030704535.1_ASM3070453v1_translated_cds.db"
REFSEQ_FA = R + "ASM3070453v1/GCA_030704535.1/GCF_030704535.1_ASM3070453v1_translated_cds.faa"
ST10 = os.path.join(os.path.dirname(__file__), "../results/Supplementary_Table_8.xlsx")
MISS = 25

def fasta(path, keep=None):
    seqs, name = {}, None
    for l in open(path):
        if l.startswith(">"):
            name = l[1:].split()[0]; seqs[name] = [] if keep is None or keep(name) else None
        elif name and seqs.get(name) is not None: seqs[name].append(l.strip())
    return {k: "".join(v) for k, v in seqs.items() if v is not None}

rows = list(openpyxl.load_workbook(ST10, read_only=True)["Domain_architecture"].iter_rows(values_only=True))
genes, names = [], {}
for r in rows[3:]:
    if r[0] and re.match(r"Vitvi\d", str(r[0])) and r[0] not in names: genes.append(r[0]); names[r[0]] = r[2] or ""
for g in re.findall(r"(Vitvi\w+) and (Vitvi\w+) none", rows[0][0])[0]:   # the two genes without Pfam hits, named in the sheet header
    genes.append(g); names[g] = ""
assert len(genes) == 43, len(genes)
v4 = fasta(V4, keep=lambda n: n.split("_")[0] in names and n.endswith("_P001"))
assert len(v4) == 43, len(v4)

def blast(seqs, db, extra):
    with tempfile.NamedTemporaryFile("w", suffix=".faa", delete=False) as f:
        for k, s in seqs.items(): f.write(f">{k}\n{s}\n")
    out = subprocess.run(["blastp", "-query", f.name, "-db", db, "-evalue", "1e-6", *extra, "-outfmt",
                          "6 qseqid sseqid pident length qstart qend qlen sstart send slen evalue bitscore"],
                         capture_output=True, text=True, check=True).stdout
    os.unlink(f.name)
    return [dict(zip("q s pid len qs qe ql ss se sl e bit".split(), l.split("\t"))) for l in out.strip().splitlines()]

def best_tair(seqs):
    best = {}
    for h in blast(seqs, TAIR_DB, ["-max_target_seqs", "5"]):
        k = (float(h["e"]), -float(h["bit"]))
        if h["q"] not in best or k < best[h["q"]][0]: best[h["q"]] = (k, h)
    res = {}
    for q, (_, h) in best.items():
        qs, qe, ql, ss, se, sl = (int(h[x]) for x in ("qs", "qe", "ql", "ss", "se", "sl"))
        res[q] = dict(tair=h["s"], pid=float(h["pid"]), ql=ql, sl=sl, cov=f"{ss}-{se}",
                      n_missing=(ss - 1) - (qs - 1), c_missing=(sl - se) - (ql - qe))
    return res

t4 = best_tair(v4)
# same-locus RefSeq model: >= 98% identity over >= 80% of the v4 model; keep the longest such protein
same = {}
for h in blast(v4, REFSEQ_DB, ["-max_target_seqs", "20"]):
    if float(h["pid"]) >= 98 and int(h["len"]) >= 0.8 * int(h["ql"]):
        if h["q"] not in same or int(h["sl"]) > same[h["q"]][1]: same[h["q"]] = (h["s"], int(h["sl"]))
rs_ids = {v[0] for v in same.values()}
rs = fasta(REFSEQ_FA, keep=lambda n: n in rs_ids)
trs = best_tair(rs)

cols = ["gene", "name", "v4_len", "best_Athal_hit", "Athal_len", "pident", "Athal_aligned", "N_missing", "C_missing", "flag",
        "RefSeq_model", "RefSeq_len", "RefSeq_N_missing", "RefSeq_C_missing", "RefSeq_resolves"]
print("\t".join(cols))
for g in genes:
    q = g + "_P001"; a = t4.get(q)
    if a is None: print("\t".join([g, names[g], str(len(v4[q]))] + ["NA"] * 12)); continue
    flag = ";".join(x for x, v in (("N", a["n_missing"]), ("C", a["c_missing"])) if v > MISS) or "ok"
    r, rl = same.get(q, ("NA", "NA")); b = trs.get(r)
    resolves = "NA" if b is None or flag == "ok" else str(all(
        b[k] <= MISS for k, e in (("n_missing", "N"), ("c_missing", "C")) if e in flag))
    print("\t".join(map(str, [g, names[g], a["ql"], a["tair"], a["sl"], a["pid"], a["cov"], a["n_missing"], a["c_missing"], flag,
                              r, rl, b["n_missing"] if b else "NA", b["c_missing"] if b else "NA", resolves])))
