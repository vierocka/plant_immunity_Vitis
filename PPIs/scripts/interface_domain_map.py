"""2026-10-19: domain of every predicted interface residue in SOBIR1-LRIP1 AlphaFold2-Multimer models (reviewer: an interface on
the SOBIR1 kinase would be topologically impossible for an extracellular LRR-only partner). Contacts as in the Methods:
inter-chain C-alpha-C-alpha < 8 A between residues with pLDDT (B-factor column) >= 50.
Domain map of SOBIR1 Vitvi17g00964_P001 (632 aa; genome-wide hmmscan + Kyte-Doolittle hydropathy): signal-peptide-like 1-27,
extracellular incl. LRRs 28-276 (Pfam LRR hits 100-214), transmembrane helix 277-313 (GxxxGxxxG motifs 283-291, 287-295),
juxtamembrane 314-339, kinase 340-628 (Pfam Pkinase). LRIP1 Vitvi09g04475_P001 (410 aa) = LRRs 5-377, no signal peptide or TM in this
gene model (aligns to Arabidopsis LRIP1 NP_181039.1 residues 363-770; the Arabidopsis protein has an N-terminal signal peptide).
Usage: python3 PPIs/scripts/interface_domain_map.py SOBIR1_CHAIN LRIP1_CHAIN model1.pdb [model2.pdb ...] > out.tsv"""
import sys, math, collections
DOM = [(1, 27, "signal_peptide"), (28, 276, "extracellular_LRR"), (277, 313, "transmembrane"), (314, 339, "juxtamembrane"), (340, 628, "kinase"), (629, 632, "C_terminus")]
dom = lambda i: next((n for a, b, n in DOM if a <= i <= b), "outside_map")
def ca(path):
    res = collections.defaultdict(dict)
    for l in open(path):
        if l.startswith("ATOM") and l[12:16].strip() == "CA":
            res[l[21]][int(l[22:26])] = (float(l[30:38]), float(l[38:46]), float(l[46:54]), float(l[60:66]))
    return res
cs, cl = sys.argv[1], sys.argv[2]
print("model\tSOBIR1_residue\tSOBIR1_domain\tLRIP1_residue\tdistance_A\tpLDDT_SOBIR1\tpLDDT_LRIP1")
summ = collections.Counter()
for m in sys.argv[3:]:
    r = ca(m); s, t = r[cs], r[cl]
    assert s and t, f"chains {cs}/{cl} not found in {m}: {list(r)}"
    assert max(s) <= 632 and max(t) <= 410, "chain lengths do not match SOBIR1 (632) / LRIP1 (410) - check chain order"
    for i, (x1, y1, z1, p1) in s.items():
        for j, (x2, y2, z2, p2) in t.items():
            dd = math.dist((x1, y1, z1), (x2, y2, z2))
            if dd < 8 and p1 >= 50 and p2 >= 50:
                print(f"{m}\t{i}\t{dom(i)}\t{j}\t{dd:.2f}\t{p1:.1f}\t{p2:.1f}"); summ[(m, dom(i))] += 1
for (m, d), n in sorted(summ.items()): print(f"# {m}\t{d}\t{n} contacts", file=sys.stderr)
