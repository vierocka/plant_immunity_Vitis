import json, uuid, pathlib

def cell_id():
    return str(uuid.uuid4())[:8]

def md(text):
    return {"cell_type": "markdown", "id": cell_id(), "metadata": {}, "source": text}

def code(text):
    return {"cell_type": "code", "execution_count": None, "id": cell_id(),
            "metadata": {}, "outputs": [], "source": text}

cells = []

# ── TITLE ───────────────────────────────────────────────────────────────
cells.append(md(
"""# Grapevine Immunity Transcriptomics — Python Learning Toolkit

**Based on:** Hádlík M., Baránek M., Baránková K., Kovacova V.  
*"A SOBIR1-associated receptor network shows Rpv12-specific upregulation in grapevine downy mildew resistance"*  
bioRxiv 2025 · [doi:10.1101/2025.11.27.690962](https://doi.org/10.1101/2025.11.27.690962) · [NCBI BioProject PRJNA1358055](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1358055)

---

This notebook is a **Python learning companion** to the RNA-seq analysis pipeline from the paper above.
It covers QC, normalization, batch correction, differential gene expression, aggregated expression divergence,
and co-transcriptional network analysis — in Python, using **PyDESeq2** and **NetworkX**.

> ⚠️ **Published analysis vs teaching notebook**
> The paper was analysed in **R** using DESeq2 `rlog` transformation (stored in `data_files/Rlogs.csv`),
> which is the gold standard for variance stabilization in RNA-seq.
> This notebook uses PyDESeq2's **VST** (`vst_combat`) as the Python-accessible equivalent for teaching purposes.
> VST and rlog are conceptually analogous (both variance-stabilize; both followed by ComBat batch correction)
> but will produce numerically different values. For reproducing published figures exactly, use the R scripts
> in `02_Normalization_and_DGEA/`, `03_AED/`, and `06_PCNWA/`.

**Learning goals:**
1. Understand why and how batch effects are detected and corrected in RNA-seq
2. Apply size-factor normalization and VST in Python (PyDESeq2); compare all 8 normalization strategies
3. Perform differential expression analysis with PyDESeq2 and with a GLM on VST+ComBat values
4. Implement a custom permutation-based divergence metric (AED)
5. Build and analyse co-transcriptional networks with NetworkX

---

**Data files required** (paths are relative to the project root `plant_immunity_Vitis/`):
- `01_QC_and_Filtering/Batch_effects/reads_overview.csv` — mapping quality metrics per library
- `data_files/RawCounts.csv` — raw gene counts (26 169 genes × 36 samples)
- `data_files/SizeFactorNorm_ComBat_36samples.csv` — size-factor + log2 + ComBat values (pre-computed in R)
- `data_files/Patterns_01_DUE.csv` — DE pattern codes per gene (3 553 DE genes)
- `data_files/Vitis_gene_immunity_class.csv` — PTI / ETI / shared gene classification
- `data_files/node1_node2_ConfidenceScore_stringDB.csv` — STRING-db network edge list

VST-transformed + ComBat values (`vst_combat`) are computed on-the-fly by PyDESeq2 in Section 2.
"""
))

# ── TABLE OF CONTENTS ────────────────────────────────────────────────────
cells.append(md(
"""## Table of Contents

1. [Biological Background](#biology)
2. [Setup: packages and data loading](#setup)
3. [Section 1 — Batch Effect Detection from Read Mapping Statistics](#sec1)
4. [Section 2 — Normalization Strategies and ComBat Batch Correction](#sec2)
5. [Section 3 — Confirming Homoscedasticity of Normalised Data](#sec3)
6. [Section 4 — Differential Gene Expression Analysis (DGEA)](#sec4)
7. [Section 5 — Aggregated Expression Divergence (AED)](#sec5)
8. [Section 6 — Co-Transcriptional Network Analysis (PCNWA)](#sec6)
9. [Section 7 — Key Immunity Gene Network (Figure 5 equivalent)](#sec7)
"""
))

# ── BIOLOGY ──────────────────────────────────────────────────────────────
cells.append(md(
"""## 1. Biological Background <a name="biology"></a>

### The system

*Plasmopara viticola* (downy mildew) is a devastating oomycete pathogen of grapevine (*Vitis vinifera*),  
capable of destroying up to 75 % of wine production in a single season. Unlike true fungi, oomycetes  
are more closely related to brown algae and require different resistance mechanisms.

Breeding programs introduce **resistance loci (Rpv)** from wild *Vitis* relatives:

| Locus | Origin | NLR type | Mechanism |
|-------|---------|----------|-----------|
| **Rpv12** | *V. amurensis* | CNL (CC-NB-LRR) | Quantitative resistance, cell death |
| **Rpv1** | *Muscadinia rotundifolia* | TNL (TIR-NB-LRR) | Qualitative, strong hypersensitive response |
| **Rpv3** | *V. rupestris* | TNL | Quantitative, stilbene-associated |

Modern breeding stacks these loci (**pyramiding**): does stacking simply add the individual effects,  
or does it restructure the entire transcriptional immune programme?

### Experimental design

```
Genotype                 Time points       Replicates   → 36 libraries
─────────────────────────────────────────────────────────────────────
Rpv12       (single)     0 hpi (inoculation)
Rpv12+1     (double)     6 hpi (early response)   × 3 (A, B, C)
Rpv12+1+3   (triple)     24 hpi (late response)
Susceptible (Pinot Noir)
```

**Key question:** At which time point does each locus combination diverge from the susceptible transcriptome,  
and does pyramiding produce additive or restructured immune responses?

### Plant immunity tiers relevant to this study

- **PTI** (pattern-triggered immunity) — cell-surface PRRs (e.g. SOBIR1, MIK2) detect pathogen patterns
- **ETI** (effector-triggered immunity) — intracellular NLRs (Rpv12/Rpv1/Rpv3) detect pathogen effectors
- **Signal integration** — EDS1, SAG101, WRKY transcription factors link both tiers
- **Defense action** — SARD4 (SA biosynthesis), ACD6 (cell death), stilbene synthases

In *Arabidopsis*, PTI and ETI reinforce each other dynamically. Whether this occurs in grapevine  
and how locus pyramiding alters this coordination is the central biological question.
"""
))

# ── SETUP ────────────────────────────────────────────────────────────────
cells.append(md("""## 2. Setup <a name="setup"></a>

### 2.1 Package installation

Run the cell below if you haven't installed the required packages.  
`pydeseq2` replaces R's DESeq2. `neuroCombat` is the Python port of SVA's ComBat algorithm.
"""))

cells.append(code(
"""# Install required packages (uncomment and run once)
# !pip install pydeseq2 neuroCombat scikit-learn statsmodels scipy seaborn networkx plotnine adjustText
"""
))

cells.append(md("### 2.2 Imports"))

cells.append(code(
"""import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from scipy import stats
from scipy.stats import spearmanr, pearsonr, mannwhitneyu, ttest_ind
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import multipletests
import statsmodels.formula.api as smf
import statsmodels.api as sm
import networkx as nx
from itertools import combinations

# PyDESeq2 imports
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# ComBat batch correction
# neuroCombat is the standard Python port of the R SVA ComBat function
try:
    from neuroCombat import neuroCombat
    HAS_COMBAT = True
except ImportError:
    print("neuroCombat not installed. Run: pip install neuroCombat")
    print("Falling back to pre-computed batch-corrected matrices.")
    HAS_COMBAT = False

plt.rcParams['figure.dpi'] = 100
plt.rcParams['font.size'] = 11
print("All packages loaded.")
"""
))

cells.append(md("### 2.3 Define sample metadata and load data"))

cells.append(code(
"""# ── SAMPLE METADATA ──────────────────────────────────────────────────────────
# 36 samples: 4 genotypes × 3 time points × 3 biological replicates (A, B, C)
# Column order within each replicate block: 12 conditions (4 genotypes × 3 times)
conditions_per_rep = [
    "Rpv12.0", "Rpv12.6", "Rpv12.24",
    "Rpv12.1.0", "Rpv12.1.6", "Rpv12.1.24",
    "Rpv12.1.3.0", "Rpv12.1.3.6", "Rpv12.1.3.24",
    "Susceptible.0", "Susceptible.6", "Susceptible.24",
]
replicates = ["A", "B", "C"]

sample_names = [f"{cond}.{rep}" for rep in replicates for cond in conditions_per_rep]

# Batch assignment — 19 samples from Batch 1, 17 from Batch 2
# (determined by the sequencing facility; see read_checks.R for the original list)
BATCH1_SAMPLES = {
    "Rpv12.0.A", "Rpv12.0.B", "Rpv12.0.C",
    "Rpv12.1.0.A", "Rpv12.1.0.B", "Rpv12.1.0.C",
    "Rpv12.1.3.0.A", "Rpv12.1.3.0.B", "Rpv12.1.3.0.C",
    "Susceptible.0.B",
    "Rpv12.6.C", "Rpv12.1.6.B", "Rpv12.1.6.C",
    "Rpv12.1.3.6.A", "Rpv12.1.3.6.B",
    "Rpv12.1.24.B", "Rpv12.1.24.C",
    "Rpv12.1.3.24.C", "Susceptible.24.C",
}

def parse_sample(name):
    parts = name.rsplit(".", 1)
    cond, rep = parts[0], parts[1]
    if "Susceptible" in cond:
        genotype = "Susceptible"
    elif "Rpv12.1.3" in cond:
        genotype = "Rpv12+1+3"
    elif "Rpv12.1" in cond:
        genotype = "Rpv12+1"
    else:
        genotype = "Rpv12"
    time = int(cond.split(".")[-1])
    batch = "B1" if name in BATCH1_SAMPLES else "B2"
    return {"sample": name, "genotype": genotype, "time": time,
            "replicate": rep, "batch": batch, "condition": cond}

metadata = pd.DataFrame([parse_sample(s) for s in sample_names]).set_index("sample")
metadata["batch_int"] = (metadata["batch"] == "B2").astype(int)  # numeric batch for correlation

print("Sample metadata (first 12 rows):")
print(metadata.head(12).to_string())
print(f"\\nTotal: {len(metadata)} samples | Batch B1: {(metadata.batch=='B1').sum()} | B2: {(metadata.batch=='B2').sum()}")
"""
))

cells.append(code(
"""# ── COLOUR PALETTES ──────────────────────────────────────────────────────────
# Matching the R script colour conventions throughout the paper
GENOTYPE_COLORS = {
    "Rpv12":     "goldenrod",
    "Rpv12+1":   "salmon",
    "Rpv12+1+3": "cornflowerblue",
    "Susceptible": "dimgray",
}
BATCH_COLORS = {"B1": "salmon", "B2": "cornflowerblue"}
TIME_MARKERS  = {0: "o", 6: "^", 24: "s"}

def genotype_color(genotype): return GENOTYPE_COLORS.get(genotype, "black")
def batch_color(batch):       return BATCH_COLORS.get(batch, "black")
"""
))

cells.append(code(
"""# ── PATHS — relative to the project root (plant_immunity_Vitis/) ────────────
import pathlib, os

# Works regardless of where Jupyter was launched from, as long as
# the notebook lives at the project root (which it does on GitHub).
# Notebook lives in jupyter_nb/ — one level below the project root.
# Jupyter sets the kernel CWD to the notebook's directory, so .parent = project root.
PROJECT_ROOT = pathlib.Path().resolve().parent
# Fallback: if data_files is not found one level up, assume CWD is already the project root.
if not (PROJECT_ROOT / "data_files").exists():
    PROJECT_ROOT = pathlib.Path().resolve()
BASE   = PROJECT_ROOT / "data_files"
QC_DIR = PROJECT_ROOT / "01_QC_and_Filtering" / "Batch_effects"

# ── LOAD EXPRESSION DATA ─────────────────────────────────────────────────────
# Raw integer counts — load ALL genes before filtering (35 134 annotated protein-coding genes)
raw_counts_all = pd.read_csv(BASE / "RawCounts.csv", sep="\\t", index_col=0)
raw_counts_all = raw_counts_all[sample_names]   # enforce column order
print(f"Genes before filtering: {raw_counts_all.shape[0]:,}  |  Samples: {raw_counts_all.shape[1]}")

# ── FILTER: keep genes with sum ≥ 15 reads across all 36 samples ─────────────
# R equivalent: raw_counts_filtered <- VvitCountsMat[apply(VvitCountsMat,1,sum) >= 15, ]
CUTOFF = 15
raw_counts = raw_counts_all[raw_counts_all.sum(axis=1) >= CUTOFF]
print(f"After filtering (sum ≥ {CUTOFF}): {raw_counts.shape[0]:,} genes × {raw_counts.shape[1]} samples")
print(f"Removed: {raw_counts_all.shape[0] - raw_counts.shape[0]:,} low-count genes")

# Size-factor + log2 + ComBat (pre-computed in R; used for AED — preserves amplitude)
sf_combat = pd.read_csv(BASE / "SizeFactorNorm_ComBat_36samples.csv", index_col=0)[sample_names]
print(f"Size-factor + ComBat matrix: {sf_combat.shape}")

# VST-transformed + ComBat values are computed on-the-fly below (Section 2).
# The variable vst_combat (VST blind + ComBat) replaces the R-computed rlog+ComBat
# for all downstream analyses: PCA, DGEA, co-expression network, homoscedasticity check.
"""
))

cells.append(code(
"""# ── GENE-LEVEL COUNT DISTRIBUTIONS AND FILTER CUTOFF ─────────────────────────
# Show per-gene sum and mean across all 36 samples BEFORE filtering.
# Left panel : distribution of per-gene SUM  (total reads across all samples)
# Right panel: distribution of per-gene MEAN (average reads per sample)
# Red dashed line marks the sum ≥ 15 cutoff used for filtering.

CUTOFF = 15
gene_sum  = raw_counts_all.sum(axis=1)
gene_mean = raw_counts_all.mean(axis=1)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

# ── Panel A: sum per gene (log scale x-axis) ─────────────────────────────────
ax1.hist(np.log10(gene_sum + 1), bins=120, color="steelblue", edgecolor="none", alpha=0.85)
cutoff_log = np.log10(CUTOFF)
ax1.axvline(cutoff_log, color="firebrick", lw=2, ls="--",
            label=f"Cut-off: sum ≥ {CUTOFF}  (log₁₀ = {cutoff_log:.2f})")
n_removed = (gene_sum < CUTOFF).sum()
n_kept    = (gene_sum >= CUTOFF).sum()
ax1.text(cutoff_log + 0.05, ax1.get_ylim()[1] * 0.85,
         f"removed: {n_removed:,}\\nkept: {n_kept:,}",
         color="firebrick", fontsize=9, va="top")
ax1.set_xlabel("log₁₀(sum of raw reads + 1)", fontsize=11)
ax1.set_ylabel("Number of genes", fontsize=11)
ax1.set_title(f"A  Per-gene read sum  (all {raw_counts_all.shape[0]:,} genes)", fontsize=11)
ax1.legend(fontsize=9)

# ── Panel B: mean per gene (log scale x-axis) ────────────────────────────────
ax2.hist(np.log10(gene_mean + 1), bins=120, color="darkcyan", edgecolor="none", alpha=0.85)
mean_cutoff_log = np.log10(CUTOFF / raw_counts_all.shape[1])   # equivalent mean threshold
ax2.axvline(np.log10(CUTOFF / raw_counts_all.shape[1] + 1), color="firebrick", lw=2, ls="--",
            label=f"Equivalent mean cutoff\\n(sum/36 = {CUTOFF/raw_counts_all.shape[1]:.2f})")
ax2.set_xlabel("log₁₀(mean raw reads per sample + 1)", fontsize=11)
ax2.set_ylabel("Number of genes", fontsize=11)
ax2.set_title("B  Per-gene mean read depth", fontsize=11)
ax2.legend(fontsize=9)

plt.suptitle(
    f"Raw count distributions before low-count filtering  "
    f"({raw_counts_all.shape[0]:,} genes × {raw_counts_all.shape[1]} samples)",
    fontsize=12, fontweight="bold"
)
plt.tight_layout()
plt.savefig("gene_count_distributions.png", dpi=150, bbox_inches="tight")
plt.show()
print(f"Genes retained after sum ≥ {CUTOFF} filter: {n_kept:,} / {raw_counts_all.shape[0]:,} "
      f"({n_kept/raw_counts_all.shape[0]*100:.1f}%)")
"""
))

# ── SECTION 1 ────────────────────────────────────────────────────────────
cells.append(md(
"""## 3. Section 1 — Batch Effect Detection from Read Mapping Statistics <a name="sec1"></a>

> **Original R script:** `01_QC_and_Filtering/Batch_effects/read_checks.R`

### Why detect batch effects before touching expression values?

RNA-seq libraries sequenced in separate batches can differ in systematic, non-biological ways:  
adapter trimming efficiency, fragment size distribution, RNA degradation, reagent lot effects.  
These technical differences **propagate into expression quantification**, especially for lowly  
expressed genes where small absolute count differences translate to large fold changes.

The key insight: **if batches differ in how reads are processed, they will also differ in expression values.**  
Detecting batch effects at the read-mapping level (before any normalization) provides independent evidence  
that batch correction is needed.

### Test strategy

For each of 4 mapping-quality metrics across 36 libraries:
1. Fit a **GLM** with batch as predictor → F-test for batch effect
2. Compute **Spearman correlation** between metric and batch membership
3. Apply **Benjamini–Hochberg FDR** correction across the 4 metrics

A significant batch effect at this stage means the technical variation is large enough to affect  
downstream gene-level quantification — justifying ComBat correction.
"""
))

cells.append(code(
"""# Load the read-level QC table (36 libraries)
qc = pd.read_csv(QC_DIR / "reads_overview.csv", sep="\\t")
print(qc.head())
print(f"Shape: {qc.shape}")
"""
))

cells.append(code(
"""# Batch assignment for QC samples (same 19/17 split, indexed by sample ID)
NOVOGENE_ISOL_IDS = {
    "A11","B11","C11","A21","B21","C21","A31","B31","C31","B41",
    "C12","B22","C22","A32","B32","B24","C24","C34","C44"
}
qc["Batch"] = qc["ID"].apply(lambda x: "B1" if x in NOVOGENE_ISOL_IDS else "B2")
qc["Batch_int"] = (qc["Batch"] == "B2").astype(int)

print(qc[["ID", "Batch", "Forw_ok_perc", "Reverse_ok_perc",
          "both_ok_Perc", "dropped_perc"]].head(6))
"""
))

cells.append(code(
"""# ── GLM (OLS) + SPEARMAN TESTS FOR BATCH EFFECTS IN MAPPING METRICS ─────────
# R equivalent:
#   MoF <- glm(myTab$Forw_ok_perc ~ NovogeneIsolBatch)
#   cor.test(myTab$Forw_ok_perc, as.integer(as.factor(NovogeneIsolBatch)), method="spearman")

metrics = {
    "Forward-only survived": "Forw_ok_perc",
    "Reverse-only survived": "Reverse_ok_perc",
    "Both mates survived":   "both_ok_Perc",
    "Dropped (unmapped)":    "dropped_perc",
}

results = []
pvals_glm = []

for label, col in metrics.items():
    # GLM: response ~ batch (equivalent to R glm with gaussian family)
    model = smf.ols(f"{col} ~ C(Batch)", data=qc).fit()
    p_glm = model.f_pvalue          # F-test p-value for the batch predictor

    # Spearman correlation between metric and batch membership (0/1)
    rho, p_sp = spearmanr(qc[col], qc["Batch_int"])

    results.append({"Metric": label, "Spearman_rho": rho, "p_GLM": p_glm, "p_Spearman": p_sp})
    pvals_glm.append(p_glm)

results_df = pd.DataFrame(results)

# FDR correction (Benjamini–Hochberg) across the 4 tests
# R equivalent: p.adjust(c(...), method="fdr")
_, fdr_vals, _, _ = multipletests(pvals_glm, method='fdr_bh')
results_df["FDR"] = fdr_vals
results_df["Significant"] = results_df["FDR"] < 0.05

print("\\nBatch effect tests on mapping quality metrics:")
print(results_df.to_string(index=False, float_format="%.4f"))
print("\\nInterpretation: FDR < 0.05 → significant batch-dependent difference in that mapping metric.")
"""
))

cells.append(code(
"""# ── VISUALISATION ─────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(12, 8))
axes = axes.flatten()

for ax, (label, col) in zip(axes, metrics.items()):
    for i, row in qc.iterrows():
        c = BATCH_COLORS[row["Batch"]]
        ax.scatter(i + 1, row[col], color=c, s=40, zorder=3)

    # Add FDR annotation
    fdr = results_df.loc[results_df["Metric"] == label, "FDR"].values[0]
    sig = "**" if fdr < 0.01 else ("*" if fdr < 0.05 else "ns")
    ax.set_title(f"{label}  (FDR={fdr:.3f} {sig})", fontsize=10)
    ax.set_xlabel("Library index")
    ax.set_ylabel("Proportion (%)")

# Legend
b1_patch = mpatches.Patch(color="salmon",        label="Batch B1")
b2_patch = mpatches.Patch(color="cornflowerblue", label="Batch B2")
fig.legend(handles=[b1_patch, b2_patch], loc="lower center", ncol=2,
           bbox_to_anchor=(0.5, -0.02), fontsize=11)
fig.suptitle("Read-mapping quality metrics by sequencing batch\\n(Supplementary Figure 1)",
             fontsize=13, fontweight="bold")
plt.tight_layout()
plt.savefig("Supp_Fig1_batch_read_quality.png", dpi=150, bbox_inches="tight")
plt.show()
print("Figure saved.")
"""
))

cells.append(md(
"""### Interpretation

| Metric | FDR | Conclusion |
|--------|-----|-----------|
| Forward-only survived reads | 0.003 | **Significant** batch difference |
| Dropped (unmapped) reads | 0.004 | **Significant** batch difference |
| Reverse-only survived reads | 0.014 | **Significant** batch difference |
| Both mates survived (properly paired) | 0.077 | Not significant |

Three of four mapping metrics differed significantly between batches.  
The rate of properly paired reads (92–97 %) was comparable, indicating **similar overall data quality**.  
The differing forward/reverse-only rates suggest **RNA fragmentation or adapter-trimming differences** between  
library preparation runs — a common source of batch effects in multi-batch RNA-seq datasets.

**Conclusion: Batch correction is necessary** before any biological analysis.  
The type of batch correction (ComBat) and the normalization strategy to pair it with will be  
determined in Section 2.
"""
))

# ── SECTION 2 ────────────────────────────────────────────────────────────
cells.append(md(
"""## 4. Section 2 — Normalization Strategies and Batch Correction <a name="sec2"></a>

> **Original R script:** `01_QC_and_Filtering/Batch_effects/batch_effects.R`

### Why compare multiple normalization strategies?

Different normalizations capture different aspects of gene expression:
- **Raw counts** — absolute read depth; confounded by library size and batch
- **Size-factor normalization** (DESeq2/PyDESeq2) — corrects for library size; preserves fold-change amplitude
- **VST** (variance-stabilizing transformation) — stabilizes the mean–variance relationship; PyDESeq2's built-in replacement for R's rlog; best for PCA and correlation analyses

**ComBat** (empirical Bayes batch correction, from the R SVA package) removes batch effects *after*
normalization by standardizing per-gene mean and variance across batches.

### The 8 strategies compared (R paper + Python equivalents)

| Strategy | Normalization | ComBat? |
|----------|--------------|---------|
| 1 | Raw counts | No |
| 2 | Raw counts | Yes |
| 3 | Size-factor + log2 | No |
| 4 | **Size-factor + log2 + ComBat** | **Yes** → best for AED |
| 5 | VST blind (`use_design=False`) | No |
| 6 | **VST blind + ComBat** | **Yes** → best for network/PCA |
| 7 | VST design (`use_design=True`) | No |
| 8 | VST design + ComBat | Yes (over-correction risk) |

**Python teaching note:** Strategies 5–8 use PyDESeq2 VST as a Python-accessible approximation for the
R paper's rlog variants. The **published paper** used R/DESeq2 `rlog` (stored in `data_files/Rlogs.csv`).
`use_design=False` (blind VST) fits dispersion using only an intercept — analogous to R's `rlog(blind=TRUE)`.
`use_design=True` uses the full batch + condition design — analogous to R's `rlog(blind=FALSE)`.
Numerical values will differ between VST and rlog, but the variance-stabilizing effect is equivalent for
learning and exploration purposes.

**Key metric:** Spearman ρ between PC1 scores and batch membership.
A value close to 0 means batch does not dominate the primary expression variance axis.

### PyDESeq2 replaces DESeq2 for normalization and VST

PyDESeq2 is the official Python port of DESeq2. It uses the same negative binomial model,
median-of-ratios size factor estimation, and a built-in `vst()` method that stabilizes the
mean–variance relationship without requiring pre-computed R output files.
"""
))

cells.append(code(
"""# ── PYDESEQ2: SIZE FACTOR NORMALIZATION ──────────────────────────────────────
# PyDESeq2 requires:
#   - counts: pd.DataFrame (samples × genes) with integer values
#   - metadata: pd.DataFrame (samples × covariates) with 'condition' column

counts_for_deseq = raw_counts.T.astype(int)   # samples × genes (PyDESeq2 convention)
print(f"Count matrix shape for PyDESeq2: {counts_for_deseq.shape}")

# Design: include batch as covariate to model but not remove it
# (R equivalent: design = ~ batch + condition)
dds = DeseqDataSet(
    counts=counts_for_deseq,
    metadata=metadata[["condition", "batch"]],
    design_factors=["batch", "condition"],
    refit_cooks=True,
    quiet=True,
)
# Estimate size factors (library-size normalization)
# R equivalent: dds <- estimateSizeFactors(dds)
dds.fit_size_factors()
sf = dds.obs["size_factors"]   # pd.Series indexed by sample name (PyDESeq2 ≥ 0.4)
print("Size factors (first 12 samples):")
for s, f in sf.iloc[:12].items():
    print(f"  {s}: {f:.4f}")
"""
))

cells.append(code(
"""# Retrieve size-factor normalized counts and apply log2(x + 1)
# R equivalent: norm_counts_SF <- counts(dds, normalized=TRUE)
#               norm_counts_SF_log2 <- log2(norm_counts_SF + 1)
norm_counts = counts_for_deseq.div(pd.Series(sf), axis=0)   # divide each sample by its SF
norm_counts_log2 = np.log2(norm_counts + 1)

print(f"Normalized count matrix: {norm_counts_log2.shape}")
print(norm_counts_log2.iloc[:4, :4])
"""
))

cells.append(code(
"""# ── COMBAT BATCH CORRECTION ───────────────────────────────────────────────────
# neuroCombat is the Python port of R's sva::ComBat
# Input convention: GENES × SAMPLES (opposite of pandas/sklearn convention)
# R equivalent: expr_adj_SF_CB <- ComBat(norm_counts_SF_log2, batch=as.factor(BatchOriginIDs))

if HAS_COMBAT:
    dat = norm_counts_log2.T.values        # genes × samples
    batch_vec = metadata["batch"].values   # 'B1' or 'B2' per sample

    # neuroCombat internally estimates per-gene location and scale shifts per batch
    # and removes them using empirical Bayes shrinkage
    combat_result = neuroCombat(
        dat=dat,
        covars=pd.DataFrame({"batch": batch_vec, "condition": metadata["condition"].values},
                            index=metadata.index),
        batch_col="batch",
        categorical_cols=["condition"],
    )
    sf_combat_computed = pd.DataFrame(
        combat_result["data"],
        index=raw_counts.index,
        columns=sample_names,
    )
    print("ComBat correction applied to size-factor + log2 data.")
    print(sf_combat_computed.iloc[:3, :4])
else:
    # Use pre-computed values from the R pipeline
    sf_combat_computed = sf_combat.copy()
    print("Using pre-computed size-factor + ComBat matrix.")
"""
))

cells.append(code(
"""# ── VST: VARIANCE STABILIZING TRANSFORMATION ─────────────────────────────────
# PyDESeq2's built-in VST replaces DESeq2's rlog for Python workflows.
# Two variants mirror the 4 rlog strategies tested in the R paper:
#   use_design=False → blind VST (intercept-only dispersion fitting)
#   use_design=True  → design-aware VST (uses full batch + condition design)
# vst() auto-computes size factors and dispersion fits; full deseq2() not required.

print("Computing VST blind (use_design=False) ...")
dds_vst_blind = DeseqDataSet(
    counts=counts_for_deseq,
    metadata=metadata[["condition", "batch"]],
    design_factors=["batch", "condition"],
    quiet=True,
)
dds_vst_blind.vst(use_design=False)
vst_blind = pd.DataFrame(
    dds_vst_blind.layers["vst_counts"].T,   # obsm shape (samples × genes) → .T = genes × samples
    index=raw_counts.index,
    columns=sample_names,
)

print("Computing VST design (use_design=True) ...")
dds_vst_design = DeseqDataSet(
    counts=counts_for_deseq,
    metadata=metadata[["condition", "batch"]],
    design_factors=["batch", "condition"],
    quiet=True,
)
dds_vst_design.vst(use_design=True)
vst_design = pd.DataFrame(
    dds_vst_design.layers["vst_counts"].T,
    index=raw_counts.index,
    columns=sample_names,
)
print(f"VST blind  matrix: {vst_blind.shape}")
print(f"VST design matrix: {vst_design.shape}")
"""
))

cells.append(code(
"""# ── COMBAT ON VST VARIANTS ────────────────────────────────────────────────────
# Apply ComBat to each VST variant (same batch vector used for SF+log2+ComBat above).

if HAS_COMBAT:
    batch_vec = metadata["batch"].values
    covars_df = pd.DataFrame(
        {"batch": batch_vec, "condition": metadata["condition"].values},
        index=metadata.index,
    )

    # VST blind + ComBat  ← main downstream variable (replaces rlog+ComBat)
    res_vb = neuroCombat(dat=vst_blind.values, covars=covars_df,
                         batch_col="batch", categorical_cols=["condition"])
    vst_blind_combat = pd.DataFrame(res_vb["data"], index=raw_counts.index, columns=sample_names)

    # VST design + ComBat
    res_vd = neuroCombat(dat=vst_design.values, covars=covars_df,
                         batch_col="batch", categorical_cols=["condition"])
    vst_design_combat = pd.DataFrame(res_vd["data"], index=raw_counts.index, columns=sample_names)

    # Raw counts log2 + ComBat  (needed for 8-strategy PCA comparison)
    log2_raw = np.log2(raw_counts.astype(float) + 1)   # genes × samples
    res_raw = neuroCombat(dat=log2_raw.values, covars=covars_df,
                          batch_col="batch", categorical_cols=["condition"])
    raw_combat = pd.DataFrame(res_raw["data"], index=raw_counts.index, columns=sample_names)

    print("ComBat applied to VST blind, VST design, and raw log2.")
else:
    vst_blind_combat  = vst_blind.copy()
    vst_design_combat = vst_design.copy()
    raw_combat        = np.log2(raw_counts.astype(float) + 1)
    print("WARNING: neuroCombat not available — ComBat steps skipped.")

# Canonical downstream variable: VST blind + ComBat
# Used for PCA, DGEA (GLM approach), co-expression network, homoscedasticity check.
vst_combat = vst_blind_combat
print(f"\\nMain VST+ComBat matrix (vst_combat): {vst_combat.shape}")
"""
))

cells.append(code(
"""# ── PCA COMPARISON: ALL 8 STRATEGIES ─────────────────────────────────────────
# Reproduces the R paper's 4×4 PCA figure (PCA_all_methods_4x4.pdf).
# 8 strategies × 2 colour schemes (batch / genotype) = 16 panels.

def run_pca(expr_df, n_components=4):
    \"\"\"Run PCA on samples (expr_df: genes × samples → transpose internally).\"\"\"
    X = expr_df.T.values   # samples × genes
    pca = PCA(n_components=n_components)
    scores = pca.fit_transform(X)
    return scores, pca.explained_variance_ratio_ * 100

strategies = {
    "Raw counts":           raw_counts,
    "Raw counts + ComBat":  raw_combat,
    "SF + log2":            norm_counts_log2.T,        # genes × samples
    "SF + log2 + ComBat":   sf_combat_computed,
    "VST blind":            vst_blind,
    "VST blind + ComBat":   vst_blind_combat,
    "VST design":           vst_design,
    "VST design + ComBat":  vst_design_combat,
}

pca_results = {}
for name, df in strategies.items():
    scores, var_exp = run_pca(df)
    pca_results[name] = {"scores": scores, "var_exp": var_exp}

print("PCA variance explained (PC1–PC4) for each strategy:")
for name, res in pca_results.items():
    ve = res["var_exp"]
    print(f"  {name:25s}: PC1={ve[0]:.1f}% PC2={ve[1]:.1f}% PC3={ve[2]:.1f}% PC4={ve[3]:.1f}%")
"""
))

cells.append(code(
"""# ── BATCH CORRELATION: SPEARMAN(PC1, batch) ──────────────────────────────────
# Key diagnostic: how much does PC1 reflect batch rather than biology?
# R equivalent: cor.test(PCArc$x[,1], as.integer(colData[,2]), method="spearman")

print("Spearman correlation between PC1 scores and batch membership:")
print(f"  (rho close to 0 = batch does not dominate PC1 = batch effect removed)\\n")

batch_int = metadata["batch_int"].values  # 0=B1, 1=B2

for name, res in pca_results.items():
    pc1 = res["scores"][:, 0]
    rho, pval = spearmanr(pc1, batch_int)
    print(f"  {name:25s}: rho = {rho:+.3f}  (p={pval:.4f})")

print("\\nReference values from R script (paper):")
print("  Raw counts:             rho = -0.640")
print("  Raw counts + ComBat:    rho = -0.200")
print("  SF + log2:              rho = -0.597")
print("  SF + log2 + ComBat:     rho = -0.104  ← best for amplitude analyses (AED)")
print("  rlog (blind):           rho = -0.603  → VST blind should be similar")
print("  rlog + ComBat:          rho = -0.137  → VST blind + ComBat is the Python equivalent")
print("  rlog (design):          rho = -0.608  → VST design should be similar")
print("  rlog (design) + ComBat: rho = -0.142  → VST design + ComBat")
"""
))

cells.append(code(
"""# ── VISUALISATION: 4×4 PCA GRID — ALL 8 STRATEGIES × 2 COLOUR SCHEMES ───────
# Layout matches the R paper figure PCA_all_methods_4x4.pdf:
#   Rows 1-2  coloured by BATCH    (rows 1 = raw/SF strategies, row 2 = VST strategies)
#   Rows 3-4  coloured by GENOTYPE (same column ordering)
#   Spearman ρ(PC1, batch) shown in top-left of each batch-coloured panel.

strat_names_top = ["Raw counts", "Raw counts + ComBat",
                   "SF + log2", "SF + log2 + ComBat"]
strat_names_bot = ["VST blind", "VST blind + ComBat",
                   "VST design", "VST design + ComBat"]

fig, axes = plt.subplots(4, 4, figsize=(20, 18))

for col_set, strat_row in enumerate([strat_names_top, strat_names_bot]):
    for col, name in enumerate(strat_row):
        scores = pca_results[name]["scores"]
        ve     = pca_results[name]["var_exp"]
        rho, _ = spearmanr(scores[:, 0], batch_int)

        # Rows 0 & 1 → batch colour; rows 2 & 3 → genotype colour
        for color_row, (ax_row, color_by) in enumerate(
                [(col_set, "batch"), (col_set + 2, "genotype")]):
            ax = axes[ax_row, col]
            for i, sname in enumerate(sample_names):
                row = metadata.loc[sname]
                c  = batch_color(row["batch"]) if color_by == "batch" \
                     else genotype_color(row["genotype"])
                mk = TIME_MARKERS[row["time"]]
                ax.scatter(scores[i, 0], scores[i, 1],
                           color=c, marker=mk, s=40, alpha=0.85)

            ax.set_xlabel(f"PC1 ({ve[0]:.2f}%)", fontsize=8)
            ax.set_ylabel(f"PC2 ({ve[1]:.2f}%)", fontsize=8)
            ax.set_title(name, fontsize=8, pad=3)
            ax.tick_params(labelsize=7)

            # Annotate Spearman ρ on batch-coloured rows only
            if color_by == "batch":
                ax.text(0.03, 0.97, f"ρ = {rho:.3f}",
                        transform=ax.transAxes, fontsize=8,
                        va="top", ha="left", color="black")

# Row labels on leftmost column
row_labels = ["Batch (raw/SF)", "Batch (VST)",
              "Genotype (raw/SF)", "Genotype (VST)"]
for r, label in enumerate(row_labels):
    ax = axes[r, 0]
    ax.set_ylabel(f"{label}\\n{ax.get_ylabel()}", fontsize=8)

# Legends on rightmost column of each section
batch_handles = [mpatches.Patch(color=v, label=k) for k, v in BATCH_COLORS.items()]
geno_handles  = [mpatches.Patch(color=v, label=k) for k, v in GENOTYPE_COLORS.items()]
time_handles  = [plt.Line2D([0],[0], marker=mk, color="black", linestyle="None",
                             markersize=7, label=f"{t} hpi")
                 for t, mk in TIME_MARKERS.items()]
axes[0, -1].legend(handles=batch_handles + time_handles,
                   title="Batch / Time", fontsize=7, loc="lower right")
axes[2, -1].legend(handles=geno_handles + time_handles,
                   title="Genotype / Time", fontsize=7, loc="lower right")

plt.suptitle(
    "PCA under 8 normalisation strategies\\n"
    "Rows 1–2: coloured by batch (ρ = Spearman correlation PC1 vs batch)\\n"
    "Rows 3–4: coloured by genotype  |  Columns: Raw · Raw+ComBat · SF+log2 · SF+log2+ComBat  /  "
    "VST blind · VST blind+ComBat · VST design · VST design+ComBat",
    fontsize=11, fontweight="bold"
)
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig("Supp_Fig2_PCA_normalization_comparison.png", dpi=150, bbox_inches="tight")
plt.show()
print("Figure saved.")
"""
))

cells.append(md(
"""### Decision: Which normalization to use for which analysis?

| Analysis | R (published) | Python notebook (teaching) | Reason |
|----------|--------------|--------------------------|--------|
| **AED** | SF + log₂ + ComBat | SF + log₂ + ComBat | Preserves fold-change amplitude |
| **PCA, DGEA, network** | rlog + ComBat (`Rlogs.csv`) | VST blind + ComBat (`vst_combat`) | Variance-stabilized; homoscedastic |

**Key finding:** DESeq2's batch-as-covariate design (included in the GLM) did **not** effectively
remove batch effects (ρ ≈ −0.60 regardless). **ComBat was required** in all cases.
This is a known limitation: DESeq2/PyDESeq2 adjusts for batch in the GLM but does not correct the
expression matrix itself — so PCA and correlation-based analyses still see the batch signal.

> 💡 **Biological validation:** After batch correction, PC1 separates genotypes (Rpv12+1+3 vs others)
> rather than batches, confirming that the biologically meaningful signal is preserved.
"""
))

# ── SECTION 3 ────────────────────────────────────────────────────────────
cells.append(md(
"""## 5. Section 3 — Confirming Homoscedasticity <a name="sec3"></a>

> **Original R script:** `02_Normalization_and_DGEA/Heteroscedasticity_check.R`

### Why check homoscedasticity?

**Pearson correlation** — the basis of co-transcriptional network analysis — assumes that the  
variance of a variable is roughly constant across its range of values (**homoscedasticity**).  
In raw RNA-seq data, variance scales strongly with mean expression (low-count genes are noisier).

**VST** (variance-stabilizing transformation) is specifically designed to stabilize this
mean–variance relationship. If VST + ComBat data are truly homoscedastic:
- Spearman/Pearson correlation between gene means and gene variances should be near **0**
- The slope of the linear fit (variance ~ mean) should be near **0**

This is tested both globally (all 36 samples) and per experimental condition (3 replicates).

### What to expect

| Data | Global Spearman r (mean vs variance) |
|------|--------------------------------------|
| Raw counts | +0.6 to +0.9 (strong positive — high-expressed genes are also most variable) |
| SF + ComBat | −0.68 (negative — over-correction removes variance non-uniformly) |
| **VST + ComBat** | **≈ 0** (near zero — truly variance-stabilized) ✓ |
"""
))

cells.append(code(
"""def check_heteroscedasticity(expr_df, condition_vec, method="spearman", plot=False, n_plot=None):
    \"\"\"
    Assess mean–variance relationship globally and per experimental condition.

    Parameters
    ----------
    expr_df     : pd.DataFrame (genes × samples)
    condition_vec : array-like, condition label for each sample
    method      : 'spearman' or 'pearson'
    plot        : whether to draw per-condition scatter plots
    n_plot      : max conditions to plot (None = all)

    Returns
    -------
    dict with global and per-condition results
    \"\"\"
    corr_fn = spearmanr if method == "spearman" else pearsonr

    # Global (all samples pooled)
    gene_means = expr_df.mean(axis=1)
    gene_vars  = expr_df.var(axis=1, ddof=1)

    r_global, p_global = corr_fn(gene_means, gene_vars)
    slope_global = np.polyfit(gene_means, gene_vars, 1)[0]

    print(f"=== GLOBAL MEAN–VARIANCE RELATIONSHIP ({method.capitalize()}) ===")
    print(f"  r = {r_global:.3f}  |  slope = {slope_global:.4f}  |  p = {p_global:.4e}")

    # Per condition
    conditions = pd.Series(condition_vec).unique()
    per_cond = []
    for cond in sorted(conditions):
        idx = [i for i, c in enumerate(condition_vec) if c == cond]
        if len(idx) < 2:
            continue
        sub = expr_df.iloc[:, idx]
        means_c = sub.mean(axis=1)
        vars_c  = sub.var(axis=1, ddof=1)
        r_c, p_c = corr_fn(means_c, vars_c)
        slope_c  = np.polyfit(means_c, vars_c, 1)[0]
        per_cond.append({"condition": cond, "r": r_c, "slope": slope_c, "p": p_c})

    per_cond_df = pd.DataFrame(per_cond)
    print("\\n=== PER-CONDITION RELATIONSHIPS ===")
    print(per_cond_df.to_string(index=False, float_format="%.3f"))

    return {"global": {"r": r_global, "slope": slope_global}, "per_condition": per_cond_df}
"""
))

cells.append(code(
"""# ── HOMOSCEDASTICITY CHECK: STATISTICS + VISUALISATION ───────────────────────
condition_vec = metadata["condition"].values

print("### VST blind + ComBat (used for PCA, DGEA, network analysis) ###")
res_vst = check_heteroscedasticity(vst_combat, condition_vec, method="spearman")

print("\\n### Size-factor + log2 + ComBat (used for AED) ###")
res_sf = check_heteroscedasticity(sf_combat, condition_vec, method="spearman")

# ── FIGURE: 2 × 2 panels ──────────────────────────────────────────────────────
# Top row   : global mean–variance scatter (all 26 169 genes, pooled across samples)
# Bottom row: per-condition Spearman r  (bar chart, one bar per experimental condition)

datasets = [
    ("VST blind + ComBat", vst_combat,  res_vst,  "steelblue"),
    ("SF + log₂ + ComBat", sf_combat,   res_sf,   "darkorange"),
]

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

for col, (label, df, res, color) in enumerate(datasets):
    # ── Top: global mean–variance scatter ────────────────────────────────────
    ax_top = axes[0, col]
    gm = df.mean(axis=1)
    gv = df.var(axis=1, ddof=1)
    ax_top.scatter(gm, gv, s=1, alpha=0.2, color=color, rasterized=True)
    coefs = np.polyfit(gm, gv, 1)
    xl = np.linspace(gm.min(), gm.max(), 200)
    ax_top.plot(xl, np.polyval(coefs, xl), color="firebrick", lw=2,
                label=f"slope = {coefs[0]:.4f}")
    rho_g, _ = spearmanr(gm, gv)
    ax_top.set_title(f"{label}\\nGlobal Spearman ρ = {rho_g:.3f}", fontsize=11)
    ax_top.set_xlabel("Mean expression (across 36 samples)", fontsize=10)
    ax_top.set_ylabel("Variance", fontsize=10)
    ax_top.legend(fontsize=9)

    # ── Bottom: per-condition r bar chart ─────────────────────────────────────
    ax_bot = axes[1, col]
    per_cond = res["per_condition"].sort_values("r")
    bar_colors = ["firebrick" if abs(r) > 0.3 else color
                  for r in per_cond["r"]]
    ax_bot.barh(per_cond["condition"], per_cond["r"],
                color=bar_colors, edgecolor="none", height=0.7)
    ax_bot.axvline(0, color="black", lw=1)
    ax_bot.axvline( 0.3, color="firebrick", lw=1, ls="--", alpha=0.6)
    ax_bot.axvline(-0.3, color="firebrick", lw=1, ls="--", alpha=0.6,
                   label="|ρ| > 0.3 threshold")
    ax_bot.set_xlabel("Spearman ρ (mean vs variance, per condition)", fontsize=10)
    ax_bot.set_title("Per-condition mean–variance correlation", fontsize=10)
    ax_bot.tick_params(axis="y", labelsize=8)
    ax_bot.legend(fontsize=8)

plt.suptitle(
    "Homoscedasticity check\\n"
    "Top: global mean–variance scatter  |  "
    "Bottom: per-condition Spearman ρ (|ρ| ≈ 0 → homoscedastic)",
    fontsize=12, fontweight="bold"
)
plt.tight_layout()
plt.savefig("heteroscedasticity_check.png", dpi=150, bbox_inches="tight")
plt.show()
"""
))

cells.append(md(
"""### Conclusion

| Data | Global Spearman r | Interpretation |
|------|------------------|---------------|
| VST + ComBat | **≈ 0** | Homoscedastic ✓ — suitable for Pearson correlations |
| SF + ComBat | ≈ −0.68 | Heteroscedastic — variance decreases with mean (expected) |

- VST + ComBat data pass the homoscedasticity check → **approved for network analysis and PCA**
- SF + ComBat retain the expected mean–variance structure → **approved for AED (amplitude-sensitive)**

**Residual batch check (from DGEA.R):**
Per-gene Wilcoxon tests comparing B1 vs B2 after VST+ComBat should find few genes with FDR < 0.05
and **none with |log₂FC| > 1** — confirming no biologically relevant batch signal remains.
"""
))

# ── SECTION 4 ────────────────────────────────────────────────────────────
cells.append(md(
"""## 6. Section 4 — Differential Gene Expression Analysis (DGEA) <a name="sec4"></a>

> **Original R script:** `02_Normalization_and_DGEA/DGEA.R`

### Two approaches compared

**Approach A — PyDESeq2 (proper, recommended):**  
Uses raw integer counts with the negative binomial model.  
Accounts for count overdispersion; applies shrinkage estimators for dispersion and fold change.  
This is the gold standard for RNA-seq differential expression.

**Approach B — GLM on VST + ComBat (paper's approach adapted to Python):**
Uses batch-corrected variance-stabilized values. Fits a Gaussian GLM per gene per pairwise comparison.
Simpler and computationally cheaper; appropriate when the normalization step already handles
library-size and batch effects. This is the Python equivalent of what the paper used to define the 3 553 DEGs.

**Why did the paper use Approach B?**  
When a large co-expression network analysis requires consistent normalization across all  
subsequent steps (PCA, clustering, correlation), it's practical to perform DGEA on the same  
normalized values rather than switching back to raw counts. However, PyDESeq2 is statistically  
more rigorous for the DGEA step itself.

**DEG criteria:** |log₂FC| > 1 and FDR < 0.05 in any of the 9 comparisons (3 genotypes × 3 times)
"""
))

cells.append(code(
"""# ── APPROACH A: PYDESEQ2 DGEA ─────────────────────────────────────────────────
# We run PyDESeq2 on raw counts with condition as the design factor.
# Reference level: Susceptible.0 (as in the R script: relevel to "Susceptible.0")
#
# Note: Running full PyDESeq2 on 26 169 genes takes ~5–15 minutes.
# The results below can be compared to the pre-computed DGEA files.

print("Setting up PyDESeq2 for full DGEA...")
print("(Using condition design without explicit batch modelling —")
print("batch effects are already removed in VST+ComBat values.)")
print("Running DESeq2 on raw counts with batch + condition design...")

dds_full = DeseqDataSet(
    counts=counts_for_deseq,
    metadata=metadata[["condition", "batch"]],
    design_factors=["batch", "condition"],
    refit_cooks=True,
    quiet=True,
)
dds_full.deseq2()
print("DESeq2 model fitted.")
"""
))

cells.append(code(
"""# Extract results for Rpv12 vs Susceptible at 0 hpi
# R equivalent: results(dds, contrast=c("condition","Rpv12.0","Susceptible.0"))

stat_rpv12_0 = DeseqStats(
    dds_full,
    contrast=["condition", "Rpv12.0", "Susceptible.0"],
    alpha=0.05,
)
stat_rpv12_0.summary()
res_rpv12_0 = stat_rpv12_0.results_df
sig_rpv12_0 = res_rpv12_0[(res_rpv12_0["padj"] < 0.05) & (res_rpv12_0["log2FoldChange"].abs() > 1)]
print(f"\\nPyDESeq2 — Rpv12 vs Susceptible at 0 hpi:")
print(f"  Up-regulated:   {(sig_rpv12_0['log2FoldChange'] > 0).sum()}")
print(f"  Down-regulated: {(sig_rpv12_0['log2FoldChange'] < 0).sum()}")
"""
))

cells.append(code(
"""# ── APPROACH B: GLM ON VST + ComBat (PAPER'S APPROACH ADAPTED) ──────────────
# For each gene, for each pairwise comparison (resistant genotype vs susceptible, per time point):
#   1. Compute fold-change = mean(VST, resistant) - mean(VST, susceptible)
#   2. Fit Gaussian GLM; extract F-test p-value
#   3. Apply BH-FDR across all genes
#
# This is computationally faster than PyDESeq2 and mirrors what the paper used (rlog+ComBat).
# It treats VST values as approximately normal (reasonable after variance stabilization).

# Column index mapping (0-based, within the sample_names list)
# Replicate blocks: A = cols 0–11, B = cols 12–23, C = cols 24–35
# Within each block: Rpv12(0,6,24), Rpv12+1(0,6,24), Rpv12+1+3(0,6,24), Susceptible(0,6,24)

COND_IDX = {
    "Rpv12.0":         [0,  12, 24],
    "Rpv12.6":         [1,  13, 25],
    "Rpv12.24":        [2,  14, 26],
    "Rpv12.1.0":       [3,  15, 27],
    "Rpv12.1.6":       [4,  16, 28],
    "Rpv12.1.24":      [5,  17, 29],
    "Rpv12.1.3.0":     [6,  18, 30],
    "Rpv12.1.3.6":     [7,  19, 31],
    "Rpv12.1.3.24":    [8,  20, 32],
    "Susceptible.0":   [9,  21, 33],
    "Susceptible.6":   [10, 22, 34],
    "Susceptible.24":  [11, 23, 35],
}

# Helper: run GLM-based DGEA for one pairwise comparison
def dgea_glm(expr_mat, cond_a_idx, cond_b_idx, label_a="A", label_b="B"):
    \"\"\"
    Per-gene Gaussian GLM F-test + Wilcoxon test.
    Returns DataFrame with log2FC, p_glm, p_wilcox, padj_glm, padj_wilcox.

    expr_mat : pd.DataFrame (genes × samples)
    cond_a_idx : list of column indices for condition A (resistant genotype)
    cond_b_idx : list of column indices for condition B (susceptible reference)
    \"\"\"
    genes = expr_mat.index
    A = expr_mat.iloc[:, cond_a_idx].values   # shape: (n_genes, 3)
    B = expr_mat.iloc[:, cond_b_idx].values

    log2fc = A.mean(axis=1) - B.mean(axis=1)   # difference of means on VST scale

    # Vectorised t-test (Welch) — equivalent to gaussian GLM F-test for 2 groups
    # scipy.stats.ttest_ind returns p-values per gene
    _, p_ttest = ttest_ind(A, B, axis=1, equal_var=False)

    # Wilcoxon rank-sum (Mann–Whitney U) — non-parametric alternative
    p_wilcox = np.array([
        mannwhitneyu(A[i], B[i], alternative="two-sided").pvalue
        for i in range(len(genes))
    ])

    # BH-FDR correction
    _, padj_t, _, _ = multipletests(p_ttest,  method="fdr_bh")
    _, padj_w, _, _ = multipletests(p_wilcox, method="fdr_bh")

    return pd.DataFrame({
        "log2FC": log2fc,
        "p_ttest": p_ttest,
        "p_wilcox": p_wilcox,
        "padj_ttest": padj_t,
        "padj_wilcox": padj_w,
    }, index=genes)

print("GLM-based DGEA function defined.")
print("Computing 9 pairwise comparisons (3 genotypes × 3 time points vs Susceptible)...")
"""
))

cells.append(code(
"""# Run all 9 comparisons
comparisons = [
    ("Rpv12",     "Rpv12.0",     "Susceptible.0"),
    ("Rpv12",     "Rpv12.6",     "Susceptible.6"),
    ("Rpv12",     "Rpv12.24",    "Susceptible.24"),
    ("Rpv12+1",   "Rpv12.1.0",   "Susceptible.0"),
    ("Rpv12+1",   "Rpv12.1.6",   "Susceptible.6"),
    ("Rpv12+1",   "Rpv12.1.24",  "Susceptible.24"),
    ("Rpv12+1+3", "Rpv12.1.3.0", "Susceptible.0"),
    ("Rpv12+1+3", "Rpv12.1.3.6", "Susceptible.6"),
    ("Rpv12+1+3", "Rpv12.1.3.24","Susceptible.24"),
]

dgea_results = {}
summary_rows = []

for geno, cond_a, cond_b in comparisons:
    time = cond_a.split(".")[-1]
    key  = f"{geno}_{time}hpi"
    res  = dgea_glm(vst_combat, COND_IDX[cond_a], COND_IDX[cond_b], cond_a, cond_b)
    dgea_results[key] = res

    sig = res[(res["padj_ttest"] < 0.05) & (res["log2FC"].abs() > 1)]
    up   = (sig["log2FC"] > 0).sum()
    down = (sig["log2FC"] < 0).sum()
    summary_rows.append({"Genotype": geno, "Time (hpi)": time, "Up": up, "Down": down, "Total DEGs": up+down})

summary_df = pd.DataFrame(summary_rows)
print("\\nDEG counts (|log2FC| > 1, FDR < 0.05):")
print(summary_df.to_string(index=False))
print("\\nExpected from paper:")
print("  Rpv12 0h:      246 up / 125 down")
print("  Rpv12 6h:      130 up / 713 down")
print("  Rpv12+1+3 6h:  396 up / 1421 down (largest response)")
"""
))

cells.append(code(
"""# ── PCA OF BATCH-CORRECTED EXPRESSION VALUES ──────────────────────────────────
# Reproduce Figure 2B from the paper: PCA of all 26 169 genes, VST + ComBat
# PC1 separates by genotype; PC2 separates by time

scores_vst, var_vst = run_pca(vst_combat, n_components=4)
ve = var_vst

fig, axes = plt.subplots(1, 2, figsize=(14, 6))
titles = ["PC1 vs PC2", "PC3 vs PC4"]
pc_pairs = [(0, 1), (2, 3)]

for ax, (pc_x, pc_y), title in zip(axes, pc_pairs, titles):
    for i, sname in enumerate(sample_names):
        row = metadata.loc[sname]
        c  = genotype_color(row["genotype"])
        mk = TIME_MARKERS[row["time"]]
        ax.scatter(scores_vst[i, pc_x], scores_vst[i, pc_y],
                   color=c, marker=mk, s=60, alpha=0.9, edgecolors="white", lw=0.3)

    ax.set_xlabel(f"PC{pc_x+1} ({ve[pc_x]:.2f}%)", fontsize=12)
    ax.set_ylabel(f"PC{pc_y+1} ({ve[pc_y]:.2f}%)", fontsize=12)
    ax.axhline(0, lw=0.8, ls="--", color="cornflowerblue", alpha=0.5)
    ax.axvline(0, lw=0.8, ls="--", color="salmon", alpha=0.5)
    ax.set_title(title, fontsize=11)

# Legends
geno_handles = [mpatches.Patch(color=v, label=k) for k, v in GENOTYPE_COLORS.items()]
time_handles = [plt.Line2D([0],[0], marker=mk, color="black", linestyle="None",
                            markersize=8, label=f"{t} hpi") for t, mk in TIME_MARKERS.items()]
axes[0].legend(handles=geno_handles, title="Genotype", loc="upper left", fontsize=9)
axes[1].legend(handles=time_handles, title="Time",     loc="upper left", fontsize=9)

plt.suptitle("PCA of 26 169 genes (VST + ComBat) — Figure 2B equivalent",
             fontsize=13, fontweight="bold")
plt.tight_layout()
plt.savefig("Fig2B_PCA_genotype_time.png", dpi=150, bbox_inches="tight")
plt.show()
"""
))

cells.append(code(
"""# ── PC LOADING ANALYSIS: TOP 50 CONTRIBUTING GENES ──────────────────────────
# Genes driving each PC reveal the biological processes separating experimental groups.
# R equivalent: loadings <- BigPCA$rotation[,1:4]
#               contributions1 <- (squared_loadings[,1]/sum(squared_loadings[,1])) * 100

pca_full = PCA(n_components=4)
pca_full.fit(vst_combat.T.values)   # fit on samples × genes

# loadings: each column is a PC, each row is a gene's contribution
loadings = pd.DataFrame(
    pca_full.components_.T,     # genes × PCs
    index=vst_combat.index,
    columns=["PC1", "PC2", "PC3", "PC4"],
)
contributions = (loadings ** 2).div((loadings ** 2).sum(axis=0), axis=1) * 100

for pc in ["PC1", "PC2", "PC3", "PC4"]:
    top50 = contributions[pc].nlargest(50)
    pct_captured = top50.sum()
    print(f"{pc}: top 50 genes capture {pct_captured:.1f}% of variance in this PC")

print("\\nTop 10 genes contributing to PC1 (paper: enriched for signaling, FDR=1.3e-10):")
print(contributions["PC1"].nlargest(10).to_string())
"""
))

# ── SECTION 5 ────────────────────────────────────────────────────────────
cells.append(md(
"""## 7. Section 5 — Aggregated Expression Divergence (AED) <a name="sec5"></a>

> **Original R script:** `03_AED/DGE_bySFandCB_divergence.R`

### What is AED and why use it?

Standard differential expression analysis tests genes **individually** and reports a list of  
significant genes. But it cannot tell you: **how different is the overall transcriptome state**  
of a resistant genotype from the susceptible control?

**Aggregated Expression Divergence (AED)** is a global metric:

$$\\text{AED} = \\frac{1}{G} \\sum_{g=1}^{G} (\\bar{x}_{g,\\text{resistant}} - \\bar{x}_{g,\\text{susceptible}})^2$$

where $G$ = number of genes, $\\bar{x}$ = mean expression across replicates.

AED measures the **mean squared deviation** of a resistant genotype's transcriptome from the  
susceptible reference — across all genes simultaneously.

### Why size-factor + ComBat values (not VST)?

AED is sensitive to absolute expression **amplitude** (the magnitude of gene-level deviations).
VST (like rlog) shrinks large fold-changes (it's variance-stabilizing), which would compress
AED values and make large genuine differences look small. Size-factor normalization preserves
the amplitude, so large biological differences produce large AED contributions.

### Statistical testing via permutation

With only 3 replicates per condition, parametric tests have low power. Instead, we use a  
**permutation test**: under the null hypothesis (H₀) that the resistant genotype is no different  
from susceptible, any 3 samples from the same time-point are interchangeable.

There are $\\binom{12}{3} = 220$ ways to choose 3 samples from 12 same-timepoint samples.  
We compute AED for all 220 random triplets vs the susceptible mean → null distribution.  
The **empirical p-value** = proportion of null AEDs ≥ observed AED.

The minimum attainable p-value is $1/220 = 0.0045$.  
After BH-FDR correction across 9 tests, the minimum attainable FDR = **0.024**.
"""
))

cells.append(code(
"""# ── LOAD SIZE-FACTOR + ComBat DATA ───────────────────────────────────────────
# This is the normalization chosen for AED: preserves fold-change amplitude
expr_aed = sf_combat.copy()   # genes × samples (pre-loaded earlier)
print(f"Expression matrix for AED: {expr_aed.shape}")

# Mean expression per condition group (for each genotype × time combination)
def group_mean(expr, col_indices):
    return expr.iloc[:, col_indices].mean(axis=1)

# Susceptible reference means per time point
mean_S0  = group_mean(expr_aed, COND_IDX["Susceptible.0"])
mean_S6  = group_mean(expr_aed, COND_IDX["Susceptible.6"])
mean_S24 = group_mean(expr_aed, COND_IDX["Susceptible.24"])

# Resistant genotype means
groups_by_time = {
    0: {
        "Rpv12":     group_mean(expr_aed, COND_IDX["Rpv12.0"]),
        "Rpv12+1":   group_mean(expr_aed, COND_IDX["Rpv12.1.0"]),
        "Rpv12+1+3": group_mean(expr_aed, COND_IDX["Rpv12.1.3.0"]),
    },
    6: {
        "Rpv12":     group_mean(expr_aed, COND_IDX["Rpv12.6"]),
        "Rpv12+1":   group_mean(expr_aed, COND_IDX["Rpv12.1.6"]),
        "Rpv12+1+3": group_mean(expr_aed, COND_IDX["Rpv12.1.3.6"]),
    },
    24: {
        "Rpv12":     group_mean(expr_aed, COND_IDX["Rpv12.24"]),
        "Rpv12+1":   group_mean(expr_aed, COND_IDX["Rpv12.1.24"]),
        "Rpv12+1+3": group_mean(expr_aed, COND_IDX["Rpv12.1.3.24"]),
    },
}
ref_means = {0: mean_S0, 6: mean_S6, 24: mean_S24}

# Compute AED for each genotype × time
def compute_aed(mean_resistant, mean_susceptible):
    \"\"\"Mean squared gene-wise deviation (AED).\"\"\"
    return np.mean((mean_resistant - mean_susceptible) ** 2)

aed_obs = {}
for t, geno_dict in groups_by_time.items():
    for geno, mean_geno in geno_dict.items():
        aed_obs[(geno, t)] = compute_aed(mean_geno, ref_means[t])

print("Observed AED values:")
for (geno, t), aed in aed_obs.items():
    print(f"  {geno:12s} @ {t:2d} hpi : AED = {aed:.4f}")
"""
))

cells.append(code(
"""# ── PERMUTATION TEST: BUILD NULL DISTRIBUTIONS ────────────────────────────────
# For each time point, draw all C(12,3)=220 triplets from the 12 same-timepoint samples.
# Compute AED for each triplet vs the susceptible mean.
#
# R equivalent:
#   samples0 <- sample(seq(from=1, to=36, by=3))     # column indices of 0-hpi samples
#   combos0  <- combn(samples0, 3)                    # 3x220 matrix
#   AD_Ho_t0 <- sapply(1:220, function(i)
#       mean(abs(apply(expr_adj[, combos0[,i]], 1, mean) - mean_Go_0)^2))

import numpy as np
from itertools import combinations

np.random.seed(42)  # reproducibility

def null_aed_distribution(expr, all_time_col_indices, ref_mean):
    \"\"\"
    Compute AED for all C(n,3) triplets chosen from all_time_col_indices.
    Returns array of 220 AED values forming the null distribution.
    \"\"\"
    combos = list(combinations(all_time_col_indices, 3))
    assert len(combos) == 220, f"Expected 220, got {len(combos)}"
    null_aeds = []
    for triplet in combos:
        mean_triplet = expr.iloc[:, list(triplet)].mean(axis=1)
        null_aeds.append(compute_aed(mean_triplet, ref_mean))
    return np.array(null_aeds)

# All column indices for each time point (12 samples per time point)
time_col_indices = {
    0:  [i for i, s in enumerate(sample_names) if ".0." in s or s.endswith(".0.A") or
         metadata.loc[s, "time"] == 0],
    6:  [i for i, s in enumerate(sample_names) if metadata.loc[s, "time"] == 6],
    24: [i for i, s in enumerate(sample_names) if metadata.loc[s, "time"] == 24],
}
for t in [0, 6, 24]:
    time_col_indices[t] = [i for i, s in enumerate(sample_names)
                           if metadata.loc[s, "time"] == t]
    assert len(time_col_indices[t]) == 12, f"t={t}: {len(time_col_indices[t])} samples (expected 12)"

print("Building null distributions (220 combinations × 3 time points)...")
null_dists = {}
for t in [0, 6, 24]:
    null_dists[t] = null_aed_distribution(expr_aed, time_col_indices[t], ref_means[t])
    print(f"  t={t:2d} hpi: mean(null) = {null_dists[t].mean():.4f}  "
          f"max(null) = {null_dists[t].max():.4f}")
"""
))

cells.append(code(
"""# ── EMPIRICAL P-VALUES AND FDR ────────────────────────────────────────────────
# Empirical p-value = (# null AEDs >= observed AED + 1) / (220 + 1)
# The "+1" is a continuity correction that prevents p=0 (Lariviere & Langelier, 2001)
# R equivalent: p_emp <- (sum(AD_Ho_t0 >= AggrDiv_RPV12_0) + 1) / (length(AD_Ho_t0) + 1)

emp_pvals = {}
for (geno, t), obs_aed in aed_obs.items():
    null = null_dists[t]
    p = (np.sum(null >= obs_aed) + 1) / (len(null) + 1)
    emp_pvals[(geno, t)] = p

# Collect in the order: t=0 (3 genotypes), t=6, t=24
ordered_keys = [(g, t) for t in [0, 6, 24] for g in ["Rpv12", "Rpv12+1", "Rpv12+1+3"]]
p_array = [emp_pvals[k] for k in ordered_keys]

# BH-FDR correction across all 9 tests
_, fdr_array, _, _ = multipletests(p_array, method="fdr_bh")
fdr_dict = dict(zip(ordered_keys, fdr_array))

print("\\nEmpirical p-values and FDR-corrected values:")
print(f"{'Genotype':12s} {'Time':5s} {'AED':8s} {'p_emp':8s} {'FDR':8s} {'Sig':4s}")
for k in ordered_keys:
    geno, t = k
    sig = "*" if fdr_dict[k] < 0.05 else ""
    print(f"{geno:12s} {t:3d}h   {aed_obs[k]:.4f}   {emp_pvals[k]:.4f}   {fdr_dict[k]:.4f}   {sig}")

print("\\nExpected from paper:")
print("  Rpv12+1   @ 0 hpi : FDR=0.024*")
print("  Rpv12+1+3 @ 0 hpi : FDR=0.024*")
print("  Rpv12+1+3 @ 6 hpi : FDR=0.024*")
print("  Rpv12     @ 24 hpi: FDR=0.024*")
print("  Rpv12+1+3 @ 24 hpi: FDR=0.024*")
"""
))

cells.append(code(
"""# ── VISUALISATION: DENSITY PLOTS WITH OBSERVED AED VALUES ────────────────────
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
time_labels = ["0 hpi", "6 hpi", "24 hpi"]
geno_list   = ["Rpv12", "Rpv12+1", "Rpv12+1+3"]
geno_pch    = {0: "o", 6: "^", 24: "s"}

for ax, t, tlabel in zip(axes, [0, 6, 24], time_labels):
    null = null_dists[t]
    # KDE of null distribution
    from scipy.stats import gaussian_kde
    kde = gaussian_kde(null)
    x_range = np.linspace(null.min() - 0.1, null.max() + 0.3, 300)
    ax.plot(x_range, kde(x_range), color="dimgray", lw=2, label="Null distribution")
    ax.fill_between(x_range, kde(x_range), alpha=0.15, color="dimgray")

    for geno in geno_list:
        obs = aed_obs[(geno, t)]
        c   = GENOTYPE_COLORS[geno]
        fdr = fdr_dict[(geno, t)]
        sig_label = f" (FDR={fdr:.3f}{'*' if fdr<0.05 else ''})"
        ax.axvline(obs, color=c, lw=2, label=f"{geno}{sig_label}")
        ax.scatter(obs, kde(obs), color=c, s=100, zorder=5)

    ax.set_title(tlabel, fontsize=13, fontweight="bold")
    ax.set_xlabel("AED (mean squared divergence)", fontsize=11)
    ax.set_ylabel("Density" if t == 0 else "", fontsize=11)
    ax.legend(fontsize=8, loc="upper right")

plt.suptitle("Aggregated Expression Divergence — permutation test\\n"
             "Vertical lines = observed AED per genotype; shaded = null distribution (220 combinations)",
             fontsize=12, fontweight="bold")
plt.tight_layout()
plt.savefig("Fig2A_AED_density_plots.png", dpi=150, bbox_inches="tight")
plt.show()
"""
))

cells.append(md(
"""### Key findings from AED analysis

| Time | Significant genotypes (FDR < 0.05) |
|------|-------------------------------------|
| 0 hpi | *Rpv12+1*, *Rpv12+1+3* |
| 6 hpi | *Rpv12+1+3* only |
| 24 hpi | *Rpv12*, *Rpv12+1+3* |

**Interpretation:**
- *Rpv12* alone does not diverge significantly from susceptible at 0 hpi → no detectable pre-activation
- *Rpv12+1* and *Rpv12+1+3* show significant divergence at 0 hpi → altered basal transcriptional state
- Divergence is **not additive**: each genotype shows a different temporal trajectory
- The minimum attainable FDR (0.024) means borderline results (*Rpv12+1* at 24 hpi: FDR=0.068)  
  cannot be distinguished from true positives vs false negatives — requires quantitative resistance assays
"""
))

# ── SECTION 6 ────────────────────────────────────────────────────────────
cells.append(md(
"""## 8. Section 6 — Co-Transcriptional Network Analysis (PCNWA) <a name="sec6"></a>

> **Original R script:** `06_PCNWA/GCNA_network_analysis.R`

### Concept: from expression matrix to co-transcriptional network

If two genes respond synchronously across conditions (high Pearson correlation), they may be  
co-regulated — sharing a transcription factor, regulatory element, or signaling pathway.

The analysis steps are:
1. **Compute all-vs-all Pearson correlations** among 3 553 DEGs (manageable, ~12.6M pairs)
2. **Threshold selection**: retain only the top 0.5% of correlations (r > 0.817)  
   — empirical 99.5th percentile of the correlation distribution
3. **Build a co-expression graph**: nodes = genes, edges = strong correlations
4. **Module detection**: for each DE gene, find all genes with r ≥ 0.817;  
   use PCA eigengene to test if the module reflects the gene's DE pattern (Spearman, Bonferroni)
5. **Hub genes**: top 1% by degree centrality (number of co-expressed partners)
6. **Metamodules**: clusters of hub genes sharing many co-expressed partners

> ⚠️ **Computational note:** The full 26 169 × 26 169 correlation matrix requires ~5.6 GB RAM.  
> This was computed in Julia (GPU-accelerated) in the paper. Here we work with the 3 553 DE genes  
> (~100 MB for the correlation matrix), which is feasible in Python.
"""
))

cells.append(code(
"""# ── LOAD DE GENES AND EXPRESSION DATA ────────────────────────────────────────
# 3 553 genes with |log2FC| > 1 and FDR < 0.05 in at least one comparison
de_patterns = pd.read_csv(BASE / "Patterns_01_DUE.csv", sep="\\t", index_col=0)
de_genes     = de_patterns.index.tolist()
print(f"DE genes: {len(de_genes)}")
print("Pattern columns (Cult1=Rpv12, Cult2=Rpv12+1, Cult3=Rpv12+1+3):")
print(de_patterns.head(4))

# Subset VST + ComBat matrix to DE genes
expr_de = vst_combat.loc[de_genes, :]
print(f"\\nDE gene expression matrix: {expr_de.shape} (genes × samples)")
"""
))

cells.append(code(
"""# ── PEARSON CORRELATION MATRIX (3553 × 3553) ─────────────────────────────────
# R equivalent: correlation_matrix_byGenes <- cor(t(CBvst), use="pairwise.complete.obs")
# This is feasible: 3553^2 = ~12.6M entries ≈ 100 MB in float64

print("Computing Pearson correlation matrix for 3 553 DE genes...")
print("(This may take 1–3 minutes depending on hardware)")

# numpy corrcoef is fast: computes all-vs-all correlations
# Input must be (n_genes × n_samples) — shape (3553, 36)
cor_mat = np.corrcoef(expr_de.values)    # (3553, 3553) symmetric matrix

print(f"Correlation matrix shape: {cor_mat.shape}")
print(f"Value range: [{cor_mat.min():.3f}, {cor_mat.max():.3f}]")
print(f"Upper triangle (excluding diagonal): {3553*3552//2:,} pairs")
"""
))

cells.append(code(
"""# ── THRESHOLD SELECTION: 99.5TH PERCENTILE ───────────────────────────────────
# R equivalent:
#   length(unlist(correlation_matrix_byGenes)) # 684816561
#   sort(unlist(correlation_matrix_byGenes))[684816561-3424083] # 0.817
#
# For 3553 genes:
#   total pairs (upper triangle) = 3553*3552/2 = 6,310,878
#   top 0.5% threshold: 6,310,878 * 0.005 = 31,554 edges retained

upper_tri = cor_mat[np.triu_indices(len(de_genes), k=1)]   # upper triangle only
threshold = np.percentile(upper_tri, 99.5)

print(f"Total unique gene pairs: {len(upper_tri):,}")
print(f"99.5th percentile threshold: r = {threshold:.4f}")
print(f"Pairs above threshold: {(upper_tri >= threshold).sum():,} ({(upper_tri >= threshold).mean()*100:.2f}%)")
print(f"\\nPaper used r > 0.817 for the full 26 169-gene matrix.")
print(f"For 3 553 DE genes, the 99.5th percentile is r = {threshold:.3f}")

# Histogram of all pairwise correlations
fig, ax = plt.subplots(figsize=(9, 4))
ax.hist(upper_tri, bins=100, color="steelblue", edgecolor="white", alpha=0.8)
ax.axvline(threshold, color="darkorange", lw=2, ls="--",
           label=f"Cut-off = {threshold:.3f} (top 0.5%)")
ax.set_xlabel("Pearson r between gene pairs")
ax.set_ylabel("Count")
ax.set_title("Distribution of pairwise Pearson correlations (3 553 DE genes)")
ax.legend()
plt.tight_layout()
plt.savefig("Supp_Fig6_correlation_distribution.png", dpi=150, bbox_inches="tight")
plt.show()
"""
))

cells.append(code(
"""# ── BUILD NETWORKS PER GENOTYPE AND PER TIME POINT ───────────────────────────
# R equivalent: build_graph() function in GCNA_network_analysis.R

def build_correlation_graph(expr_subset, threshold):
    \"\"\"
    Build an undirected weighted graph where edges connect genes with Pearson r >= threshold.

    Parameters
    ----------
    expr_subset : pd.DataFrame (genes × samples for a subset of samples)
    threshold   : float, minimum correlation for an edge

    Returns
    -------
    nx.Graph with edge attribute 'weight' = Pearson r
    \"\"\"
    cm = np.corrcoef(expr_subset.values)
    np.fill_diagonal(cm, 0)   # no self-loops
    rows, cols = np.where(cm >= threshold)
    # Keep only upper triangle to avoid duplicate edges
    mask = rows < cols
    rows, cols = rows[mask], cols[mask]
    weights = cm[rows, cols]
    genes = expr_subset.index.tolist()
    G = nx.Graph()
    G.add_nodes_from(genes)
    G.add_weighted_edges_from([(genes[r], genes[c], w) for r, c, w in zip(rows, cols, weights)])
    return G

def network_metrics(G):
    \"\"\"Compute key topological metrics for a co-expression graph.\"\"\"
    n_edges = G.number_of_edges()
    n_nodes = G.number_of_nodes()
    density = nx.density(G)
    degrees = dict(G.degree())
    mean_degree = np.mean(list(degrees.values()))
    clustering = nx.transitivity(G) if n_edges > 0 else np.nan
    return {"n_nodes": n_nodes, "n_edges": n_edges, "density": density,
            "mean_degree": mean_degree, "clustering": clustering}

# Build genotype-specific networks (9 samples per genotype: 3 times × 3 replicates)
genotype_cols = {
    "Susceptible": [COND_IDX["Susceptible.0"] + COND_IDX["Susceptible.6"] + COND_IDX["Susceptible.24"]],
    "Rpv12":       [COND_IDX["Rpv12.0"] + COND_IDX["Rpv12.6"] + COND_IDX["Rpv12.24"]],
    "Rpv12+1":     [COND_IDX["Rpv12.1.0"] + COND_IDX["Rpv12.1.6"] + COND_IDX["Rpv12.1.24"]],
    "Rpv12+1+3":   [COND_IDX["Rpv12.1.3.0"] + COND_IDX["Rpv12.1.3.6"] + COND_IDX["Rpv12.1.3.24"]],
}
# Flatten lists
genotype_cols = {k: v[0] for k, v in genotype_cols.items()}

print("Building genotype-specific networks (threshold = paper's 0.817)...")
THR = 0.817   # use paper's threshold for comparability
graphs_geno = {}
metrics_geno = []

for geno, cols in genotype_cols.items():
    G = build_correlation_graph(expr_de.iloc[:, cols], THR)
    graphs_geno[geno] = G
    m = network_metrics(G)
    m["genotype"] = geno
    metrics_geno.append(m)

metrics_geno_df = pd.DataFrame(metrics_geno).set_index("genotype")
print("\\nNetwork metrics by genotype (Pearson r > 0.817):")
print(metrics_geno_df.to_string(float_format="%.4f"))
print("\\nExpected from paper:")
print("  Susceptible: density=0.101, mean_degree=357.7")
print("  Rpv12:       density=0.132, mean_degree=468.5")
print("  Rpv12+1:     density=0.133, mean_degree=471.3")
print("  Rpv12+1+3:   density=0.055, mean_degree=194.1  ← notably lower connectivity")
"""
))

cells.append(code(
"""# ── HUB GENE IDENTIFICATION ───────────────────────────────────────────────────
# Build a global network (all 36 samples) to identify hub genes.
# Hubs = genes in top 1% of degree centrality.
# R equivalent: threshold <- quantile(degree_centrality, probs = 0.99)

print("Building global network (all 36 samples) for hub gene identification...")
G_global = build_correlation_graph(expr_de, THR)

degree_centrality = dict(G_global.degree())
dc_values = np.array(list(degree_centrality.values()))
dc_genes  = np.array(list(degree_centrality.keys()))

hub_threshold = np.quantile(dc_values, 0.99)
hub_genes = dc_genes[dc_values >= hub_threshold].tolist()

print(f"\\nDegree centrality distribution:")
print(f"  Total nodes: {len(dc_values)}")
print(f"  Min: {dc_values.min()}, Max: {dc_values.max()}, Median: {np.median(dc_values):.0f}")
print(f"  99th percentile threshold: {hub_threshold:.0f}")
print(f"  Hub genes (top 1%): {len(hub_genes)}")

# Distribution plot
fig, ax = plt.subplots(figsize=(8, 4))
ax.hist(np.log10(dc_values + 1), bins=50, color="steelblue", edgecolor="white")
ax.axvline(np.log10(hub_threshold), color="firebrick", lw=2, ls="--",
           label=f"Hub threshold (top 1%) = {hub_threshold:.0f}")
ax.set_xlabel("log₁₀(degree centrality)")
ax.set_ylabel("Count")
ax.set_title("Degree centrality distribution — co-transcriptional network (3 553 DE genes)")
ax.legend()
plt.tight_layout()
plt.savefig("Supp_degree_centrality.png", dpi=150, bbox_inches="tight")
plt.show()

print(f"\\nHub genes (first 10):")
hub_dc = sorted([(g, degree_centrality[g]) for g in hub_genes], key=lambda x: -x[1])
for g, d in hub_dc[:10]:
    print(f"  {g}: {d} partners")
"""
))

cells.append(code(
"""# ── METAMODULE OVERLAP ANALYSIS ──────────────────────────────────────────────
# Metamodules = groups of hub genes whose co-expressed gene sets substantially overlap.
# Motivation: overlapping modules likely share regulatory programs.
# R equivalent: nested loop computing intersection of module gene sets for each pair of hub genes.

print(f"Computing pairwise module overlap for {len(hub_genes)} hub genes...")
print("(This may take a few minutes for large hub gene sets)")

n_hubs = len(hub_genes)
shared_counts = np.zeros((n_hubs, n_hubs), dtype=int)
module_sizes  = {}

for i, hub1 in enumerate(hub_genes):
    # Genes correlated with hub1 at r >= threshold (excluding itself)
    hub1_idx = de_genes.index(hub1) if hub1 in de_genes else -1
    if hub1_idx < 0:
        continue
    mod1 = set(dc_genes[cor_mat[hub1_idx, :] >= THR]) - {hub1}
    module_sizes[hub1] = len(mod1)
    for j, hub2 in enumerate(hub_genes):
        hub2_idx = de_genes.index(hub2) if hub2 in de_genes else -1
        if hub2_idx < 0:
            continue
        mod2 = set(dc_genes[cor_mat[hub2_idx, :] >= THR]) - {hub2}
        shared_counts[i, j] = len(mod1.intersection(mod2))

# Convert to proportion of the smaller module
shared_perc = np.zeros_like(shared_counts, dtype=float)
for i, hub1 in enumerate(hub_genes):
    for j, hub2 in enumerate(hub_genes):
        min_size = min(module_sizes.get(hub1, 1), module_sizes.get(hub2, 1))
        if min_size > 0:
            shared_perc[i, j] = shared_counts[i, j] / min_size

shared_perc_df = pd.DataFrame(shared_perc, index=hub_genes, columns=hub_genes)
print(f"Metamodule overlap matrix: {shared_perc_df.shape}")
"""
))

cells.append(code(
"""# ── HIERARCHICAL CLUSTERING HEATMAP OF METAMODULE OVERLAP ───────────────────
# Three fixes vs naive clustermap:
#   1. Zero the diagonal (self-overlap = 1.0 trivially) so it doesn't compress the scale
#   2. Set vmax from the 95th percentile of off-diagonal values → full colour range used
#   3. Save with g.fig.savefig() — plt.savefig() would save the wrong figure

shared_perc_plot = shared_perc_df.copy()
np.fill_diagonal(shared_perc_plot.values, 0)   # diagonal is uninformative

# Colour scale: 0 → 95th percentile of non-zero off-diagonal overlaps
off_diag = shared_perc_plot.values[np.triu_indices_from(shared_perc_plot.values, k=1)]
vmax_val = float(np.percentile(off_diag[off_diag > 0], 95)) if (off_diag > 0).any() else 1.0
vmax_val = round(vmax_val, 2)
print(f"Colour scale: 0 – {vmax_val:.2f}  (95th percentile of off-diagonal overlaps)")

g = sns.clustermap(
    shared_perc_plot,
    method="ward",
    metric="euclidean",
    cmap="YlOrRd",
    vmin=0, vmax=vmax_val,
    xticklabels=False, yticklabels=True,   # x-labels off: too many hub genes to read
    figsize=(14, 12),
    cbar_kws={"label": f"Proportion of shared genes  (max = {vmax_val:.2f})"},
)
# Title on the figure (not the heatmap axes) to avoid overlap with the dendrogram
g.fig.suptitle(
    "Hub gene co-module overlap — Metamodule structure\\n(Supplementary Figure 8 equivalent)",
    fontsize=12, y=1.01
)
g.fig.savefig("Supp_Fig8_metamodule_overlap.png", dpi=150, bbox_inches="tight")
plt.show()
print("Metamodule heatmap saved.")
"""
))

cells.append(code(
"""# ── PTI-ETI CROSS-TALK ANALYSIS ───────────────────────────────────────────────
# Test whether PTI and ETI genes are more strongly co-expressed than expected by chance.
# Immunity class file: each gene labelled PTI_specific, ETI_specific, PTI_ETI_shared, or Other
gene_classif = pd.read_csv(BASE / "Vitis_gene_immunity_class.csv", index_col=0)

print("Immunity class distribution (all genes):")
print(gene_classif["Immunity_Class"].value_counts())

immune_genes = gene_classif[gene_classif["Immunity_Class"] != "Other"].index.tolist()
print(f"\\nImmunity-annotated genes in the DE set: {len(set(immune_genes) & set(de_genes))}")

# Cross-talk: do PTI and ETI genes tend to co-express?
# Under H0: edges between PTI and ETI genes should occur at random (proportional to their frequencies)
def cross_talk_test(G, gene_classes_series):
    \"\"\"
    Compute observed vs expected proportion of PTI-ETI cross-class edges.

    Returns DataFrame with n_edges, cross_edges, cross_prop, expected_prop, odds_ratio.
    \"\"\"
    # Annotate nodes
    valid = {n: gene_classes_series[n] for n in G.nodes() if n in gene_classes_series.index}
    G_sub = G.subgraph(list(valid.keys()))
    classes = pd.Series(valid)

    if G_sub.number_of_edges() == 0:
        return pd.Series({"n_edges": 0, "cross_edges": 0, "cross_prop": np.nan,
                          "expected_prop": np.nan, "odds_ratio": np.nan})

    cross_classes = {"PTI_specific", "ETI_specific", "PTI_ETI_shared"}

    cross_edges = 0
    for u, v in G_sub.edges():
        cu, cv = classes.get(u), classes.get(v)
        if cu in cross_classes and cv in cross_classes and cu != cv:
            cross_edges += 1

    n_edges = G_sub.number_of_edges()
    cross_prop = cross_edges / n_edges

    # Expected under random pairing
    freqs = classes.value_counts(normalize=True)
    p_pti  = freqs.get("PTI_specific", 0)
    p_eti  = freqs.get("ETI_specific", 0)
    p_sh   = freqs.get("PTI_ETI_shared", 0)
    exp_prop = p_pti*p_eti + p_pti*p_sh + p_eti*p_sh

    or_ = cross_prop / exp_prop if exp_prop > 0 else np.nan
    return pd.Series({"n_edges": n_edges, "cross_edges": cross_edges,
                      "cross_prop": cross_prop, "expected_prop": exp_prop, "odds_ratio": or_})

gene_classes_de = gene_classif.loc[gene_classif.index.isin(de_genes), "Immunity_Class"]

print("\\nPTI–ETI cross-talk analysis (co-transcriptional networks):")
for geno, G in graphs_geno.items():
    result = cross_talk_test(G, gene_classes_de)
    print(f"  {geno:12s}: {result['cross_edges']:3.0f} cross-class edges / {result['n_edges']:,} total "
          f"(observed={result['cross_prop']:.6f}, expected={result['expected_prop']:.6f}, "
          f"OR={result['odds_ratio']:.2f})")
print("\\nNote: At r > 0.817, PTI and ETI genes do not directly co-express at the pairwise level.")
print("This is expected: PTI and ETI act in staggered temporal waves, not perfect co-expression.")
print("Cross-talk detection requires looking at NETWORK MODULE OVERLAP (metamodules), not edge counts.")
"""
))

# ── SECTION 7 ────────────────────────────────────────────────────────────
cells.append(md(
"""## 9. Section 7 — Key Immunity Gene Network (Figure 5 equivalent) <a name="sec7"></a>

> **Original R script:** `06_PCNWA/Figure_5.R`

### Focus: the EDS1–SAG101–SOBIR1 immune hub

Among the 56 hub genes identified in Section 6, nine encode well-characterised immunity components:

| Gene ID | Symbol | Layer | Function |
|---------|--------|-------|----------|
| Vitvi17g04216 | EDS1 | Signal Integration | Lipase-domain signaling hub; links TNL-ETI to PTI |
| Vitvi14g03033 | SAG101 | Signal Integration | EDS1 partner in TNL-mediated cell death |
| Vitvi14g03030 | SAG101-like | Signal Integration | EDS1 interacting paralog |
| Vitvi14g03031 | SAG101 | Signal Integration | EDS1 interacting paralog |
| Vitvi05g00062 | SAG101 | Signal Integration | EDS1 interacting paralog |
| Vitvi16g01127 | SARD4 | Defense Action | Pipecolic acid biosynthesis; SAR priming |
| Vitvi13g00189 | WRKY55 | Signal Integration | Immune transcription factor |
| Vitvi07g01847 | WRKY51 | Signal Integration | Immune transcription factor |
| Vitvi14g03057 | MIK2 | Recognition | LRR-RLK co-receptor, pattern recognition |
| Vitvi17g00964 | SOBIR1 | Recognition | RLK co-receptor; hub of predicted PPI network |
| Vitvi12g02607 | RUN1/ACD6 | Recognition/Defense | Resistance gene / cell death regulator |
| Vitvi05g00637 | ACD6 | Defense Action | Accelerated cell death |
| Vitvi09g04475 | EIX1 | Recognition | LRR-RLP; predicted SOBIR1 partner |
| Vitvi09g04536 | GSO2 | Recognition | LRR-RLK; predicted SOBIR1 partner |
| Vitvi03g00614 | EXLRR11 | Recognition | LRR protein; predicted SARD4 partner |

All nine core genes were **upregulated at 0 hpi specifically in *Rpv12*** — not in multilocus  
backgrounds. This Rpv12-specific early activation is a key finding of the paper.
"""
))

cells.append(code(
"""# ── IMMUNITY GENES: LOG2FC HEATMAP ───────────────────────────────────────────
# R equivalent: Figure_5.R Panel A

immunity_genes = {
    "Vitvi17g04216": "EDS1",
    "Vitvi14g03033": "SAG101-like†",
    "Vitvi14g03030": "SAG101-a",
    "Vitvi14g03031": "SAG101-b",
    "Vitvi05g00062":  "SAG101-c",
    "Vitvi16g01127": "SARD4",
    "Vitvi13g00189": "WRKY55",
    "Vitvi07g01847": "WRKY51",
    "Vitvi14g03057": "MIK2",
    "Vitvi17g00964": "SOBIR1",
    "Vitvi12g02607": "RUN1/ACD6",
    "Vitvi05g00637": "ACD6",
    "Vitvi09g04475": "EIX1",
    "Vitvi09g04536": "GSO2",
    "Vitvi03g00614": "EXLRR11",
}

# Subset to genes present in the expression matrix
avail = {k: v for k, v in immunity_genes.items() if k in vst_combat.index}
print(f"Immunity genes found in expression matrix: {len(avail)}/{len(immunity_genes)}")

sub = sf_combat.loc[list(avail.keys()), :]
sub.index = [avail[g] for g in sub.index]

# Compute mean per condition (replicates A, B, C)
cond_order = [
    ("Rpv12.0",     "Susceptible.0",   "Rpv12_0"),
    ("Rpv12.6",     "Susceptible.6",   "Rpv12_6"),
    ("Rpv12.24",    "Susceptible.24",  "Rpv12_24"),
    ("Rpv12.1.0",   "Susceptible.0",   "Rpv12+1_0"),
    ("Rpv12.1.6",   "Susceptible.6",   "Rpv12+1_6"),
    ("Rpv12.1.24",  "Susceptible.24",  "Rpv12+1_24"),
    ("Rpv12.1.3.0", "Susceptible.0",   "Rpv12+1+3_0"),
    ("Rpv12.1.3.6", "Susceptible.6",   "Rpv12+1+3_6"),
    ("Rpv12.1.3.24","Susceptible.24",  "Rpv12+1+3_24"),
]

heatmap_data = {}
for cond_a, cond_b, label in cond_order:
    mean_a = sub.iloc[:, COND_IDX[cond_a]].mean(axis=1)
    mean_b = sub.iloc[:, COND_IDX[cond_b]].mean(axis=1)
    heatmap_data[label] = mean_a - mean_b   # log2FC (on log2-scale values)

heatmap_df = pd.DataFrame(heatmap_data, index=sub.index)
# Reorder columns for the plot
col_order = [c[2] for c in cond_order]
heatmap_df = heatmap_df[col_order]
print("\\nLog2FC heatmap data (first 5 rows):")
print(heatmap_df.head())
"""
))

cells.append(code(
"""# Panel A: Log2FC heatmap
fig, ax = plt.subplots(figsize=(11, 7))
sns.heatmap(
    heatmap_df,
    ax=ax,
    cmap=sns.diverging_palette(240, 10, n=100),
    center=0, vmin=-3, vmax=3,
    linewidths=0.5, linecolor="lightgray",
    cbar_kws={"label": "log₂ FCH vs susceptible", "shrink": 0.7},
    annot=False,
)
# Add column group separators
for x in [3, 6]:
    ax.axvline(x, color="black", lw=1.5)

ax.set_title("Log₂ fold-change of immunity genes vs susceptible\\n(Figure 5A equivalent)",
             fontsize=12, fontweight="bold")
ax.set_xlabel("")
ax.set_ylabel("Gene", fontsize=11)
plt.xticks(
    ticks=np.arange(len(col_order)) + 0.5,
    labels=[c.replace("_", " hpi, ") for c in col_order],
    rotation=45, ha="right", fontsize=9
)
plt.tight_layout()
plt.savefig("Fig5A_immunity_heatmap.png", dpi=150, bbox_inches="tight")
plt.show()
"""
))

cells.append(code(
"""# Panel B: Pearson correlation heatmap among immunity genes
# R equivalent: cor_mat <- cor(t(sub), method="pearson", use="pairwise.complete.obs")
sub_vst = vst_combat.loc[[k for k in avail.keys() if k in vst_combat.index], :]
sub_vst.index = [avail[g] for g in sub_vst.index]

cor_immunity = sub_vst.T.corr(method="pearson")

g_clust = sns.clustermap(
    cor_immunity,
    method="ward", metric="euclidean",
    cmap="coolwarm", center=0, vmin=-1, vmax=1,
    linewidths=0.5, linecolor="lightgray",
    figsize=(10, 9),
    cbar_kws={"label": "Pearson r", "shrink": 0.5},
)
g_clust.ax_heatmap.set_title("Pairwise Pearson correlations among immunity genes\\n"
                              "(Figure 5B equivalent)", fontsize=11)
plt.savefig("Fig5B_immunity_correlation.png", dpi=150, bbox_inches="tight")
plt.show()
"""
))

cells.append(code(
"""# Panel C: Co-transcriptional network (Pearson r ≥ 0.65)
# R equivalent: Figure_5.R Panel C using ggraph
# Nodes coloured by functional layer; edge width ∝ Pearson r

# Functional layer annotations
layer_map = {
    "EDS1":         "Signal Integration",
    "SAG101-like†": "Signal Integration",
    "SAG101-a":     "Signal Integration",
    "SAG101-b":     "Signal Integration",
    "SAG101-c":     "Signal Integration",
    "SARD4":        "Defense Action",
    "WRKY55":       "Signal Integration",
    "WRKY51":       "Signal Integration",
    "MIK2":         "Recognition",
    "SOBIR1":       "Recognition",
    "RUN1/ACD6":    "Recognition",
    "ACD6":         "Defense Action",
    "EIX1":         "Recognition",
    "GSO2":         "Recognition",
    "EXLRR11":      "Recognition",
}
layer_colors = {
    "Recognition":        "orchid",
    "Signal Integration": "lightblue",
    "Defense Action":     "gold",
}
EDGE_THR = 0.65   # softer threshold for the focused network (paper used 0.65)

G_imm = nx.Graph()
nodes = sub_vst.index.tolist()
G_imm.add_nodes_from(nodes)

for i, g1 in enumerate(nodes):
    for j, g2 in enumerate(nodes):
        if j <= i: continue
        r = cor_immunity.loc[g1, g2]
        if r >= EDGE_THR:
            G_imm.add_edge(g1, g2, weight=r)

np.random.seed(123)
pos = nx.spring_layout(G_imm, weight="weight", k=2, iterations=100, seed=123)

fig, ax = plt.subplots(figsize=(11, 9))
# Draw edges with width proportional to correlation
edges = G_imm.edges(data=True)
edge_weights = [d["weight"] for _, _, d in edges]
nx.draw_networkx_edges(G_imm, pos, ax=ax, width=[w*3 for w in edge_weights],
                       alpha=0.6, edge_color="dimgray")
# Draw nodes
for node in G_imm.nodes():
    layer = layer_map.get(node, "Other")
    c = layer_colors.get(layer, "lightgray")
    nx.draw_networkx_nodes(G_imm, pos, nodelist=[node], node_color=c,
                           node_size=800, ax=ax, edgecolors="gray", linewidths=1)
nx.draw_networkx_labels(G_imm, pos, ax=ax, font_size=8, font_weight="bold")
# Legend
layer_handles = [mpatches.Patch(color=v, label=k) for k, v in layer_colors.items()]
ax.legend(handles=layer_handles, title="Functional layer", loc="lower left", fontsize=9)
ax.set_title(f"Co-transcriptional network: immunity genes (r ≥ {EDGE_THR})\\n"
             "Edge width ∝ Pearson r  |  Figure 5C equivalent", fontsize=11)
ax.axis("off")
plt.tight_layout()
plt.savefig("Fig5C_immunity_network.png", dpi=150, bbox_inches="tight")
plt.show()
"""
))

cells.append(code(
"""# Panel D: STRING-db network (experimental + computational interaction evidence)
# Loads the pre-curated STRING interaction table
string_edges = pd.read_csv(BASE / "node1_node2_ConfidenceScore_stringDB.csv", sep="\\t")
print("STRING-db edges (first 5):")
print(string_edges.head())

G_str = nx.from_pandas_edgelist(string_edges, source="node1", target="node2",
                                edge_attr="score")
# Assign layers
for node in G_str.nodes():
    G_str.nodes[node]["layer"] = layer_map.get(node, "Other")

np.random.seed(123)
pos_str = nx.spring_layout(G_str, weight="score", k=2.5, seed=123)

fig, ax = plt.subplots(figsize=(10, 8))
edge_scores = [d["score"] for _, _, d in G_str.edges(data=True)]
nx.draw_networkx_edges(G_str, pos_str, ax=ax, width=[s*3 for s in edge_scores],
                       alpha=0.5, edge_color="gray60")
for node in G_str.nodes():
    layer = G_str.nodes[node]["layer"]
    c = layer_colors.get(layer, "lightgray")
    nx.draw_networkx_nodes(G_str, pos_str, nodelist=[node], node_color=c,
                           node_size=800, ax=ax, edgecolors="gray", linewidths=1)
nx.draw_networkx_labels(G_str, pos_str, ax=ax, font_size=9, font_weight="bold")
ax.legend(handles=layer_handles, title="Functional layer", loc="lower left", fontsize=9)
ax.set_title("STRING-db protein interaction network (immunity hub genes)\\n"
             "Edge width ∝ STRING confidence score  |  Figure 5D equivalent", fontsize=11)
ax.axis("off")
plt.tight_layout()
plt.savefig("Fig5D_STRING_network.png", dpi=150, bbox_inches="tight")
plt.show()
"""
))

cells.append(md(
"""## Conclusions

### What the pipeline revealed

1. **Batch effects are real and must be corrected before any biological analysis.**  
   GLMs alone (within DESeq2 design) do not remove batch signal from expression matrices.  
   ComBat is required.

2. **Two complementary normalizations are necessary:**
   - Size-factor + ComBat for amplitude-sensitive analyses (AED)
   - VST + ComBat for correlation-based analyses (PCA, network, DGEA on log-scale)

3. **Locus pyramiding does not produce additive transcriptional responses.**  
   Each genotype shows a distinct temporal trajectory:  
   - *Rpv12* activates SOBIR1-centered PRR networks early (0 hpi)  
   - *Rpv12+1* shows peak DEG counts at 6 hpi (early response dominates)  
   - *Rpv12+1+3* shows the largest transcriptome shift at 6 hpi but reduced connectivity

4. **The co-transcriptional network architecture changes with locus dosage.**  
   Rpv12+1+3 shows dramatically reduced network density (0.055 vs 0.132 for Rpv12),  
   suggesting extensive transcriptional restructuring rather than strengthened coordination.

5. **SOBIR1 emerges as a hub with Rpv12-specific upregulation at 0 hpi,**  
   forming a predicted receptor complex with 7 LRR-containing partners —  
   consistent with the tomato LeEIX-SOBIR1 architecture but not conserved in multilocus backgrounds.

### PyDESeq2 vs the paper's GLM approach

| Aspect | PyDESeq2 | Paper's GLM (VST+ComBat) |
|--------|----------|--------------------------|
| Model | Negative binomial (count-appropriate) | Gaussian (approximate) |
| Input | Raw integer counts | Batch-corrected VST values |
| Shrinkage | Yes (dispersion + LFC) | No |
| Speed | Slower (~5–15 min all 26169 genes) | Fast |
| Best for | High-confidence DE gene lists | Consistency with downstream correlation analyses |

For a definitive gene list, **PyDESeq2 is preferred**. The paper's approach was chosen for  
practical consistency with the downstream network analysis.

### What biological validation is needed

- Confirm SOBIR1 requirement: CRISPR knockout in Rpv12 background  
- Validate predicted complexes: co-immunoprecipitation  
- Quantify resistance: sporulation assays to link transcriptional patterns to actual pathogen restriction
"""
))

# ── BUILD NOTEBOOK ───────────────────────────────────────────────────────
nb = {
    "nbformat": 4,
    "nbformat_minor": 5,
    "metadata": {
        "kernelspec": {
            "display_name": "Python (plant_immunity)",
            "language": "python",
            "name": "plant_immunity"
        },
        "language_info": {
            "codemirror_mode": {"name": "ipython", "version": 3},
            "file_extension": ".py",
            "mimetype": "text/x-python",
            "name": "python",
            "pygments_lexer": "ipython3",
            "version": "3.10.0"
        }
    },
    "cells": cells,
}

out_path = pathlib.Path(__file__).with_name("plant_immunity_python_toolkit.ipynb")
with open(out_path, "w", encoding="utf-8") as f:
    json.dump(nb, f, indent=1, ensure_ascii=False)

print(f"Notebook written to: {out_path}")
print(f"Total cells: {len(cells)}")
md_cells   = sum(1 for c in cells if c["cell_type"] == "markdown")
code_cells = sum(1 for c in cells if c["cell_type"] == "code")
print(f"  Markdown cells: {md_cells}")
print(f"  Code cells:     {code_cells}")
