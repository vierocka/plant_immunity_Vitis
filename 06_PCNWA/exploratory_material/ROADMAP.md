# ROADMAP: exploratory_material

**AIM.** Hold superseded network builds and pre-rebuild figure drafts, out
of the active folders, without deleting them.

**MOTIVATION.** Nothing here is read by any current script or cited in the
manuscript.

## Contents and why each item is here
- **`canonical_DESeq2_GCNA/` + `recompute_GCNA_from_canonical_DESeq2.R`** —
  built on a hardcoded r=0.817, superseded by the data-driven r=0.8077
  rebuild in `../combat_protected/`. Nothing downstream reads this folder
  (checked directly).
- **`GCNA_module_reconstruction/` + `GCNA_module_reconstruction_both_ComBat.R`** —
  anchored on the historical 3,553-gene t-test panel; superseded for
  current purposes by the canonical 9,459-gene panel used throughout
  `../combat_protected/` and `../network_robustness/`.
- **`SF7_network_density_by_time_r0.8077.R`** — the original 1-panel
  Supplementary Figure 7, superseded by the 2-panel
  `../combat_protected/scripts/SF7_rebuild.R`.
- **Pre-rebuild Figure 4/5 material** (`Figure_4_new.*`, `Figure_4_panel*.png`,
  `Figure_5.R`, `Figure_5.jpeg`, `Figure_5.tiff`, `Figure5.tiff`,
  `Figure_5_new*.R/.png/.tiff`) — earlier script/figure generations (April–
  November, predating the August 2026 `_rebuild` versions in
  `../combat_protected/`).
- **`155modules_size_distribution.jpg`, `ChangedExpression_3553genes_fullInfo.csv`** —
  pre-canonical (2024–2025), historical-panel-based.
