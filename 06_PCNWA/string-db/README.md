# string-db

`scripts/SF6_network_robustness_rebuild.R` depends on
`../combat_protected/results/full_gene_correlation_protected_ComBat_dist.rds`
(the cached correlation matrix) and the module-count CSVs in
`../combat_protected/correlation_cutoffs/`. Raw STRING data in `raw_data/`
needs no rebuilding — it's a static v12.0 download.

See `ROADMAP.md` for the STRING coverage caveat.
