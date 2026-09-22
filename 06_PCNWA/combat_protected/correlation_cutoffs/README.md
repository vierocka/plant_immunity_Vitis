# correlation_cutoffs

Run order: `GCNA_rebuild_9459panel_r0.757.R` and `GCNA_rebuild_9459panel_r0.8839.R`
independently (each rebuilds the module set at its own threshold) ->
`hub_genes_1283_r0.8839_bonf0.001.R`, `network_density_recompute_r0.8839.R`.
Both depend on `../scripts/save_full_correlation_matrix_compact.R`'s cached
correlation matrix having already been built.

See `ROADMAP.md` for what each threshold represents.
