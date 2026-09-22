# combat_protected

Run order: `scripts/save_full_correlation_matrix_compact.R` (cached matrix,
needed by several others) -> `scripts/GCNA_rebuild_9459panel_FDR_datadriven_rcutoff.R`
(primary module build) -> `scripts/hub_genes_6600_r0.8077_FDR.R` ->
`scripts/create_immunity_classes.R` -> `scripts/Figure3_pheatmap_MMs.R`,
`Figure_4_rebuild_18genes.R`, `Figure_5_rebuild.R`, `SF7_rebuild.R`
(figures, can run independently once the module build exists).

See `ROADMAP.md` for what each piece is and why.
