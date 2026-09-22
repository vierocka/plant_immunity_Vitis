###############################################################################
# AIM: save the full 26,169 x 26,169 gene-gene Pearson correlation matrix
# (protected ComBat) in the smallest practical, R-reusable form.
# MOTIVATION: a dense matrix would be ~5.48 GB, symmetric with a trivial
# diagonal; storing only the lower triangle (class "dist") + xz-compression
# shrinks this substantially.
# Reuse: d <- readRDS(...); cor_mat <- as.matrix(d); diag(cor_mat) <- 1
# (as.matrix(dist) zeroes the diagonal, fix it back to 1).
###############################################################################

expr_df <- read.delim("data_files/Rlogs_ComBat_protected.csv", check.names = FALSE)
expr_mat <- as.matrix(expr_df[, -1, drop = FALSE])
storage.mode(expr_mat) <- "double"
rownames(expr_mat) <- as.character(expr_df[[1]])
stopifnot(!anyNA(expr_mat), ncol(expr_mat) == 36, nrow(expr_mat) == 26169)

message("Computing full correlation matrix...")
t0 <- Sys.time()
full_cor <- cor(t(expr_mat), use = "everything", method = "pearson")
message("Done in ", round(as.numeric(Sys.time() - t0, units = "secs"), 1), "s.")

d <- as.dist(full_cor)   # lower triangle only, class "dist", attr(d,"Labels") = gene IDs
rm(full_cor); invisible(gc())

out_path <- "06_PCNWA/full_gene_correlation_protected_ComBat_dist.rds"
message("Writing (xz-compressed) to ", out_path, " ...")
t1 <- Sys.time()
saveRDS(d, out_path, compress = "xz")
message("Done in ", round(as.numeric(Sys.time() - t1, units = "secs"), 1), "s. File size: ",
        round(file.size(out_path) / 1e9, 2), " GB.")
