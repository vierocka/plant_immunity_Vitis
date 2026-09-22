###############################################################################
# Updated Supplementary Figure 7 (2026-08-28): network density across time,
# rebuilt on the 9,459-gene canonical-DESeq2 DE panel, rlog + condition-
# protected ComBat, r >= 0.8077 (99.5th percentile) -- replaces the old
# 3,553-gene / r>0.817 version. Values read directly from
# network_density_by_timepoint_r0.8077_protected_9459panel.csv (already
# verified against the live rebuild -- see help.txt).
###############################################################################

suppressPackageStartupMessages(library(ggplot2))

out_dir <- "../../draft/Sections_by_VK/MPMI/Supplementary_figures"

dens <- read.csv("network_density_by_timepoint_r0.8077_protected_9459panel.csv",
                  stringsAsFactors = FALSE)
dens$hpi <- c(0, 6, 24)

p <- ggplot(dens, aes(x = hpi, y = density)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 4) +
  scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25)) +
  labs(x = "Time post infection (hpi)",
       y = expression("Network density (" * r >= 0.8077 * ")")) +
  theme_bw(base_size = 14)

ggsave(file.path(out_dir, "Supplementary_Figure_7_PROTECTED.png"), p, width = 8, height = 6.8, dpi = 300)
ggsave(file.path(out_dir, "Supplementary_Figure_7_PROTECTED.pdf"), p, width = 8, height = 6.8)

cat("Density values plotted:\n")
print(dens[, c("hpi", "density")])
cat("\nDone.\n")
