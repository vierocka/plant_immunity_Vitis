############ DATA LOADING AND PREPARATION ####################
myData <- read.table("data_files/Rlogs.csv", header = TRUE, sep="\t")
CBrlogs <- as.matrix(myData[,c(2:37)])

myData2 <- read.csv("data_files/SizeFactorNorm_ComBat_36samples.csv")
SFcombat <- as.matrix(myData2[,c(2:37)])

rownames(CBrlogs) <- myData[,1]
rownames(SFcombat) <- myData2[,2]

# define vectors for conditions (info: introduced loci and time points)
condition <- as.factor(rep(c("Rpv12.0","Rpv12.6","Rpv12.24","Rpv12.1.0","Rpv12.1.6","Rpv12.1.24","Rpv12.1.3.0","Rpv12.1.3.6","Rpv12.1.3.24", "Susceptible.0","Susceptible.6","Susceptible.24"),3))
# only time point info
myTime <- as.factor(rep(c(0,6,24),12))
# resistance genotype
myResistance <- as.factor(rep(c(rep("Rpv1",3), rep("Rpv1.12",3), rep("Rpv1.12.3",3), rep("Susceptible",3)),3))

######################################### HETEROSCEDASCITY #######################################

############ RLOG + COMBAT ##############

######################### SPEARMAN ##############
## assumption: mean–variance relationships in expression data are nonlinear
## Does variance increase systematically with mean expression?
## Spearman!! because Spearman correlation is rank-based, so it measures monotonic trends

check_heteroscedasticity <- function(expr_mat, condition, plot = TRUE) {
  # expr_mat: numeric matrix, genes in rows, samples in columns (e.g. rlog + ComBat)
  # condition: factor or vector with same length as ncol(expr_mat)
  # plot: whether to draw diagnostic plots
  
  stopifnot(ncol(expr_mat) == length(condition))
  condition <- as.factor(condition)
  
  cat("=== GLOBAL MEAN–VARIANCE RELATIONSHIP ===\n")
  gene_means <- rowMeans(expr_mat)
  gene_vars  <- apply(expr_mat, 1, var)
  
  cor_global <- cor(gene_means, gene_vars, method = "spearman")
  lm_global  <- lm(gene_vars ~ gene_means)
  slope_global <- coef(lm_global)[2]
  
  cat(sprintf("Spearman r = %.3f | slope = %.4f\n", cor_global, slope_global))
  if (plot) {
    plot(gene_means, gene_vars, pch = 16, cex = 0.4,
         xlab = "Mean expression (rlog scale)",
         ylab = "Variance (rlog scale)",
         main = sprintf("Global mean–variance (r=%.2f)", cor_global))
    abline(lm_global, col = "red", lwd = 2)
  }
  
  cat("\n=== PER-CONDITION RELATIONSHIPS ===\n")
  results <- data.frame(
    condition = levels(condition),
    Spearman_r = NA,
    slope = NA,
    stringsAsFactors = FALSE
  )
  
  for (cond in levels(condition)) {
    idx <- which(condition == cond)
    if (length(idx) < 2) next  # skip groups with <2 samples
    
    gene_means <- rowMeans(expr_mat[, idx, drop = FALSE])
    gene_vars  <- apply(expr_mat[, idx, drop = FALSE], 1, var)
    
    r <- cor(gene_means, gene_vars, method = "spearman")
    s <- coef(lm(gene_vars ~ gene_means))[2]
    results[results$condition == cond, c("Spearman_r", "slope")] <- c(r, s)
    
    cat(sprintf("%-15s | r = %.3f | slope = %.4f\n", cond, r, s))
    
    if (plot) {
      plot(gene_means, gene_vars, pch = 16, cex = 0.4,
           xlab = sprintf("Mean (%s)", cond),
           ylab = "Variance",
           main = sprintf("%s (r=%.2f)", cond, r))
      abline(lm(gene_vars ~ gene_means), col = "red")
    }
  }
  
  invisible(list(global = list(r = cor_global, slope = slope_global),
                 per_condition = results))
}

check_heteroscedasticity(CBrlogs, condition)
# === GLOBAL MEAN–VARIANCE RELATIONSHIP ===
#  Spearman r = -0.030 | slope = -0.0014

# === PER-CONDITION RELATIONSHIPS ===
# Rpv12.0         | r = -0.053 | slope = -0.0011
# Rpv12.1.0       | r = -0.149 | slope = -0.0026
# Rpv12.1.24      | r = -0.008 | slope = -0.0000
# Rpv12.1.3.0     | r = -0.058 | slope = -0.0024
# Rpv12.1.3.24    | r = -0.090 | slope = -0.0016
# Rpv12.1.3.6     | r = 0.007 | slope = -0.0027
# Rpv12.1.6       | r = -0.047 | slope = -0.0014
# Rpv12.24        | r = -0.250* | slope = -0.0018 
# Rpv12.6         | r = 0.021 | slope = 0.0002
# Susceptible.0   | r = -0.080 | slope = -0.0014
# Susceptible.24  | r = -0.007 | slope = -0.0003
# Susceptible.6   | r = -0.388* | slope = -0.0022

check_heteroscedasticity(CBrlogs, myTime)
# === PER-CONDITION RELATIONSHIPS ===
# 0               | r = -0.031 | slope = 0.0013
# 6               | r = -0.064 | slope = -0.0035
# 24              | r = -0.036 | slope = -0.0034

check_heteroscedasticity(CBrlogs, myResistance)
# === PER-CONDITION RELATIONSHIPS ===
# Rpv1            | r = 0.019 | slope = 0.0018
# Rpv1.12         | r = -0.020 | slope = 0.0012
# Rpv1.12.3       | r = -0.079 | slope = -0.0023
# Susceptible     | r = -0.057 | slope = -0.0012

######################### PEARSON ##############
## Does variance increase linearly with mean?
## Pearson checks residual heteroscedasticity under a linear assumption.
## assumption: data are already well-normalized (rlog)

check_heteroscedasticity_Pearson <- function(expr_mat, condition, plot = TRUE) {
  # expr_mat: numeric matrix, genes in rows, samples in columns (e.g. rlog + ComBat)
  # condition: factor or vector with same length as ncol(expr_mat)
  # plot: whether to draw diagnostic plots
  
  stopifnot(ncol(expr_mat) == length(condition))
  condition <- as.factor(condition)
  
  cat("=== GLOBAL MEAN–VARIANCE RELATIONSHIP ===\n")
  gene_means <- rowMeans(expr_mat)
  gene_vars  <- apply(expr_mat, 1, var)
  
  cor_global <- cor(gene_means, gene_vars, method = "pearson")
  lm_global  <- lm(gene_vars ~ gene_means)
  slope_global <- coef(lm_global)[2]
  
  cat(sprintf("Pearson r = %.3f | slope = %.4f\n", cor_global, slope_global))
  if (plot) {
    plot(gene_means, gene_vars, pch = 16, cex = 0.4,
         xlab = "Mean expression (rlog scale)",
         ylab = "Variance (rlog scale)",
         main = sprintf("Global mean–variance (r=%.2f)", cor_global))
    abline(lm_global, col = "red", lwd = 2)
  }
  
  cat("\n=== PER-CONDITION RELATIONSHIPS ===\n")
  results <- data.frame(
    condition = levels(condition),
    Spearman_r = NA,
    slope = NA,
    stringsAsFactors = FALSE
  )
  
  for (cond in levels(condition)) {
    idx <- which(condition == cond)
    if (length(idx) < 2) next  # skip groups with <2 samples
    
    gene_means <- rowMeans(expr_mat[, idx, drop = FALSE])
    gene_vars  <- apply(expr_mat[, idx, drop = FALSE], 1, var)
    
    r <- cor(gene_means, gene_vars, method = "pearson")
    s <- coef(lm(gene_vars ~ gene_means))[2]
    results[results$condition == cond, c("Pearson_r", "slope")] <- c(r, s)
    
    cat(sprintf("%-15s | r = %.3f | slope = %.4f\n", cond, r, s))
    
    if (plot) {
      plot(gene_means, gene_vars, pch = 16, cex = 0.4,
           xlab = sprintf("Mean (%s)", cond),
           ylab = "Variance",
           main = sprintf("%s (r=%.2f)", cond, r))
      abline(lm(gene_vars ~ gene_means), col = "red")
    }
  }
  
  invisible(list(global = list(r = cor_global, slope = slope_global),
                 per_condition = results))
}

check_heteroscedasticity_Pearson(CBrlogs, condition)
# === GLOBAL MEAN–VARIANCE RELATIONSHIP ===
#  Pearson r = -0.011 | slope = -0.0014

# === PER-CONDITION RELATIONSHIPS ===
# Rpv12.0         | r = -0.058  | slope = -0.0011
# Rpv12.1.0       | r = -0.049  | slope = -0.0026
# Rpv12.1.24      | r = -0.000  | slope = -0.0000
# Rpv12.1.3.0     | r = -0.086  | slope = -0.0024
# Rpv12.1.3.24    | r = -0.028  | slope = -0.0016
# Rpv12.1.3.6     | r = -0.038  | slope = -0.0027
# Rpv12.1.6       | r = -0.023  | slope = -0.0014
# Rpv12.24        | r = -0.175* | slope = -0.0018
# Rpv12.6         | r = 0.001   | slope = 0.0002
# Susceptible.0   | r = -0.013  | slope = -0.0014
# Susceptible.24  | r = -0.002  | slope = -0.0003
# Susceptible.6   | r = -0.243* | slope = -0.0022

check_heteroscedasticity_Pearson(CBrlogs, myTime)
# === PER-CONDITION RELATIONSHIPS ===
# 0               | r = 0.010  | slope = 0.0013
# 6               | r = -0.030 | slope = -0.0035
# 24              | r = -0.023 | slope = -0.0034

check_heteroscedasticity_Pearson(CBrlogs, myResistance)
# === PER-CONDITION RELATIONSHIPS ===
# Rpv1            | r = 0.017 | slope = 0.0018
# Rpv1.12         | r = 0.013 | slope = 0.0012
# Rpv1.12.3       | r = -0.047 | slope = -0.0023
# Susceptible     | r = -0.012 | slope = -0.0012


############ SIZE FACTOR + COMBAT ##############

############ SPEARMAN ######
check_heteroscedasticity(SFcombat, condition)
# === GLOBAL MEAN–VARIANCE RELATIONSHIP ===
#  Spearman r = -0.681** | slope = -0.1500

# === PER-CONDITION RELATIONSHIPS ===
# Rpv12.0         | r = -0.310* | slope = -0.1155
# Rpv12.1.0       | r = -0.451* | slope = -0.0806
# Rpv12.1.24      | r = -0.475* | slope = -0.1168
# Rpv12.1.3.0     | r = -0.295* | slope = -0.0839
# Rpv12.1.3.24    | r = -0.559** | slope = -0.0784
# Rpv12.1.3.6     | r = -0.479* | slope = -0.1027
# Rpv12.1.6       | r = -0.520** | slope = -0.0881
# Rpv12.24        | r = -0.540** | slope = -0.0716
# Rpv12.6         | r = -0.419* | slope = -0.1125
# Susceptible.0   | r = -0.521** | slope = -0.1244
# Susceptible.24  | r = -0.451* | slope = -0.1319
# Susceptible.6   | r = -0.611** | slope = -0.0726

check_heteroscedasticity(SFcombat, myTime)
# === PER-CONDITION RELATIONSHIPS ===
# 0               | r = -0.643 | slope = -0.1473
# 6               | r = -0.656 | slope = -0.1446
# 24              | r = -0.666 | slope = -0.1621

check_heteroscedasticity(SFcombat, myResistance)
# === PER-CONDITION RELATIONSHIPS ===
# Rpv1            | r = -0.600 | slope = -0.1160
# Rpv1.12         | r = -0.598 | slope = -0.1016
# Rpv1.12.3       | r = -0.662 | slope = -0.0967
# Susceptible     | r = -0.641 | slope = -0.1196

############ PEARSON #######
check_heteroscedasticity_Pearson(SFcombat, condition)
# === GLOBAL MEAN–VARIANCE RELATIONSHIP ===
#  Pearson r = -0.365 | slope = -0.1500

# === PER-CONDITION RELATIONSHIPS ===
# Rpv12.0         | r = -0.271 | slope = -0.1155
# Rpv12.1.0       | r = -0.448 | slope = -0.0806
# Rpv12.1.24      | r = -0.277 | slope = -0.1168
# Rpv12.1.3.0     | r = -0.399 | slope = -0.0839
# Rpv12.1.3.24    | r = -0.339 | slope = -0.0784
# Rpv12.1.3.6     | r = -0.353 | slope = -0.1027
# Rpv12.1.6       | r = -0.332 | slope = -0.0881
# Rpv12.24        | r = -0.500 | slope = -0.0716
# Rpv12.6         | r = -0.287 | slope = -0.1125
# Susceptible.0   | r = -0.325 | slope = -0.1244
# Susceptible.24  | r = -0.261 | slope = -0.1319
# Susceptible.6   | r = -0.520 | slope = -0.0726

check_heteroscedasticity_Pearson(SFcombat, myTime)
# === PER-CONDITION RELATIONSHIPS ===
# 0               | r = -0.329 | slope = -0.1473
# 6               | r = -0.339 | slope = -0.1446
# 24              | r = -0.363 | slope = -0.1621

check_heteroscedasticity_Pearson(SFcombat, myResistance)
# === PER-CONDITION RELATIONSHIPS ===
# Rpv1            | r = -0.367 | slope = -0.1160
# Rpv1.12         | r = -0.330 | slope = -0.1016
# Rpv1.12.3       | r = -0.471 | slope = -0.0967
# Susceptible     | r = -0.384 | slope = -0.1196

