#!/usr/bin/env Rscript
# DIABLO (Data Integration Analysis for Biomarker discovery using Latent
# cOmponents, Singh et al. 2019): a supervised model that uses all three
# omics blocks together to predict subtype, rather than any one block
# alone. This is the actual integration step; script 02 is its baseline.
#
# Run from the repo root, after 02_single_omics_baseline.R:
#   Rscript scripts/03_diablo_integration.R

suppressMessages({
  library(mixOmics)
  library(ggplot2)
})

breast.TCGA <- readRDS("data/processed/breast_TCGA.rds")
blocks <- c("mrna", "mirna", "protein")

X_train <- list(
  mrna = breast.TCGA$data.train$mrna,
  mirna = breast.TCGA$data.train$mirna,
  protein = breast.TCGA$data.train$protein
)
Y_train <- breast.TCGA$data.train$subtype

# The real test set has no protein (RPPA) block at all (see README). DIABLO
# can predict from a partial set of blocks, which is the actual point of
# this test: can mRNA + miRNA alone, run through a model trained on all
# three blocks, still predict subtype on samples missing the third block.
X_test <- list(
  mrna = breast.TCGA$data.test$mrna,
  mirna = breast.TCGA$data.test$mirna
)
Y_test <- breast.TCGA$data.test$subtype

# Design: every pair of blocks is connected with a moderate correlation
# weight (0.1), the standard starting point in the mixOmics DIABLO
# vignette; it is not tuned here (see Limitations in the README).
design <- matrix(0.1, nrow = 3, ncol = 3, dimnames = list(blocks, blocks))
diag(design) <- 0

set.seed(42)
# Sparse DIABLO (block.splsda): keepX forces each component to use only a
# handful of features per block, which is the actual point of the method
# (a small candidate biomarker signature per omics layer, not "everything,
# ranked"). 15 per component per block is a reasonable fixed value, not
# tuned via tune.block.splsda's grid search, which is costly and listed
# under Future improvements instead.
keepX <- list(mrna = c(15, 15), mirna = c(15, 15), protein = c(15, 15))
diablo_model <- block.splsda(X_train, Y_train, ncomp = 2, keepX = keepX, design = design)

perf_res <- perf(diablo_model, validation = "Mfold", folds = 5, nrepeat = 10, progressBar = FALSE)
cv_ber <- perf_res$WeightedVote.error.rate$max.dist["Overall.BER", "comp2"]
cat(sprintf("DIABLO (3 blocks integrated): CV balanced error rate = %.3f\n", cv_ber))

test_pred <- predict(diablo_model, newdata = X_test, dist = "max.dist")
test_class <- test_pred$WeightedVote$max.dist[, "comp2"]
test_acc <- mean(test_class == Y_test)
cat(sprintf("DIABLO (3 blocks integrated): test accuracy = %.3f\n", test_acc))

integrated_row <- data.frame(
  model = "All 3 blocks integrated (DIABLO)", n_blocks = 3,
  cv_balanced_error_rate = cv_ber, test_accuracy = test_acc
)

baseline <- read.csv("results/tables/single_omics_baseline_performance.csv")
performance_summary <- rbind(baseline, integrated_row)
write.csv(performance_summary, "results/tables/model_performance_summary.csv", row.names = FALSE)
cat("\nFull single-block vs integrated comparison:\n")
print(performance_summary)

# Confusion matrix on the held-out test set for the integrated model.
confusion <- table(true = Y_test, predicted = test_class)
write.csv(as.data.frame(confusion), "results/tables/diablo_test_confusion_matrix.csv", row.names = FALSE)

# Sample plot: DIABLO's shared latent component 1, one panel per block,
# colored by true subtype.
png("results/figures/diablo_sample_plot.png", width = 1200, height = 500, res = 120)
plotIndiv(diablo_model, legend = TRUE, title = "DIABLO component 1, by block")
dev.off()

# Which features in each block actually drive the integrated model.
loadings_rows <- list()
for (b in blocks) {
  loadings <- selectVar(diablo_model, block = b, comp = 1)[[b]]$value
  loadings_rows[[b]] <- data.frame(
    block = b, feature = rownames(loadings), loading = loadings$value.var, component = 1
  )
}
top_features <- do.call(rbind, loadings_rows)
top_features <- top_features[order(top_features$block, -abs(top_features$loading)), ]
write.csv(top_features, "results/tables/diablo_top_features_comp1.csv", row.names = FALSE)

cat(sprintf(
  "\n%d features selected on component 1 across the 3 blocks (%s)\n",
  nrow(top_features), paste(table(top_features$block), collapse = "/")
))

saveRDS(diablo_model, "data/processed/diablo_model.rds")
cat("\nWrote performance summary, confusion matrix, sample plot, and top features to results/\n")
