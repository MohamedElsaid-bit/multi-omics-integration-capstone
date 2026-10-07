#!/usr/bin/env Rscript
# Baseline: what does each omics block achieve on its own, before any
# integration? PCA for visualization, then a PLS-DA classifier per block
# with cross-validated performance, so "integration helps" (script 03) has
# an honest number to beat rather than just a visual impression.
#
# Run from the repo root, after 01_explore_data.R:
#   Rscript scripts/02_single_omics_baseline.R

suppressMessages({
  library(mixOmics)
  library(ggplot2)
})

breast.TCGA <- readRDS("data/processed/breast_TCGA.rds")
blocks <- c("mrna", "mirna", "protein")
block_labels <- c(mrna = "mRNA expression", mirna = "miRNA expression", protein = "Protein (RPPA)")

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

set.seed(42)

pca_coords_all <- list()
baseline_rows <- list()

for (b in blocks) {
  X_train <- breast.TCGA$data.train[[b]]
  X_test <- breast.TCGA$data.test[[b]]
  y_train <- breast.TCGA$data.train$subtype
  y_test <- breast.TCGA$data.test$subtype

  # --- PCA, for visualization only, train + test combined ------------
  # (Unsupervised and descriptive, not part of the predictive pipeline
  # below, so combining splits here is not a leakage concern the way it
  # would be for the classifiers.)
  if (is.null(X_test)) {
    X_all <- X_train
    y_all <- y_train
    split_all <- rep("train", nrow(X_train))
  } else {
    X_all <- rbind(X_train, X_test)
    y_all <- factor(c(as.character(y_train), as.character(y_test)), levels = levels(y_train))
    split_all <- c(rep("train", nrow(X_train)), rep("test", nrow(X_test)))
  }

  pca_res <- pca(X_all, ncomp = 2, center = TRUE, scale = TRUE)
  pca_df <- data.frame(
    sample_id = rownames(X_all), block = block_labels[[b]],
    PC1 = pca_res$variates$X[, 1], PC2 = pca_res$variates$X[, 2],
    subtype = y_all, split = split_all
  )
  pca_coords_all[[b]] <- pca_df

  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = subtype)) +
    geom_point(size = 2.5, alpha = 0.8) +
    labs(
      title = sprintf("%s alone: PCA by subtype", block_labels[[b]]),
      x = sprintf("PC1 (%.1f%%)", 100 * pca_res$prop_expl_var$X[1]),
      y = sprintf("PC2 (%.1f%%)", 100 * pca_res$prop_expl_var$X[2])
    ) +
    theme_bw(base_size = 13)
  ggsave(sprintf("results/figures/pca_%s.png", b), p, width = 5.5, height = 4.5, dpi = 150)

  # --- PLS-DA baseline classifier, cross-validated -------------------
  plsda_res <- plsda(X_train, y_train, ncomp = 2)
  perf_res <- perf(plsda_res, validation = "Mfold", folds = 5, nrepeat = 10, progressBar = FALSE)
  cv_ber <- perf_res$error.rate$BER["comp2", "max.dist"]

  # The real test set has no protein block at all (see README); there is
  # no held-out test accuracy to report for that block alone.
  if (is.null(X_test)) {
    test_acc <- NA
    cat(sprintf(
      "%-18s CV balanced error rate = %.3f, test accuracy = NA (no protein block in the test set)\n",
      block_labels[[b]], cv_ber
    ))
  } else {
    test_pred <- predict(plsda_res, X_test, dist = "max.dist")
    test_class <- test_pred$class$max.dist[, 2]
    test_acc <- mean(test_class == y_test)
    cat(sprintf(
      "%-18s CV balanced error rate = %.3f, test accuracy = %.3f\n",
      block_labels[[b]], cv_ber, test_acc
    ))
  }

  baseline_rows[[b]] <- data.frame(
    model = block_labels[[b]], n_blocks = 1,
    cv_balanced_error_rate = cv_ber, test_accuracy = test_acc
  )
}

pca_coords_long <- do.call(rbind, pca_coords_all)
write.csv(pca_coords_long, "results/tables/pca_coordinates_by_block.csv", row.names = FALSE)

baseline_summary <- do.call(rbind, baseline_rows)
write.csv(baseline_summary, "results/tables/single_omics_baseline_performance.csv", row.names = FALSE)

cat("\nWrote PCA coordinates and single-omics baseline performance to results/tables/\n")
