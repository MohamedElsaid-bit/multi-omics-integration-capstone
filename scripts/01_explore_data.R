#!/usr/bin/env Rscript
# Load the real TCGA breast cancer multi-omics dataset bundled with mixOmics
# (Singh et al. 2019, the DIABLO paper) and report basic dimensions: three
# omics blocks (mRNA, miRNA, protein/RPPA) measured on the same patients,
# labeled by molecular subtype.
#
# Run from the repo root: Rscript scripts/01_explore_data.R

suppressMessages(library(mixOmics))

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

data(breast.TCGA)

blocks <- c("mrna", "mirna", "protein")
cat("Training set:\n")
for (b in blocks) {
  d <- dim(breast.TCGA$data.train[[b]])
  cat(sprintf("  %-8s %4d samples x %5d features\n", b, d[1], d[2]))
}
cat("\nTest set:\n")
for (b in blocks) {
  if (is.null(breast.TCGA$data.test[[b]])) {
    cat(sprintf("  %-8s not present in the test set (see README)\n", b))
    next
  }
  d <- dim(breast.TCGA$data.test[[b]])
  cat(sprintf("  %-8s %4d samples x %5d features\n", b, d[1], d[2]))
}

cat("\nTraining set subtype counts:\n")
print(table(breast.TCGA$data.train$subtype))
cat("\nTest set subtype counts:\n")
print(table(breast.TCGA$data.test$subtype))

# All three blocks must share the same sample IDs within a split; confirm
# that plainly rather than assuming it.
train_ids_match <- identical(rownames(breast.TCGA$data.train$mrna), rownames(breast.TCGA$data.train$mirna)) &&
  identical(rownames(breast.TCGA$data.train$mrna), rownames(breast.TCGA$data.train$protein))
cat(sprintf("\nTraining set: sample IDs identical across all 3 blocks: %s\n", train_ids_match))

# Sample level metadata export (also used as a Power BI dimension table).
metadata <- data.frame(
  sample_id = c(rownames(breast.TCGA$data.train$mrna), rownames(breast.TCGA$data.test$mrna)),
  subtype = c(as.character(breast.TCGA$data.train$subtype), as.character(breast.TCGA$data.test$subtype)),
  split = c(rep("train", nrow(breast.TCGA$data.train$mrna)), rep("test", nrow(breast.TCGA$data.test$mrna)))
)
write.csv(metadata, "results/tables/sample_metadata.csv", row.names = FALSE)

saveRDS(breast.TCGA, "data/processed/breast_TCGA.rds")
cat(sprintf("\nWrote sample metadata for %d samples to results/tables/sample_metadata.csv\n", nrow(metadata)))
