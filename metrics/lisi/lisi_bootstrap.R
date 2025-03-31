library(reticulate)
use_condaenv("r_eval", required = TRUE)

#!/usr/bin/env Rscript
library(anndata)
library(lisi)
library(Matrix)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Please provide four file paths (.h5ad) for Raw, scGen, FedscGen, FedscGen-SMPC, and an output directory.")
}

set.seed(42)
raw <- args[1]
scgen <- args[2]
fedscgen <- args[3]
fedscgen_smpc <- args[4]
output_dir <- args[5]
files <- c(raw=raw, scgen=scgen, fedscgen=fedscgen, "fedscgen-smpc"=fedscgen_smpc)

for (file in files) {
  if (!file.exists(file)) {
    stop(paste("File", file, "does not exist."))
  }
}

results_df <- data.frame()
n_boot <- 10  # Adjust to 50 or 100 for more precision if desired
for (name in names(files)) {
  file_path <- files[name]
  adata <- read_h5ad(file_path)
  myPCA <- adata$obsm$pca_20
  lisi_meta_data <- adata$obs[, c('batch', 'cell_type')]
  n_samples <- nrow(myPCA)

  boot_lisi <- replicate(n_boot, {
    idx <- sample(n_samples, replace=TRUE)
    compute_lisi(myPCA[idx,], lisi_meta_data[idx,], c('batch', 'cell_type'))
  }, simplify=FALSE)

  lisi_res <- do.call(rbind, boot_lisi)
  lisi_res$file_name <- name
  results_df <- rbind(results_df, lisi_res)
}

write.csv(results_df, file.path(output_dir, "lisi_results.csv"), row.names=FALSE)