args <- commandArgs(trailingOnly = TRUE)

# Ensure we have at least one argument (the directory)
if (length(args) < 3) {
  stop("Please provide the directory path containing the .h5ad files.")
}

set.seed(42)

# Directory path
raw <- args[1]
scgen <- args[2]
fedscgen <- args[3]
fedscgen_smpc <- args[4]
output_dir= args[5]


files <- c(raw = raw, scgen = scgen, fedscgen = fedscgen)

if (fedscgen_smpc != "none") {
  files <- c(files, fedscgen_smpc = fedscgen_smpc)
}


# Check if files exist in the directory
for(file in files) {
  if(!file.exists(file)) {
    stop(paste("File", file, "does not exist in the provided directory."))
  }
}

# library(zellkonverter)
library(anndata)
library(Matrix)
library(lisi)

results_df <- data.frame()

for(name in names(files)) {
  file_path <- files[name]
  adata <- read_h5ad(file_path)
  myPCA <- adata$obsm$pca_20
  lisi_meta_data <- adata$obs[, c('batch', 'cell_type')]

  lisi_label = c('batch', 'cell_type')

  lisi_res <- lisi::compute_lisi(myPCA, lisi_meta_data, c('batch', 'cell_type'))
  print(head(lisi_res))

  # Add an identifier column to lisi_res
  lisi_res$file_name <- name

  # Append the results to results_df
  results_df <- rbind(results_df, lisi_res)
}

# Write the results to a CSV file
write.csv(results_df, file.path(output_dir, "lisi_results.csv"), row.names = TRUE)
