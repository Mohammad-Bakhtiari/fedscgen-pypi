args <- commandArgs(trailingOnly = TRUE)

# Ensure we have at least one argument (the directory)
if (length(args) < 1) {
  stop("Please provide the directory path containing the .h5ad files.")
}

# Directory path
directory_path <- args[1]

# Check if directory exists
if(!dir.exists(directory_path)) {
  stop("Provided directory does not exist.")
}

files <- list.files(directory_path, pattern = "\\.h5ad$", full.names = TRUE)

# Check if files exist in the directory
for(file in files) {
  if(!file.exists(file)) {
    stop(paste("File", file, "does not exist in the provided directory."))
  }
}

# Check if files exist in the directory
for(file in files) {
  if(!file.exists(file)) {
    stop(paste("File", file, "does not exist in the provided directory."))
  }
}

library(kBET)
library(anndata)
library(Matrix)

# Define a function to get the summary for a given file and k0
get_kBET_summary <- function(adata, k_value) {
  result <- kBET(as.matrix(adata$obsm$pca_20), adata$obs$batch, k0=k_value, plot=FALSE, do.pca=FALSE)
  return(result$stats)
}

get_kBET_summary_on_pca <- function(adata, k_value) {
  result <- kBET(as.matrix(adata$X), adata$obs$batch, k0=k_value, plot=FALSE, do.pca=TRUE, dim.pca=20)
  return(result$stats)
}

get_approach_from_filename <- function(file_path) {
  file_name <- basename(file_path)
  file_name_no_ext <- tools::file_path_sans_ext(file_name)
  parts <- strsplit(file_name_no_ext, "_")[[1]]
  return(tail(parts, n = 1))
}

results_df <- data.frame()

for(i in seq_along(files)) {
  file <- files[i]
  print(file)
  adata <- read_h5ad(file)
  ds_size <- nrow(adata$X)
  k <- c(0.05, 0.10, 0.15, 0.20, 0.25)
  # Check the type and value of ds_size
  print(paste("ds_size type:", class(ds_size)))
  print(paste("ds_size value:", ds_size))

  # Compute k_values ensuring no overflow
  k_values <- sapply(k, function(perc) {
    val <- ds_size * perc
    if (is.infinite(val) || is.na(val)) {
      stop(paste("Overflow encountered for percentage:", perc))
    }
    as.integer(val)
  })

  for(index in 1:length(k)) {
    perc <- k[index]
    # k0 <- as.integer64(k_values[index])
    k0 <- k_values[index]

    print(k0)
    print("without pca")
    result_summary <- get_kBET_summary(adata, k0)

    temp_df <- data.frame(
      filename = basename(file),
      sample_size = k0,
      k_value = perc,
      ds_size = ds_size,
      approach = get_approach_from_filename(file)
    )
    temp_df <- cbind(temp_df, result_summary)

    results_df <- rbind(results_df, temp_df)
  }

}

output_file <- file.path(directory_path, "kBET_summary_results.csv")
write.csv(results_df, output_file)
