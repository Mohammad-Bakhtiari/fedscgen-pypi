library(Matrix)
filename <- "./dataset8/dropviz_and_nuclei_combined_filtered_UMI.RDS"
# Read the RDS file
data <- readRDS(filename)

# Directly check if the data is a sparse matrix without converting it
if (inherits(data, "dgCMatrix")) {
    # Save the sparse matrix in Matrix Market format
    writeMM(data, sub(".RDS", ".mtx", filename))
} else {
    # If data is not a sparse matrix, consider other saving methods or throw an error
    stop("Data is not a sparse matrix. Cannot use writeMM.")
}
