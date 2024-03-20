
# Check for BiocManager and install if it's not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install zellkonverter using BiocManager
install.packages("anndata")

# Install Matrix from CRAN
install.packages("Matrix")

# Install devtools from CRAN
install.packages("devtools")

# Use devtools to install lisi from GitHub
devtools::install_github("immunogenomics/lisi")

library(devtools)
install_github('theislab/kBET')