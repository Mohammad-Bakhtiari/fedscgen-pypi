# Set CRAN mirror explicitly
options(repos = c(CRAN = "https://cloud.r-project.org/"))

# Check for BiocManager and install if it's not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install anndata from CRAN
install.packages("anndata")

# Install Matrix from CRAN
install.packages("Matrix")

# Install devtools from CRAN
install.packages("devtools")

# Load devtools
library(devtools)

# Install lisi from GitHub using devtools
install_github("immunogenomics/lisi")

# Install kBET from GitHub using devtools
install_github("theislab/kBET")