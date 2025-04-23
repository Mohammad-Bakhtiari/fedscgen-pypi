# Set CRAN mirror explicitly (reliable and fast)
options(repos = c(CRAN = "https://cloud.r-project.org/"))

# Install CRAN packages if not already installed
cran_packages <- c("anndata", "Matrix", "devtools", "reticulate")
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# Install BiocManager if missing
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Load devtools to install from GitHub
library(devtools)

# Install GitHub packages if not already installed
if (!requireNamespace("lisi", quietly = TRUE)) {
  install_github("immunogenomics/lisi")
}
if (!requireNamespace("kBET", quietly = TRUE)) {
  install_github("theislab/kBET")
}
