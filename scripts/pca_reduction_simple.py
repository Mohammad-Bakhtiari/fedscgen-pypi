import __init__
import argparse
import os
from fedscgen.utils import calc_obsm_pca, set_seed
import anndata

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate and Plot NMI for a set of adata files.')

    # Add arguments
    parser.add_argument('--path', required=True, help='Path to the raw adata file.')
    parser.add_argument('--n_components', type=int, default=50, help="Number of components for PCA")
    parser.add_argument('--output_dir', required=True, help='path to save the files.')
    # Parse the arguments
    args = parser.parse_args()
    set_seed()
    adata_files = calc_obsm_pca({"a": args.path}, args.n_components)
    adata_files["a"].write(args.output_dir)
