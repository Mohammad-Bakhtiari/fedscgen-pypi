import argparse
import os
from fedscgen.utils import calc_obsm_pca

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate and Plot NMI for a set of adata files.')

    # Add arguments
    parser.add_argument('--raw', required=True, help='Path to the raw adata file.')
    parser.add_argument('--centralized', required=True, help='Path to the centralized adata file.')
    parser.add_argument('--federated', required=True, help='Path to the federated adata file.')
    parser.add_argument('--n_components', type=int, default=50, help="Number of components for PCA")
    parser.add_argument("--common_space", action="store_true", default=False, help="PCA on common space")
    parser.add_argument('--output_dir', required=True, help='path to save the files.')

    # Parse the arguments
    args = parser.parse_args()

    adata_file_paths = {
        'Raw': args.raw,
        'ScGen': args.centralized,
        'FedScGen': args.federated
    }

    adata_files = calc_obsm_pca(adata_file_paths, args.n_components, args.common_space)
    for name, adata in adata_files.items():
        adata.write(os.path.join(args.output_dir, f"{name}.h5ad"))
