import argparse
import ast
import anndata
from fedscgen.utils import get_cuda_device, set_seed
from fedscgen.scenarios.centralized import run_centralized


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--model_path", type=str, default="./models/centralized/HumanPancreas/init_model")
    parser.add_argument("--data_path", type=str, default="./datasets/reproduce/pancreas-cleaned.h5ad")
    parser.add_argument("--output_path", type=str, default="./results/mapping/pancreas/reproduce")
    parser.add_argument("--epoch", type=int, default=1)
    parser.add_argument("--batch_key", type=str, default="batch")
    parser.add_argument("--cell_key", type=str, default="cell_type")
    parser.add_argument("--early_stopping_kwargs", type=ast.literal_eval,
                        default='{"early_stopping_metric": "val_loss", "patience": 20, "threshold": 0,'
                                ' "reduce_lr": True, "lr_patience": 13, "lr_factor": 0.1}')
    parser.add_argument("--just_umap", action='store_true', default=False)
    parser.add_argument("--hidden_layers_sizes", type=str, default="800,800")
    parser.add_argument("--z_dim", type=int, default=10)
    parser.add_argument("--gpu", type=int, default=0)
    parser.add_argument("--lr", type=float, default=0.001)
    parser.add_argument("--batch_size", type=int, default=50)
    parser.add_argument("--overwrite", action='store_true', default=False)
    parser.add_argument("--target_batches", type=str, default="")
    parser.add_argument("--ref_model", type=str, default="ref_model")
    parser.add_argument("--remove_cell_types", type=str, default="")
    parser.add_argument("--combine", action='store_true', default=False)
    parser.add_argument("--seed", type=int, default=42)

    return parser.parse_args()

def main():
    args = parse_args()
    set_seed(args.seed)

    args.target_batches = [b.strip() for b in args.target_batches.strip().split(",")] if args.target_batches else []
    args.hidden_layers_sizes = [int(s.strip()) for s in args.hidden_layers_sizes.split(",")]
    args.remove_cell_types = [c.strip() for c in args.remove_cell_types.strip().split(",")]
    args.ref_model = f"{args.output_path}/{args.ref_model}"
    args.device = get_cuda_device(args.gpu)

    adata = anndata.read_h5ad(args.data_path)
    run_centralized(adata, args)


if __name__ == '__main__':
    main()
