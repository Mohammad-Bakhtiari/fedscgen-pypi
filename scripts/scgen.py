import __init__
import argparse
import ast
import anndata
from fedscgen.FedScGen import ScGen
from fedscgen.utils import remove_cell_types, combine_cell_types, get_cuda_device, set_seed


def main(adata, args):
    if args.combine:
        adata = combine_cell_types(adata, args.remove_cell_types, args.cell_key)
    else:
        adata = remove_cell_types(adata, args.remove_cell_types, args.cell_key)

    network = ScGen(init_model_path=args.model_path,
                    ref_model_path=args.ref_model,
                    adata=adata,
                    hidden_layer_sizes=args.hidden_layers_sizes,
                    z_dimension=args.z_dim,
                    batch_key=args.batch_key,
                    cell_key=args.cell_key,
                    lr=args.lr,
                    epoch=args.epoch,
                    batch_size=args.batch_size,
                    stopping=args.early_stopping_kwargs,
                    overwrite=args.overwrite,
                    device=args.device
                    )
    if len(args.target_batches) > 0:
        network.split(target_batches=args.target_batches)
    network.train_centralized(args.output_path)


if __name__ == '__main__':
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

    args = parser.parse_args()
    set_seed(args.seed)
    args.target_batches = [batch.strip() for batch in args.target_batches.strip().split(",")] if args.target_batches else []
    args.ref_model = f"{args.output_path}/{args.ref_model}"
    args.hidden_layers_sizes = [int(s.strip()) for s in args.hidden_layers_sizes.split(",")]
    args.remove_cell_types = [c.strip() for c in args.remove_cell_types.strip().split(",")]
    args.device = get_cuda_device(args.gpu)
    adata = anndata.read_h5ad(args.data_path)
    main(adata, args)
