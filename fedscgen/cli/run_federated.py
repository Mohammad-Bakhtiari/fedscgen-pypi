import argparse
import ast
import anndata
from fedscgen.utils import get_cuda_device, set_seed
from fedscgen.scenarios.federated import run_federated


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--init_model_path", type=str, default="./models/centralized/HumanPancreas")
    parser.add_argument("--adata", type=str, default="./datasets/HumanPancreas.h5ad")
    parser.add_argument("--output", type=str, default="./results/federated/HumanPancreas")
    parser.add_argument("--epoch", type=int, default=1)
    parser.add_argument("--cell_key", type=str, default="cell_type")
    parser.add_argument("--batch_key", type=str, default="batch")
    parser.add_argument("--batches", type=str, default="0,1,2,3,4")
    parser.add_argument("--lr", type=float, default=0.001)
    parser.add_argument("--batch_size", type=int, default=50)
    parser.add_argument("--hidden_size", type=str, default="800,800")
    parser.add_argument("--z_dim", type=int, default=10)
    parser.add_argument("--gpu", type=int, default=0)
    parser.add_argument("--early_stopping_kwargs", type=ast.literal_eval,
                        default='{"early_stopping_metric": "val_loss", "patience": 20, "threshold": 0,'
                                ' "reduce_lr": True, "lr_patience": 13, "lr_factor": 0.1}')
    parser.add_argument("--batch_out", type=int, default=0)
    parser.add_argument("--n_clients", type=int, default=5)
    parser.add_argument("--n_rounds", type=int, default=1)
    parser.add_argument("--ref_model", type=str, default="ref_model")
    parser.add_argument("--remove_cell_types", type=str, default="")
    parser.add_argument("--combine", action='store_true', default=False)
    parser.add_argument("--per_round_snapshots", action='store_true', default=False)
    parser.add_argument("--smpc", action='store_true', default=False)
    parser.add_argument("--debug", action='store_true', default=False)
    parser.add_argument("--aggregation", type=str, default="fedavg", choices=["fedavg", "weighted_fedavg"])
    parser.add_argument("--seed", type=int, default=42)
    return parser.parse_args()


def main():
    args = parse_args()
    set_seed(args.seed)
    args.hidden_size = [int(x.strip()) for x in args.hidden_size.split(",")]
    args.batches = [x.strip() for x in args.batches.split(",")]
    args.ref_model = f"{args.output}/{args.ref_model}"
    args.remove_cell_types = [x.strip() for x in args.remove_cell_types.strip().split(",")]
    args.device = get_cuda_device(args.gpu)

    adata = anndata.read_h5ad(args.adata)
    run_federated(adata, args)


if __name__ == '__main__':
    main()
