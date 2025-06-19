from fedscgen.FedScGen import ScGen
from fedscgen.utils import remove_cell_types, combine_cell_types


def run_centralized(adata, args):
    """Run centralized training using ScGen."""
    if args.combine:
        adata = combine_cell_types(adata, args.remove_cell_types, args.cell_key)
    else:
        adata = remove_cell_types(adata, args.remove_cell_types, args.cell_key)

    network = ScGen(
        init_model_path=args.model_path,
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
        device=args.device,
    )

    if len(args.target_batches) > 0:
        network.split(target_batches=args.target_batches)

    network.train_centralized(args.output_path)
