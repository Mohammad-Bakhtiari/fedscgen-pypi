import __init__
import copy
import argparse
import anndata
from functools import partial
import ast
import numpy as np
from fedscgen.FedScGen import FedScGen, ScGen
from fedscgen.utils import testset_combination, aggregate, aggregate_batch_sizes, remove_cell_types, combine_cell_types, \
    get_cuda_device
from fedscgen.plots import translate, plot_all_umaps, single_plot


def update_clients(clients, g_weights):
    """
    Update the clients with the global weights and return the updated weights and the number of samples in each client.
    Parameters
    ----------
    clients: list
    g_weights: weights of the global model

    Returns
    -------

    """
    weights = []
    n_samples = []
    for c in clients:
        # train on combined source and train part of target (both integrated)
        w, n_s = c.local_update(g_weights)
        weights.append(w)
        n_samples.append(n_s)
    return weights, n_samples


def main(args):
    adata = anndata.read_h5ad(args.adata)
    if args.combine:
        adata = combine_cell_types(adata, args.remove_cell_types, args.cell_key)
    else:
        adata = remove_cell_types(adata, args.remove_cell_types, args.cell_key)
    kwargs = {"ref_model_path": args.ref_model,
              "adata": adata,
              "hidden_layer_sizes": args.hidden_size,
              "z_dimension": args.z_dim,
              "batch_key": args.batch_key,
              "cell_key": args.cell_key,
              "lr": args.lr,
              "epoch": args.epoch,
              "batch_size": args.batch_size,
              "stopping": args.early_stopping_kwargs,
              "overwrite": False,
              "device": args.device
              }

    batch_type = type(adata.obs.batch.values.tolist()[0])
    args.batches = [batch_type(b) for b in args.batches]

    print(testset_combination(args.batches, args.batch_out))
    for test_batches in testset_combination(args.batches, args.batch_out):
        test_adata = adata[adata.obs[args.batch_key].isin(test_batches)].copy()
        print(f"{test_batches} as the test batches")
        train_batches = copy.deepcopy(args.batches)
        for batch in test_batches:
            train_batches.remove(batch)
        global_model = ScGen(init_model_path=args.init_model_path, **kwargs)
        global_weights = global_model.model.state_dict()
        global_model.is_trained_ = True
        global_model.model.eval()
        clients = []
        for client in range(1, args.n_clients + 1):
            print(f"Initializing client {client}...")
            kwargs["adata"] = adata[adata.obs[args.batch_key].isin([train_batches.pop()])].copy()
            c = FedScGen(**kwargs)
            clients.append(c)
        # training
        for r in range(1, args.n_rounds + 1):
            print(f"Round {r}/{args.n_rounds} of communication...")
            local_weights, local_n_samples = update_clients(clients, global_weights)
            print("Aggregating weights...")
            global_weights = aggregate(local_weights, local_n_samples)

            if args.per_round_snapshots:
                global_model.model.load_state_dict(global_weights)
                corrected_adata = global_model.model.batch_removal(adata,
                                                                   batch_key=args.batch_key,
                                                                   cell_label_key=args.cell_key,
                                                                   return_latent=True)
                corrected_adata.write(f"{args.output}/{translate(str(test_batches))}/corrected_{r}.h5ad")
                single_plot(corrected_adata, args.batch_key, args.cell_key,
                            f"{args.output}/{translate(str(test_batches))}",
                            f"corrected_{r}.png")
        if not args.per_round_snapshots:
            global_model.model.load_state_dict(global_weights)
            corrected_adata = evaluate(adata, clients, global_model, test_adata, test_batches, args.batches,
                                       args.batch_key,
                                       args.cell_key,
                                       args.output)
            corrected_adata.write(f"{args.output}/{translate(str(test_batches))}/corrected.h5ad")
        global_model.save(f"{args.output}/{translate(str(test_batches))}/trained_model", overwrite=True)


def evaluate(adata, clients, global_model, test_adata, test_batches, batches, batch_key, cell_key, output):
    """
    Evaluate the performance of the global model bu plotting UMAPs of source, target, and all data after correction.
    Parameters
    ----------
    adata
    clients
    global_model
    test_adata
    test_batches
    batches
    batch_key
    cell_key
    output

    Returns
    -------

    """
    sca_batch_removal = partial(global_model.model.batch_removal, batch_key=batch_key, cell_label_key=cell_key,
                                return_latent=True)
    plot_func = partial(plot_all_umaps, batch_key=batch_key, cell_key=cell_key)
    train_batches = copy.deepcopy(batches)
    for batch in test_batches:
        train_batches.remove(batch)
    trns_test_batches = translate(str(test_batches))
    if len(trns_test_batches) == 0:
        trns_test_batches = "BO0"
    train_adata = adata[adata.obs[batch_key].isin(train_batches)].copy()
    # For plotting the dataset set and trainset we are not bounded by privacy issues since it will not happen in real
    # world
    plot_func(uncorrected=train_adata,
              corrected=sca_batch_removal(train_adata),
              umap_directory=f"{output}/{trns_test_batches}/trainset")

    corrected_adata = sca_batch_removal(adata)
    plot_func(uncorrected=adata,
              corrected=corrected_adata,
              umap_directory=f"{output}/{trns_test_batches}/dataset")

    if len(test_adata) > 0:
        mean_latent_features = federated_evaluation(plot_func, clients, global_model, test_adata, test_batches)
        corrected = global_model.remove_batch_effect(adata, mean_latent_features)
        corrected.write(f"{output}/{trns_test_batches}/fed_corrected.h5ad")
        # find absolute difference between the corrected and corrected_adata
        abs_diff = np.abs(corrected.X - corrected_adata.X)
        print(f"Mean absolute difference between the corrected and corrected_adata: {np.mean(abs_diff)}")
        print(f"Standard deviation of the absolute difference between the corrected and corrected_adata: {np.std(abs_diff)}")
        print(f"Maximum absolute difference between the corrected and corrected_adata: {np.max(abs_diff)}")
        print(f"Minimum absolute difference between the corrected and corrected_adata: {np.min(abs_diff)}")
        print(f"Sum of the absolute difference between the corrected and corrected_adata: {np.sum(abs_diff)}")
    return corrected_adata


def federated_evaluation(plot_func, clients, global_model, test_adata, test_batches):
    """
    Since plotting a private test data can happen in real world and relies on global cell average, we implement
    it in federated fashion We assume only one client has a test data, however, with the same fashion we can
    extend it to more clients.

    Parameters
    ----------
    plot_func
    clients
    global_model
    test_adata
    test_batches

    Returns
    -------

    """

    batch_sizes = {}
    for client_id, client in enumerate(clients):
        batch_sizes[client_id] = client.find_batch_size()
    global_cell_sizes = aggregate_batch_sizes(batch_sizes)
    global_cell_avg = {}
    for i, c in global_cell_sizes.items():
        c_avg = clients[i].avg_local_cells(c)
        global_cell_avg.update(c_avg)
    plot_func(uncorrected=test_adata,
              corrected=global_model.remove_batch_effect(test_adata, global_cell_avg),
              umap_directory=f"{args.output}/{translate(str(test_batches))}/testset")
    return global_cell_avg

if __name__ == '__main__':
    """
    Federated scenarios:    
    4 client: put zero batch out
    3 clients: put 1 batches out
    2 clients: put 2 batches out
    """
    defalt_root = "/home/bba1658"
    parser = argparse.ArgumentParser()
    parser.add_argument("--init_model_path", type=str,
                        default=f"{defalt_root}/FedScGen/models/centralized/HumanPancreas")
    parser.add_argument("--adata", type=str, default=f"{defalt_root}/FedScGen/data/datasets/HumanPancreas.h5ad")
    parser.add_argument("--output", type=str,
                        default=f"{defalt_root}/FedScGen/results/scgen/federated/HumanPancreas/all/BO0-C5")
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
    args = parser.parse_args()

    args.hidden_size = [int(num) for num in args.hidden_size.replace(" ", "").split(",")]
    args.batches = [b.strip() for b in args.batches.split(",")]
    args.ref_model = f"{args.output}/{args.ref_model}"
    args.remove_cell_types = [c.strip() for c in args.remove_cell_types.strip().split(",")]
    args.device = get_cuda_device(args.gpu)
    main(args)
