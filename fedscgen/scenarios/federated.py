import copy
import argparse
import anndata
from functools import partial
import ast
import os
from fedscgen import FedScGen
from fedscgen.utils import testset_combination, aggregate, aggregate_batch_sizes, remove_cell_types, combine_cell_types, \
    get_cuda_device, abs_diff_centrally_corrected, check_adata_nan, check_weights_nan, set_seed
from fedscgen.plots import translate, single_plot
import torch
import numpy as np


def update_clients(clients, g_weights, smpc):
    """
    Update the clients with the global weights and return the updated weights and the number of samples in each client.
    Parameters
    ----------
    aggregation
    clients: list
    g_weights: weights of the global model

    Returns
    -------

    """
    if smpc:
        return [c.local_update(g_weights) for c in clients], None
    weights = []
    n_samples = []
    for c in clients:
        # train on combined source and train part of target (both integrated)
        w, n_s = c.local_update(g_weights)
        weights.append(w)
        n_samples.append(n_s)
    return weights, n_samples



def run_federated(adata, args):
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
              "device": args.device,
              "smpc": args.smpc,
              "aggregation": args.aggregation,
              "n_total_samples": len(adata.X),
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
        global_model = FedScGen(init_model_path=args.init_model_path, **kwargs)
        global_weights = global_model.model.state_dict()
        global_model.is_trained_ = True
        global_model.model.eval()
        clients = []
        for client in range(1, args.n_clients + 1):
            print(f"Initializing client {client}...")
            kwargs["adata"] = adata[adata.obs[args.batch_key].isin([train_batches.pop()])].copy()
            c = FedScGen(**kwargs)
            clients.append(c)
        test_clients = []
        for test_batch in test_batches:
            print(f"Initializing test client for batch {test_batch}...")
            kwargs["adata"] = adata[adata.obs[args.batch_key] == test_batch].copy()
            test_clients.append(FedScGen(**kwargs))
        # training
        for r in range(1, args.n_rounds + 1):
            print(f"Round {r}/{args.n_rounds} of communication...")

            state_dicts, n_samples = update_clients(clients, global_weights, args.smpc)
            print("Aggregating weights...")
            global_weights = aggregate(state_dicts, n_samples, args.smpc, list(global_weights.keys()))
            check_weights_nan(global_weights, "after aggregation", args.debug)

            if args.per_round_snapshots:
                correction_snapshot(clients, global_weights, f"{args.output}/{translate(str(test_batches))}",
                                    filename=f"corrected_{r}")
        for client in clients:
            client.model.load_state_dict(global_weights)
            client.model.eval()
        global_model.model.load_state_dict(global_weights)
        output_dir = f"{args.output}/{translate(str(test_batches))}"
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        evaluate_correction(adata, clients, test_clients, global_model, test_batches, args.batch_key, args.cell_key,
                            output_dir)
        global_model.save(f"{args.output}/{translate(str(test_batches))}/trained_model", overwrite=True)


def correction_snapshot(clients, global_weights, path, filename):
    original_state_dicts = [client.model.state_dict() for client in clients]
    for client in clients:
        client.model.load_state_dict(global_weights)
        client.model.eval()
    mlg = post_training_correction_wf(clients)
    fed_corrected = correct_local_data(clients, mlg)
    for client, original_sd in zip(clients, original_state_dicts):
        client.model.load_state_dict(original_sd)
    fed_corrected.write(f"{path}/{filename}.h5ad")
    single_plot(fed_corrected, args.batch_key, args.cell_key, path, f"{filename}.png")


def evaluate_correction(adata, clients, test_clients, global_model, test_batches, batch_key, cell_key, output):
    print(" Centralized correction of data")
    centrally_corrected = global_model.model.batch_removal(adata, batch_key=batch_key, cell_label_key=cell_key,
                                                           return_latent=True)
    check_adata_nan(adata)
    single_plot(centrally_corrected, batch_key, cell_key, output, "_centrally_corrected.png")
    print("Evaluating the performance of the batch effect correction Without new studies...")
    mlg = post_training_correction_wf(clients)
    fed_corrected = correct_local_data(clients, mlg)
    fed_corrected.write(f"{output}/fed_corrected.h5ad")
    single_plot(fed_corrected, batch_key, cell_key, output, "_fed_corrected.png")
    fed_corrected_with_new_studies = None
    if len(test_clients) > 0:
        print("Evaluating the performance of the batch effect correction With new studies...")
        for client in test_clients:
            client.model.load_state_dict(global_model.model.state_dict())
            client.is_trained_ = True
            client.model.eval()
            clients.append(client)
        mlg = post_training_correction_wf(clients)
        fed_corrected_with_new_studies = correct_local_data(clients, mlg)
        fed_corrected_with_new_studies.write(f"{output}/fed_corrected_with_new_studies.h5ad")
        single_plot(fed_corrected_with_new_studies, batch_key, cell_key, output, "_fed_corrected_with_new_studies.png")
    abs_diff = abs_diff_centrally_corrected(centrally_corrected, fed_corrected, fed_corrected_with_new_studies)
    abs_diff.to_csv(f"{output}/abs_diff.csv", sep=",", index=True)


def correct_local_data(clients, mlg):
    fed_corrected = [client.remove_batch_effect(mlg) for client in clients]
    fed_corrected = anndata.AnnData.concatenate(*fed_corrected, batch_key="concat_batch", index_unique=None)
    if "concat_batch" in fed_corrected.obs.columns:
        del fed_corrected.obs["concat_batch"]
    return fed_corrected


def post_training_correction_wf(clients):
    batch_sizes = {}
    print("First round: finding dominant batches...")
    for client_id, client in enumerate(clients):
        batch_sizes[client_id] = client.find_batch_size()
    global_cell_sizes = aggregate_batch_sizes(batch_sizes)
    print("Second round: finding mean latent genes of dominant batches per cell type...")
    mean_latent_genes = {}
    for i, c in global_cell_sizes.items():
        c_avg = clients[i].avg_local_cells(c)
        mean_latent_genes.update(c_avg)
    return mean_latent_genes