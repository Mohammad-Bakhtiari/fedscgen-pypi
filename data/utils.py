import scanpy as sc
import pandas as pd
import os
import anndata
import sys
from pathlib import Path

# Add the parent directory to sys.path
parent_dir = str(Path(__file__).resolve().parent.parent)
if parent_dir not in sys.path:
    sys.path.append(parent_dir)
from fedscgen.plots import umap_plot, tsne_plot
from scipy.io import mmread


def normalize(adata, mtx=False):
    """
        Normalize according to https://nbviewer.org/github/Teichlab/bbknn_preprint/blob/main/pancreas.ipynb
    Parameters
    ----------
    adata

    Returns
    -------

    """
    if not mtx:
        adata.raw = adata.copy()

    # Step 1: Filter cells
    sc.pp.filter_cells(adata, min_genes=200)

    # Step 2: Get genes that are expressed in at least 3 cells across all batches
    sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
    filter_result = sc.pp.filter_genes_dispersion(adata.X, min_mean=0.0125, max_mean=2.5, min_disp=0.7)
    sc.pl.filter_genes_dispersion(filter_result)
    adata = adata[:, filter_result.gene_subset]
    sc.pp.log1p(adata)
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata)
    adata.obsm['X_pca'] *= -1
    return adata


def plot(adata, path, ds_name, batch_key="batch", cell_key="cell_type"):
    sc.settings.figdir = path
    sc.settings.set_figure_params(dpi=200, frameon=False)
    sc.set_figure_params(dpi=200)
    sc.set_figure_params(figsize=(4, 4))
    umap_plot(adata, wspace=0.6, batch_key=batch_key, cell_label_key=cell_key, plt_name=f"_{ds_name}.png")
    tsne_plot(adata, wspace=0.6, batch_key=batch_key, cell_label_key=cell_key, plt_name=f"_{ds_name}.png")


def prep(annotation: str, expression: str, path, name, adata_filename, batch_key, cell_key):
    batch_key = None if batch_key == "None" else batch_key
    if expression.__contains__(",") or annotation.__contains__(","):
        adata, annot_df = read_datasets(annotation, expression, path, unknown_batch=batch_key is None)
    else:
        adata, annot_df = read_dataset(annotation, expression, path)
    print(adata.X.shape)
    if not ".mtx" in expression:
        adata, logs = annotate(adata, annot_df)
        write_logs(logs, path)
    batch_key = "batch" if batch_key is None else batch_key
    adata.obs['org_batch'] = adata.obs[batch_key]
    adata.obs['batch'] = adata.obs[batch_key].astype("category").cat.codes.astype('str').astype('category')
    adata.obs['cell_type'] = adata.obs[cell_key]
    adata = normalize(adata)
    adata.write(os.path.join(path, adata_filename))
    # adata = anndata.read_h5ad(os.path.join(path, adata_filename))
    plot(adata, path, name)


def read_datasets(annotation, expression, path, unknown_batch=False):
    adata = []
    annot_df = []
    for batch_n, (anot, expr) in enumerate(zip(annotation.split(","), expression.split(",")), 1):
        adata_file, anot_file = read_dataset(anot, expr, path)
        if unknown_batch:
            anot_file['batch'] = batch_n
        adata.append(adata_file)
        annot_df.append(anot_file)
    adata = anndata.concat(adata, axis=1)
    annot_df = pd.concat(annot_df)
    return adata, annot_df


def write_logs(logs, path):
    f = open(os.path.join(path, "logs.txt"), "w")
    f.write(logs)


def annotate(adata, annot_df):
    adata = adata.transpose()
    logs = ""
    for column in annot_df.columns:
        unique = annot_df[column].value_counts()
        logs = log(logs, unique, column)
        try:
            adata.obs[column] = annot_df.loc[adata.obs_names, column]
        except Exception as e:
            print(f"Error processing column {column}: {e}")
            continue
    return adata, logs


def log(logs, unique, column):
    logs = f"{logs}{column} ==> number of unique values: {len(unique)}\n"
    if len(unique) < 50:
        logs = f"{logs}{unique}\n"
    return logs


def read_dataset(annotation, expression, path):
    annot_df = pd.read_csv(os.path.join(path, annotation), header=0, index_col=0, sep='\t')
    if ".mtx" in expression:
        matrix = mmread(os.path.join(path, expression)).T
        assert matrix.shape[0] == annot_df.shape[0], "Mismatch in number of cells between mtx file and df."
        adata = sc.AnnData(X=matrix, obs=annot_df, dtype='float32')
        adata.X = adata.X.tocsr()
    else:
        adata = sc.read_text(os.path.join(path, expression), delimiter='\t', first_column_names=True, dtype='float32')
    return adata, annot_df
