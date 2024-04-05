import anndata
from sklearn.metrics import f1_score, adjusted_rand_score
from collections import Counter
from sklearn.decomposition import PCA
import networkx as nx
from sklearn.neighbors import NearestNeighbors
import scanpy as sc
from scarches.metrics import asw
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from sklearn.cluster import KMeans
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests
from sklearn.metrics import silhouette_score
import pandas as pd

HEX_COLORS = [
    '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b',
    '#e377c2', '#7f7f7f', '#bcbd22', '#17becf', '#1a55FF', '#55a868',
    '#c44e52', '#8172b3', '#ccb974', '#64B5CD', '#FFB447', '#82C341',
    '#D17049', '#705696', '#FF69B4', '#BA55D3', '#CD5C5C', '#FFA07A',
    '#B0E0E6', '#FFD700', '#F08080', '#90EE90', '#20B2AA', '#778899',
    '#00FA9A', '#6B8E23', '#FF00FF', '#4682B4', '#008080', '#40E0D0',
    '#EE82EE', '#F4A460', '#DAA520', '#8B008B', '#800080', '#191970',
    '#000080', '#808000', '#FFFF00', '#00FF00', '#00FFFF', '#FF4500',
    '#FF6347', '#FFD700', '#FFA500', '#FF4500', '#DC143C', '#FF0000',
    '#B22222', '#8B0000', '#FFC0CB', '#FFB6C1', '#FF69B4', '#FF1493',
    '#C71585', '#DB7093', '#FFA07A', '#FF7F50', '#FF6347', '#FF4500',
    '#FF8C00', '#FFD700', '#FFFF00', '#ADFF2F', '#7FFF00', '#7CFC00',
    '#00FF00', '#32CD32', '#98FB98', '#00FA9A', '#00FF7F', '#3CB371',
    '#2E8B57', '#228B22', '#008000', '#006400', '#9ACD32', '#6B8E23',
    '#808000', '#556B2F', '#66CDAA', '#8FBC8F', '#20B2AA', '#008B8B',
    '#008080', '#00FFFF', '#00FFFF', '#E0FFFF', '#AFEEEE', '#7FFFD4',
    '#40E0D0', '#48D1CC', '#00CED1', '#5F9EA0', '#4682B4', '#6495ED',
    '#00BFFF', '#1E90FF', '#4169E1', '#0000FF', '#0000CD', '#00008B',
    '#000080', '#191970', '#FFF8DC', '#FFEBCD', '#FFE4C4', '#FFDEAD',
    '#F5DEB3', '#DEB887', '#D2B48C', '#BC8F8F', '#F4A460', '#DAA520',
    '#B8860B', '#CD853F', '#D2691E', '#8B4513', '#A0522D', '#A52A2A',
    '#800000', '#FFFFFF', '#FFFAFA', '#F0FFF0', '#F5FFFA', '#F0FFFF',
    '#F0F8FF', '#F8F8FF', '#F5F5F5', '#FFF5EE', '#F5F5DC', '#FDF5E6',
    '#FFFAF0', '#FFFFF0', '#FAEBD7', '#FAF0E6', '#FFF0F5', '#FFE4E1',
]

DATASETS = ["HumanDendriticCells",
            "MouseCellAtlas",
            "HumanPancreas",
            "PBMC",
            "CellLine",
            "MouseRetina",
            "MouseBrain",
            "MouseHematopoieticStemProgenitorCells"
            ]
DATASETS_ACRONYM = ["HDC", "MCA", "HP", "PBMC", "CL", "MR", "MB", "MHSPC"]
MODELS = ["Raw", "ScGen", "FedScGen"]
MODEL_COLORS = {k: v for k, v in zip(MODELS, sns.color_palette("viridis", 3))}
DATASETS_COLORS = {k: v for k, v in zip(DATASETS_ACRONYM, sns.color_palette("viridis", 8))}
DATASETS_MARKERS = {k: v for k, v in zip(DATASETS_ACRONYM, ['o', 's', '^', 'p', 'D', 'v', '*', 'H'])}
MINORITY_CLASSES = [
    "",
    "Epithelial,Dendritic,Smooth-muscle,NK",
    "stellate,endothelial,mesenchymal,macrophage,mast,epsilon,schwann,t_cell,MHC class II",
    "Plasmacytoid dendritic cell,Megakaryocyte,Hematopoietic stem cell",
    "",
    "ganglion,vascular_endothelium,horizontal,fibroblasts,microglia,pericytes,astrocytes",
    "Olfactory ensheathing cells,Choroid_plexus,Mitotic",
    "MPP,LTHSC,LMPP,Unsorted"
]


def plot_metrics_with_circles(df, plot_name):
    # Seaborn settings
    sns.set(style="whitegrid")

    # Get unique metrics and datasets
    metric_keys = df.columns[2:].tolist()
    dataset_keys = df['Dataset'].unique().tolist()
    fig, ax = plt.subplots(1, 1, figsize=(len(dataset_keys), len(metric_keys)))

    # Turn off grid lines
    ax.grid(False)
    # Define some colors, assuming up to three datasets
    colors = sns.color_palette("viridis", len(dataset_keys))

    y_labels = []

    for y, metric in enumerate(metric_keys):
        for x, dataset in enumerate(dataset_keys):
            value = df.loc[df['Dataset'] == dataset, metric].values[0]

            # Create a circle with radius proportional to the metric value
            circle = plt.Circle((x + 0.5, y + 0.5), np.sqrt(value) * 0.4,
                                color=colors[x], alpha=0.6, edgecolor='black')

            ax.add_artist(circle)

            # Add the value inside the circle
            plt.text(x + 0.5, y + 0.5, f"{value:.2f}", fontsize=12, ha='center', va='center')

        y_labels.append(metric)

    ax.set_yticks(np.arange(len(metric_keys)) + 0.5)
    ax.set_yticklabels(y_labels, fontsize=18)
    ax.set_xticks(np.arange(len(dataset_keys)) + 0.5)
    ax.set_xticklabels(dataset_keys, fontsize=18)
    ax.set_aspect('equal', adjustable='box')

    # Add horizontal grid lines
    for y in np.arange(len(metric_keys)):
        ax.axhline(y, linestyle='--', color='grey')

    # plt.title("Metrics Represented by Circles")
    plt.xlim([0, len(dataset_keys)])
    plt.ylim([0, len(metric_keys)])
    plt.tight_layout()
    plt.savefig(plot_name, dpi=300)
    plt.close()


def bar_plot(df, plot_name):
    # Plot
    fig, axs = plt.subplots(1, 9, figsize=(40, 6))
    sns.barplot(x='Dataset', y='NMI', data=df, ax=axs[0])
    axs[0].set_title('NMI')
    sns.barplot(x='Dataset', y='GC', data=df, ax=axs[1])
    axs[1].set_title('Graph Connectivity')
    sns.barplot(x='Dataset', y='ILF1', data=df, ax=axs[2])
    axs[2].set_title('Isolated Label F1')
    sns.barplot(x='Dataset', y='ARI', data=df, ax=axs[3])
    axs[3].set_title('ARI')
    sns.barplot(data=df, x='Dataset', y='EBM', ci=None, ax=axs[4])
    axs[4].set_title('EBM')
    sns.barplot(data=df, x='Dataset', y='ILS', ci=None, ax=axs[5])
    axs[5].set_title('ILS')
    sns.barplot(data=df, x='Dataset', y='KNN Acc', ci=None, ax=axs[6])
    axs[6].set_title('KNN accuracy')
    sns.barplot(data=df, x='Dataset', y='ASW_B', ci=None, ax=axs[7])
    axs[7].set_title('Batch ASW')
    sns.barplot(data=df, x='Dataset', y='ASW_B', ci=None, ax=axs[8])
    axs[8].set_title('Cell type ASW')

    plt.tight_layout()
    plt.savefig(plot_name)
    plt.close()


def ari_score(adata, label_key):
    sc.pp.neighbors(adata)
    sc.tl.louvain(adata)
    louvain_labels = adata.obs['louvain']
    # Assuming 'label_key' holds the true cell type labels
    true_labels = adata.obs[label_key]
    # Compute ARI between true labels and Louvain clustering labels
    return adjusted_rand_score(true_labels, louvain_labels)


def isolated_label_f1_score(adata, label_key, cluster_key):
    # Count instances for each label and identify isolated labels
    label_counts = Counter(adata.obs[label_key])
    min_count = min(label_counts.values())
    isolated_labels = [label for label, count in label_counts.items() if count == min_count]
    f1_scores = []
    for isolated_label in isolated_labels:
        # Find the cluster with the most instances of the isolated label
        subset = adata[adata.obs[label_key] == isolated_label]
        cluster_counts = Counter(subset.obs[cluster_key])
        main_cluster = max(cluster_counts, key=cluster_counts.get)

        # Calculate F1 score for isolated label within the main cluster
        cluster_subset = adata[adata.obs[cluster_key] == main_cluster]
        true_labels = (cluster_subset.obs[label_key] == isolated_label).astype(int)
        pred_labels = [1] * len(true_labels)
        f1 = f1_score(true_labels, pred_labels)
        f1_scores.append(f1)

    return sum(f1_scores) / len(f1_scores) if f1_scores else 0


def compute_latent_representations(raw_path, centralized_path, federated_path, n_components=50):
    """
    Compute the latent representations of the provided adata files using PCA.

    ... [rest of the docstring]
    """

    # Load adata files
    adata_files = {
        'Raw': anndata.read_h5ad(raw_path),
        'ScGen': anndata.read_h5ad(centralized_path),
        'FedScGen': anndata.read_h5ad(federated_path)
    }

    # Initialize PCA
    pca = PCA(n_components=n_components)

    # Fit PCA on raw data
    pca.fit(adata_files['Raw'].X)

    latent_representations = {}
    for key, adata in adata_files.items():
        transformed_data = pca.transform(adata.X)
        latent_adata = anndata.AnnData(transformed_data, obs=adata.obs)  # copying the observations
        latent_representations[key] = latent_adata

    return latent_representations, adata_files


def calc_obsm_pca(adata_file_paths, n_components=50, common_space=False):
    adata_files = {}
    if common_space:
        pca = PCA(n_components=n_components)
    for counter, (key, path) in enumerate(adata_file_paths.items()):
        adata = anndata.read_h5ad(path)
        if not common_space:
            pca = PCA(n_components=n_components)
            pca.fit(adata.X)
        elif counter == 0:
            pca.fit(adata.X)

        adata.obsm[f'pca_{n_components}'] = pca.transform(adata.X)
        adata_files[key] = adata

    return adata_files


def knn_graph(X, k=5):
    nbrs = NearestNeighbors(n_neighbors=k + 1, algorithm='ball_tree').fit(X)
    distances, indices = nbrs.kneighbors(X)
    G = nx.Graph()
    for i in range(len(indices)):
        for j in indices[i][1:]:  # Skip the first entry (self)
            G.add_edge(i, j)
    return G


def largest_connected_component(G):
    largest = max(nx.connected_components(G), key=len)
    return G.subgraph(largest)


def graph_connectivity_score(adata, k=5):
    total_score = 0
    cell_types = adata.obs['cell_type'].unique()
    for cell_type in cell_types:
        subset = adata[adata.obs['cell_type'] == cell_type, :]
        X = subset.X.toarray() if hasattr(subset.X, 'toarray') else subset.X  # Ensure dense representation
        k_adj = min(X.shape[0] - 1, k)  # Adjust k value
        G = knn_graph(X, k_adj)
        LCC = largest_connected_component(G)
        score = len(LCC.nodes) / len(subset)
        total_score += score
    return total_score / len(cell_types)


def compute_f1_ari(ari_batch, ari_cell_type):
    ari_batch_norm = 1 - ari_batch
    return 2 * ari_batch_norm * ari_cell_type / (ari_batch_norm + ari_cell_type)


def evaluate_ari(adata, adata_latent, label_key, batch_key):
    n_clusters = len(np.unique(adata.obs[label_key]))
    sample_idx = np.random.choice(adata.shape[0], int(adata.shape[0] * 0.8), replace=False)

    # Subsample the data
    adata_subsampled = adata[sample_idx, :]
    adata_latent_subsampled = adata_latent[sample_idx, :]

    # Perform k-means clustering on PCA components
    predicted_labels = KMeans(n_clusters=n_clusters, random_state=0).fit(adata_latent_subsampled.X).labels_

    # Compute ARI for cell type purity
    ari_cell_type = adjusted_rand_score(adata_latent_subsampled.obs[label_key], predicted_labels)

    # Only consider cells whose types are present in all batches
    common_cells = adata_subsampled.obs.groupby(batch_key).filter(
        lambda x: len(np.unique(x[label_key])) == len(adata_subsampled.obs[label_key].unique()))

    ari_batch = adjusted_rand_score(common_cells[batch_key], predicted_labels[common_cells.index])

    return ari_batch, ari_cell_type


def perform_wilcoxon_bh(data):
    # Pairwise Wilcoxon tests
    methods = data.columns
    p_values = []

    for i in range(len(methods)):
        for j in range(i + 1, len(methods)):
            _, p_value = wilcoxon(data[methods[i]], data[methods[j]])
            p_values.append(p_value)

    # Benjamini and Hochberg correction
    _, corrected_p_values, _, _ = multipletests(p_values, method='fdr_bh')
    return corrected_p_values


def run_ari_iterations(adata, adata_latent, label_key, batch_key, iterations=20):
    ari_batch_scores = []
    ari_cell_type_scores = []

    for _ in range(iterations):
        ari_batch, ari_cell_type = evaluate_ari(adata, adata_latent, label_key, batch_key)
        ari_batch_scores.append(ari_batch)
        ari_cell_type_scores.append(ari_cell_type)

    # Calculate the median ARI scores
    median_ari_batch = np.median(ari_batch_scores)
    median_ari_cell_type = np.median(ari_cell_type_scores)

    # Normalize the median ARI scores to range between 0 and 1
    f1 = compute_f1_ari(median_ari_batch, median_ari_cell_type)

    return f1, median_ari_batch, median_ari_cell_type


def compare_correction_methods(adata_list, adata_latent_list, label_key, batch_key, iterations=20):
    results = {"Method": [], "F1 Score": [], "ARI Batch": [], "ARI Cell Type": []}
    ari_scores = {}

    for adata, adata_latent in zip(adata_list, adata_latent_list):
        f1, ari_batch, ari_cell_type = run_ari_iterations(adata, adata_latent, label_key, batch_key, iterations)

        method_name = adata.uns['method_name']  # assuming a method name is stored for each adata

        results["Method"].append(method_name)
        results["F1 Score"].append(f1)
        results["ARI Batch"].append(ari_batch)
        results["ARI Cell Type"].append(ari_cell_type)

        ari_scores[method_name] = [run_ari_iterations(adata, adata_latent, label_key, batch_key, 1)[1] for _ in
                                   range(iterations)]

    results_df = pd.DataFrame(results)
    ari_scores_df = pd.DataFrame(ari_scores)

    corrected_p_values = perform_wilcoxon_bh(ari_scores_df)

    return results_df, corrected_p_values


def compute_ils(adata, cell_key, batch_key, metric='euclidean'):
    """
    Compute silhouette score over isolated labels in an AnnData object.

    adata: AnnData
        AnnData instance containing the single-cell RNA-seq data.

    cell_key: str
        Key in adata.obs to access the cell type labels.

    batch_key: str
        Key in adata.obs to access the batch labels.

    metric: string, (default='euclidean')
        The metric to use when calculating distance between instances in feature space.

    Returns: float
        Averaged silhouette score over isolated labels.
    """

    # Check if provided keys exist
    if cell_key not in adata.obs.columns or batch_key not in adata.obs.columns:
        raise ValueError("Provided cell_key or batch_key doesn't exist in adata.obs.")

    # Extract cell and batch labels
    cell_labels = adata.obs[cell_key].values
    batch_labels = adata.obs[batch_key].values

    # Identify isolated labels
    isolated_scores = []
    unique_cell_labels = np.unique(cell_labels)

    for label in unique_cell_labels:
        subset_indices = np.where(cell_labels == label)[0]

        # Check if this label is isolated (appears in one batch only)
        unique_batches_for_label = np.unique(batch_labels[subset_indices])
        if len(unique_batches_for_label) == 1:
            X_subset = adata.X[subset_indices, :]
            labels_subset = np.array([1 if l == label else 0 for l in cell_labels[subset_indices]])

            # We should ensure there's more than one label in the subset for silhouette_score
            if len(np.unique(labels_subset)) > 1:
                score = silhouette_score(X_subset, labels_subset, metric=metric)
                isolated_scores.append(score)

    return np.mean(isolated_scores) if isolated_scores else 0


def knn_accuracy(latent_adata, cell_key, k=15):
    # Extract latent space and labels
    X = latent_adata.X
    labels = latent_adata.obs[cell_key].values

    # Compute kNN
    nbrs = NearestNeighbors(n_neighbors=k + 1).fit(X)
    distances, indices = nbrs.kneighbors(X)

    # Compute accuracy for each cell
    accuracies = []
    for idx, neighbors in enumerate(indices):
        correct_neighbors = sum(1 for i in neighbors[1:] if labels[i] == labels[idx])  # exclude the cell itself
        accuracies.append(correct_neighbors / k)

    # Compute average accuracy for each cell type
    unique_labels = np.unique(labels)
    avg_accuracies = []
    for label in unique_labels:
        mask = (labels == label)
        avg_accuracies.append(np.mean([accuracies[i] for i, flag in enumerate(mask) if flag]))

    # Compute the final average
    final_avg_accuracy = np.mean(avg_accuracies)

    return final_avg_accuracy


def bar_plot_subplot(df, plot_name):
    """
    In 10 rows for 10 communication rounds and 3 columns for NMI, ARI, and EBM bar plot the values for 10 epochs
    Parameters
    find scgen in approach and plot it in horizontal line
    ----------
    df
    plot_name

    Returns
    -------

    """
    sns.set(style="whitegrid")
    fig, ax = plt.subplots(3, 10, figsize=(30, 9))
    max_value = [df['NMI'].max(), df['ARI'].max(), df['EBM'].max()]
    min_values = [df['NMI'].min(), df['ARI'].min(), df['EBM'].min()]
    scgen = df[df['Approach'] == "scGen"]
    df = df[df['Approach'] != "scGen"]
    df.drop(columns=["Approach"], inplace=True)
    for i, metric in enumerate(["NMI", "ARI", "EBM"]):
        for round in range(1, 11):
            target_ax = ax[i, round - 1]
            data = df[df['Round'] == round]
            sns.barplot(x="Epoch", y=metric, data=data, ax=target_ax, palette="viridis")
            target_ax.axhline(scgen[metric].values[0], color='r', linestyle='--')

            target_ax.set_ylabel("")
            target_ax.set_xlabel("")
            target_ax.set_ylim(np.floor(min_values[i] * 10) / 10, np.ceil(max_value[i] * 10) / 10)

            if round == 1:
                target_ax.set_ylabel(f"{metric}", fontsize=22)
                # show only three ticks
                target_ax.set_yticks([np.floor(min_values[i] * 10) / 10, np.ceil(max_value[i] * 10) / 10])
                target_ax.set_yticklabels([np.floor(min_values[i] * 10) / 10, np.ceil(max_value[i] * 10) / 10],
                                          fontsize=18)
                if i < 2:
                    target_ax.axes.get_xaxis().set_visible(False)
            if i == 2:
                target_ax.set_xticks([0, 4, 9])
                # target_ax.set_xlabel(f"{round}", fontsize=22)
                target_ax.set_xticklabels(['E1', 'E5', 'E10'], fontsize=16)
                target_ax.set_xlabel(f"Round {round}", fontsize=20)
                # if round == 5:
                #     target_ax.set_xlabel(f"Rounds", fontsize=30)
            if round > 1:
                target_ax.axes.get_yaxis().set_visible(False)
                if i < 2:
                    target_ax.axes.get_xaxis().set_visible(False)

    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.1, hspace=0.1)

    plt.savefig(plot_name, dpi=300, bbox_inches='tight')
    plt.close()
