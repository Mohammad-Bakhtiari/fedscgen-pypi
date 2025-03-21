approach=$1
dataset=$2
GPU=$3
n_clients=$4

echo "Args: approach: $approach, dataset: $dataset, GPU: $GPU, n_clients: $n_clients"

root_dir="$(dirname "$PWD")"
output="${root_dir}/results/${approach}/${dataset}/all"
if [ "$approach" == "scgen" ]; then
  adata="${output}/corrected.h5ad"
  output="${output}/classification"
else
  adata="${output}/BO0-C${n_clients}/fed_corrected.h5ad"
  output="${output}/BO0-C${n_clients}/classification"
fi
if [[ ! -f "$adata" ]]; then
  echo "ERROR: Missing adata file: $adata"
  exit 1
fi
export CUBLAS_WORKSPACE_CONFIG=:4096:8
python "${root_dir}/scripts/classification.py" \
      --adata "$adata" \
      --output "$output" \
      --cell_key "cell_type" \
      --batch_key "batch" \
      --epoch 50 \
      --lr 0.001 \
      --batch_size 128 \
      --init_model_path "${root_dir}/models/classification/${dataset}/model.pth" \
      --norm_method "min_max" \
      --model "mlp-norm" \
      --hidden_size "800,200" \
      --batch_out "1" \
      --gpu "$GPU"
