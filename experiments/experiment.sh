#!/bin/bash
AVAILABLE_GPUS="${1:-0,1,2,3}"

source "$(conda info --base)/etc/profile.d/conda.sh" || {
    echo "❌ Failed to source conda environment script."
    exit 1
}

# Activate fedscgen environment
conda activate fedscgen || {
    echo "❌ Failed to activate fedscgen environment."
    exit 1
}

chmod +x run-scgen.sh
echo "Running scgen for combined and dropped inclusion scenarios"
./run-scgen.sh "${AVAILABLE_GPUS}"

chmod +x scgen-all.sh
echo "Running scGen with different seeds for all inclusion scenario"
./scgen-all.sh "${AVAILABLE_GPUS}"

chmod +x scgen-with-batch-out.sh
echo "Running scgen on HumanPancreas with batch out"
./scgen-with-batch-out.sh 'HumanPancreas.h5ad' '' false false "0,1,2,3,4"

chmod +x run-tuning.sh
echo "Running Hyper-parameters tuning for FedscGen without SMPC"
./run-tuning.sh "${AVAILABLE_GPUS}"

chmod +x run-fedscgen.sh
echo "Running FedscGen for every inclusion scenarios"
./run-fedscgen.sh  "${AVAILABLE_GPUS}"

chmod +x run-fedscgen-smpc.sh
echo "Running FedscGen-SMPC with different seeds for the 'all' inclusion scenario"
./run-fedscgen-smpc.sh "${AVAILABLE_GPUS}"

chmod +x run-classification.sh
echo "Running centralized classification using corrected data"
./run-classification.sh "${AVAILABLE_GPUS}"