# FedscGen: privacy-aware federated batch effect correction of single-cell RNA sequencing data

## 🔧 Setup Environment

Before reproducing the results, set up the required environments.

To create the Conda environments for running FedscGen, please follow the instructions in the [`setup_env.sh`](setup_env.sh) script.

⚠️ Make sure you execute this script in a Bash shell.
```shell
chmod +x setup_env.sh
./setup_env.sh
```
It will create two Conda environments:
1. `fedscgen`: Python environment for FedscGen.
2. `r_eval`: R + Python environment for benchmarking.


## Dataset and Models
For reproducibility, please ensure the preprocessed datasets are downloaded and extracted to the [`data/datasets`](data/datasets) directory. Optionally, the initial models can also be downloaded to the [`models/`](models) directory.

* The initial PyTorch models are available at <a href="https://doi.org/10.5281/zenodo.11505054"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.11505054.svg" alt="DOI"></a>

* All preprocessed datasets used in the paper are available at <a href="https://doi.org/10.5281/zenodo.11489844"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.11489844.svg" alt="DOI"></a>

## 📊 Reproduce Results
All models are initialized using a fixed seed for reproducibility.

Navigate to the [`experiments/`](experiments) directory and run [`experiments.sh`](experiments/experiment.sh)
while providing a comma-separated list of GPU indices to use for training. For example, to use GPUs 0 and 1:

```shell
cd experiments
chmod +x experiment.sh
./experiment.sh 0,1
```

Once the experiments are complete, run the evaluation metrics by navigating to the [`metrics/`](metrics) directory and executing [`evaluate.sh`](metrics/evaluate.sh):
```shell
cd metrics
chmod +x evaluation.sh
./evaluation.sh
```
All results will be saved in the [`results/`](results) directory. Both scripts will automatically activate the appropriate Conda environment and make full use of available system resources.

## <a href="https://featurecloud.ai/app/fedscgen" target="_blank"> <img src="https://featurecloud.ai/assets/fc_logo.svg" alt="FeatureCloud Logo" width="160"/> </a> app
FedscGen is implemented for real-world federated collaboration as a FeatureCloud app with automated deployment.

* Explore the app: [FedscGen App on FeatureCloud](https://featurecloud.ai/app/fedscgen)
* Source code: [fc-fedscgen GitHub Repository](https://github.com/Mohammad-Bakhtiari/fc-fedscgen)

## License
This project is licensed under the Apache License 2.0 - see the [LICENSE](LICENSE) file for details.
