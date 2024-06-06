# FedScGen: Federated batch effect correction using single-cell generative models

## Reproducibility
The data used in the paper is available at the following link: [https://www.dropbox.com/s/ ]
The data was taken from []

To reproduce the results of the paper, please follow the instructions:
1- setup a virtual environment using conda for running FedscGen:
```bash
conda env create -f environment.yml
conda activate fedscgen
cd /experiments
```
2- Run the [experiment.sh](/experiments/experiment.sh) which includes the following experiments to reproduce the results of the paper:
- Run `./run_scgen.sh` to train the scgen model and correct the data.
- Run `./run_fedscgen.sh` to train the fedscgen model and correct the data.
- Run `./scgen-with-batch-out.sh "dataset4" "HumanPancreas.h5ad" "" false false` to train the scgen model without the left-out  batch.
- Run `./run-classification.sh` to train the classifier and evaluate the performance of the corrected and raw datasets.
- 

3- setup a `R` virtual environment using conda for running kBET and LISI for evaluation:
```bash
conda env create -f r_eval.yml
conda activate r_eval
Rscript install_packages.R
cd /metrics
```
 4- Run the [evaluation.sh](/metrics/evaluation.sh) which includes the following experiments to evaluate the performance of the corrected data.


## Dataset and Models
The datasets used in the paper are available at [https://doi.org/10.5281/zenodo.11489844](https://doi.org/10.5281/zenodo.11489844)
The initial PyTorch models are available at [ https://doi.org/10.5281/zenodo.11505054]( https://doi.org/10.5281/zenodo.11505054)

## License
This project is licensed under the Apache License 2.0 - see the [LICENSE](LICENSE) file for details.
