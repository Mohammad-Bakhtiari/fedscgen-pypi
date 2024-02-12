# FedScGen: Federated batch effect correction using single-cell generative models

## Data reproducibility
The data used in the paper is available at the following link: [https://www.dropbox.com/s/ ]
The data was taken from []

To repoduce the results of the paper, please follow the instructions:
1- setup a virtual environment using conda:
```bash
conda env create -f environment.yml
conda activate fedscgen
cd /experiments
```

2- Run the [experiments](/experiments) in the following recommended order:
- Run `./run_scgen.sh` to train the scgen model and correct the data.
- Run `./run_fedscgen.sh` to train the fedscgen model and correct the data.
- Run `./scgen-with-batch-out.sh` to train the scgen model without the left-out  batch.
- Run `./cent-scgen-with-batch-out.sh "dataset4" "HumanPancreas.h5ad" "" false false "1"` to reproduce Figure 1b.
- Run `./run_pca.sh` to reduce dimensionality of corrected and raw datasets
- Run `./run-classification.sh` to train the classifier and evaluate the performance of the corrected and raw datasets.
- 

## License
This project is licensed under the Apache License 2.0 - see the [LICENSE](LICENSE) file for details.
