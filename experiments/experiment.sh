#!/bin/bash

chmod +x run_scgen.sh
echo "Running scgen"
./run_scgen.sh

chmod +x run-tuning.sh
echo "Running tuning"
./run-tuning.sh

chmod +x run_fedscgen.sh
echo "Running fedscgen"
./run_fedscgen.sh

chmod +x scgen-with-batch-out.sh
echo "Running scgen with batch out"
./scgen-with-batch-out.sh 'HumanPancreas.h5ad' '' false false "0,1,2,3,4"



./run-classification.sh