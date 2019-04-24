# Installation
## Download data
## Install conda
## Install kipoi DeepSEA model
http://kipoi.org/models/DeepSEA/predict/
```
pip install kipoi
kipoi env create DeepSEA/predict # add --gpu to install gpu-compatible deps
source activate kipoi-DeepSEA__predict # or kipoi-gpu-DeepSEA__predict
kipoi env install DeepSEA/predict
```
## Install libraries
```
pip install biopython
```
## Run initial scripts to generate important files
```
bash extract_expression.bash
```
## Run main script :
```
bash output_predictions_500kb_optimized_TF_big_batch.bash
```