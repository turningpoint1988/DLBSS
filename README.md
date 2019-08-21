# DLBSS
Implementation of "predicting in vitro transcription factor binding sites using DNA sequence + shape"

## Requirements

+ To install Keras with Tensorflow backend, please refer to https://keras.io/#installation. 

+ Python 2.7

## Data preparation
Firstly, using encode.sh script to preprocess DNA sequences and their corresponding shape features.
+ Usage:
  ```
  bash encode.sh <pbmdata>
  ```
  **'pbmdata'** denotes the path of storing experimental data, e.g. /yourpath/pbmdata.

## Run 
**Run DeepBind_K or DeepCNN without using DNA shape information**
+ Usage: you can excute run.sh script directly, in which you should modify python command in it accordingly, e.g.:
  ```
  python train_val_test.py -datadir <data path> -run 'noshape' -model 'shallow' -batchsize 300 -k 5 -params 30 --train
  ```
  **'data path'** denotes the path of the current dataset, e.g.  the command '-run' can be a choice of {'shollow', 'deep'}, where 'shollow' means DeepBind_K, and 'deep' means DeepCNN
 
**Run DLBSS(shallow) or DLBSS(deep) using DNA shape information**
+ Usage: you can excute run.sh script directly, in which you should modify python command in it accordingly, e.g.:
  ```
  python train_val_test_hybrid.py -datadir <data path> -run 'shape' -model 'shallow' -batchsize 300 -k 5 -params 30 --train
  ```
  **data path** denotes the path of the current dataset.

+ Type the following for details on other optional arguments:
	```
  python train_val_test_hybrid.py -h
	```
