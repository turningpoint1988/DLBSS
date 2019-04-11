# DLBSS
Ｉmplementation of "predicting in vitro transcription factor binding sites using DNA sequence + shape"

## Data preparation
Firstly, using encode.py script to preprocess DNA sequences and their corresponding shape features.
+ Usage:
  ```
  bash encode.sh <pbmdata>
  ```
  <pbmdata> denotes the path of storing experimental data.

## Run 
##Run the models without using DNA shape information##
+ Usage:
  ```
  python train_val_test.py -datadir <data path> -run 'onehot' -batchsize 300 -k 5 -params 30 --train
  ```
  <data path> denotes the path of the current dataset.
 
##Run the models using DNA shape information##
+ Usage:
  ```
  python train_val_test_hybrid.py -datadir <data path> -run 'MPRH' -batchsize 300 -k 5 -params 30 --train
  ```
  <data path> denotes the path of the current dataset.

+ Type the following for details on other optional arguments:
	```
  python train_val_test_hybrid.py -h
	```
