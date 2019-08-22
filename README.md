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
+ Usage: you can excute run.sh script directly, in which you should modify python command accordingly, e.g.:
  ```
  python train_val_test.py -datadir ./pbmdata/$eachTF/data -run 'noshape' -model 'shallow' -batchsize 300 -k 5 -params 30 --train
  ```
 The command '-model' can be a choice of {'shollow', 'deep'}, where 'shollow' means DeepBind_K, and 'deep' means DeepCNN.
 
**Run DLBSS(shallow) or DLBSS(deep) using DNA shape information**
+ Usage: you can excute run.sh script directly, in which you should modify python command accordingly, e.g.:
  ```
  python train_val_test_hybrid.py -datadir ./pbmdata/$eachTF/data -run 'shape' -model 'shallow' -batchsize 300 -k 5 -params 30 --train
  ```
The command '-run' can be a choice of {'shape', 'MGW', 'ProT', 'Roll', 'HelT'}, where 'shape' means using all shape feature, 'MGW' means using MGW shape feature, and so on.<br />
The command '-model' can be a choice of {'shollow', 'deep'}, where 'shollow' means DLBSS(shallow), and 'deep' means 'DLBSS(deep)'.<br />
**Note that** you should change the ouput path in the run.sh script, the naming rule is: 'model_' + args.model + '_' + args.run.

+ Type the following for details on other optional arguments:
	```
  python train_val_test_hybrid.py -h
	```
