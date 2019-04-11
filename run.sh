#!/usr/bin/bash


# In the training and testing phase
for eachTF in `ls ./pbmdata/`
do 
	echo $eachTF
	if [ -d ./model_MGW/$eachTF ]; then
	   echo $eachTF 'has existed.'
	   continue
	fi
    # run: 'MGW', 'ProT', 'Roll', 'HelT', 'MPRH'
    python train_val_test_hybrid.py -datadir ./pbmdata/$eachTF/data -run 'MGW' -batchsize 300 -k 5 -params 30 --train
done

