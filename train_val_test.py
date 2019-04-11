#!/usr/bin/env python
# encoding: utf-8

import numpy as np, argparse, sys, h5py
from os.path import join, abspath, dirname, exists
from keras.callbacks import ModelCheckpoint, EarlyStopping
from keras.optimizers import Adadelta
from keras.models import load_model
from models import DeepCNN, DeepBind
import os
os.environ["CUDA_VISIBLE_DEVICES"] = "0"

from utils import *


def parse_args():
    parser = argparse.ArgumentParser(description="Convert sequence and target")
    parser.add_argument('-datadir', dest='datadir', type=str, help='Positive data for training, testing')
    parser.add_argument('-k', dest='k_folds', type=int, default=5, help='k-folds cross-validation')
    parser.add_argument('-batchsize', dest='batchsize', type=int, default=300, help='the size of one batch')
    parser.add_argument('-ratio', dest='ratio', type=float, default=0.2, help='the propotion of validation data over the training data')
    parser.add_argument('-params', dest='params', type=int, default=1, help='the number of paramter settings')
    parser.add_argument('--train', dest='train', action='store_true', help='train and test step')
    parser.add_argument('--no-train', dest='train', action='store_false', help='only test step.')
    parser.add_argument('-plot', dest='plot', action='store_true', default=True, help='only test step')
    parser.add_argument('-run', dest='run', type=str, default='onehot', help='three encoding methods, including onehot, kmer, gappedkmer')
    
    return parser.parse_args()

def main():

    file_path = dirname(abspath(__file__))
    args = parse_args()
    
    # print the current dataset name
    name = (args.datadir).split('/')[-2]
    print 'working on %s now' % name
    
    # convert raw data 
    if args.run == 'onehot':
       print 'only using one-hot vectors.'
       seq_path_train = args.datadir + '/train.hdf5'
       with h5py.File(seq_path_train, 'r') as f:
           seqs_data_train = np.asarray(f['data'])
           intensity_train = np.asarray(f['intensity'])
    else:
       print >> sys.stderr, 'invalid command!';sys.exit(1)
    seqs_num = seqs_data_train.shape[0]; seqs_len = seqs_data_train.shape[1]; seqs_dim = seqs_data_train.shape[2]
    print 'there are %d seqences, each of which is a %d*%d array' %(seqs_num, seqs_len, seqs_dim)
    input_shape = (seqs_len, seqs_dim)
    
    # k-folds cross-validation
    indices = np.arange(seqs_num)
    np.random.shuffle(indices)
    seqs_data_train = seqs_data_train[indices]
    intensity_train = intensity_train[indices]
    train_ids, test_ids, valid_ids = Id_k_folds(seqs_num, args.k_folds, args.ratio)

    R2 = []
    
    if not exists(file_path + '/model/%s' % name):
       print 'Building' + file_path + '/model/%s' % name
       os.mkdir(file_path + '/model/%s' % name)
    f_params = open(file_path + '/model/%s/params.txt' % name, 'w')
    for fold in range(args.k_folds):        
        x_train = seqs_data_train[train_ids[fold]]
        y_train = intensity_train[train_ids[fold]]
        
        x_test = seqs_data_train[test_ids[fold]]
        y_test = intensity_train[test_ids[fold]]
        
        x_valid = seqs_data_train[valid_ids[fold]]
        y_valid = intensity_train[valid_ids[fold]]
        
        if args.train:
           history_all = {}
           for params_num in range(args.params):
               params = RandomSample()
               print >> f_params, "the {}-th paramter setting of the {}-th folds is {}".format(params_num, fold, params)
               
               print 'Building model...'
               model = DeepBind(input_shape, params)
               #return
               checkpointer = ModelCheckpoint(filepath=file_path + '/model/%s/params%d_bestmodel_%dfold.hdf5' 
                                                % (name, params_num, fold), monitor='val_loss', verbose=1, save_best_only=True)
               earlystopper = EarlyStopping(monitor='val_loss', patience=15, verbose=1)

               print 'Training model...'
               myoptimizer = Adadelta(epsilon=params['DELTA'], rho=params['MOMENT'])
               model.compile(loss='mean_squared_error', optimizer=myoptimizer)
               History = model.fit(x_train, y_train, epochs=100, batch_size=args.batchsize, shuffle=True,
                                   validation_data=(x_valid, y_valid), callbacks=[checkpointer, earlystopper], verbose=2)
               history_all[str(params_num)] = History
           best_num = SelectBest(history_all, file_path + '/model/%s/' % name, fold, 'val_loss')
           if args.plot:
              PlotandSave(history_all[str(best_num)], file_path + '/model/%s/figure_%dfold.png' % (name, fold), fold, 'val_loss')
        print >> f_params, "\n\n" 
        f_params.flush() 
        print 'Testing model...'
        # load_model('')
        model.load_weights(file_path + '/model/%s/params%d_bestmodel_%dfold.hdf5' % (name, best_num, fold))
        results = model.evaluate(x_test, y_test)
        print results
        y_pred = model.predict(x_test, batch_size=args.batchsize, verbose=1)
        y_pred = np.asarray([y[0] for y in y_pred])
        y_real = np.asarray([y[0] for y in y_test])
        with open(file_path + '/model/%s/score_%dfold.txt' % (name, fold), 'w') as f:
           assert len(y_pred) == len(y_real), 'dismathed!'
           for i in range(len(y_pred)):
               print >> f, '{:.4f} {}'.format(y_pred[i], y_real[i])

        print 'Calculating R2...'
        coeff = ComputePCC(y_pred, y_real)
        R2.append(coeff)
       
    f_params.close()
    print "the mean R2 is {}.".format(np.mean(R2))
    outfile = file_path + '/model/%s/metrics.txt' % (name)
    with open(outfile,'w') as f:
        for i in range(len(R2)):
            print >> f, "{:.4f}".format(R2[i])
        print >> f, "{:.4f}".format(np.mean(R2))
    
if __name__ == '__main__': main()    
    
    
    
    
    
    
