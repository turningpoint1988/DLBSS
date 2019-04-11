#/usr/bin/python

import sys, os
from os.path import exists

dir_path = sys.argv[1]

if exists(dir_path + '/R2.txt'):
   os.remove(dir_path + '/R2.txt')

file_names = os.listdir(dir_path)
f_out = open(dir_path + '/R2.txt', 'w')

for name in file_names:
    with open(dir_path + '/%s/metrics.txt' % name, 'r') as f:
        lines = f.readlines()
    line = lines[-1].strip()

    print >> f_out, "{} {}".format(name, line)
    
f_out.close()
