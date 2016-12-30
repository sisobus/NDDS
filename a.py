#!/usr/bin/python
#-*- coding:utf-8 -*-
import glob
import utils

filenames = []
for i in xrange(20,121,20):
    filenames.append('result/nd_1000000_%d_u.out'%(i))

for filename in filenames:
    print filename
    with open(filename,'r') as fp:
        r = fp.read().split('\n')
    for line in r:
        if line.find('BestChild') <> -1:
            print line
