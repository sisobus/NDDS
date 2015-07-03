#!/usr/bin/python
#-*- coding:utf-8 -*-
import sys
import random
import os
import string
import commands
import glob
import utils

def generateUniformData(options):
    alphabet     = list(string.ascii_uppercase)
    dim          = options['numberOfDimension']
    size         = options['numberOfData']
    cardinality  = options['numberOfAlphabet']

    dataFileName = 'data/data_%d_%d_u_%d.txt'%(size,dim,cardinality)
    if os.path.exists(dataFileName):
        return
    cur = [ [] for i in xrange(size) ]
    for i in xrange(size):
        t = []
        for j in xrange(dim):
            t.append(alphabet[random.randrange(0,cardinality)])
        cur[i] = t
    with open(dataFileName,'w') as fp:
        fp.write('%d %d 1\n'%(size,dim))
        for c in cur:
            fp.write('0000-000000:')
            fp.write(str(c[0]))
            for j in xrange(1,len(c)):
                fp.write(','+str(c[j]))
            fp.write('\n')
    print dataFileName


def generateClusteredData(options):
    alphabet     = list(string.ascii_uppercase)
    dim          = options['numberOfDimension']
    size         = options['numberOfData']
    cardinality  = options['numberOfAlphabet']
    cluster      = int(options['distribution'][1:])
    dataFileName = 'data/data_%d_%d_%d_%d.txt'%(size,dim,cluster,cardinality)

    command = 'perl data/gdp.pl %d %d %d %d %s'%(size,cluster,dim,cardinality,dataFileName)
    print command
    utils.executeCommand(command)
    print dataFileName


def generateData(options):
    if options['distribution'] == 'u':
        generateUniformData(options)
    elif options['distribution'][0] == 'c':
        generateClusteredData(options)

if __name__ == '__main__':
    options = utils.getOptions()

    utils.createDirectory('data')
    generateData(options)    
