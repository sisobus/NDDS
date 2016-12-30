#!/usr/bin/python
#-*- coding:utf-8 -*-
import sys
import utils
import string

def calculateAllPairDistance(options):
    numberOfData = options['numberOfData']
    dim          = options['numberOfDimension']
    numberOfVP   = options['numberOfVP']
    cardinality  = options['numberOfAlphabet']
    typeOfVP     = options['typeOfVP']
    vps = utils.readDataFromFile(utils.getVPFileName(options))
    d = [ 0 for i in xrange(dim+1) ]
    s = set()
    for i in xrange(len(vps)):
        for j in xrange(i+1,len(vps)):
            dist = utils.hammingDistance(vps[i],vps[j])
            d[dist] = d[dist] + 1
            s.add(dist)
            print i,j,dist
    for i in xrange(0,dim+1):
        print i,d[i]
    print len(s)


if __name__ == '__main__':
    options = utils.getOptions()
    calculateAllPairDistance(options)
