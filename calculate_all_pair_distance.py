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
    for i in xrange(len(vps)):
        for j in xrange(i+1,len(vps)):
            print i,j,utils.hammingDistance(vps[i],vps[j])


if __name__ == '__main__':
    options = utils.getOptions()
    calculateAllPairDistance(options)
