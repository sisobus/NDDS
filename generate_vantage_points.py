#!/usr/bin/python
#-*- coding:utf-8 -*-
import utils
import string
import random

def generateUniformRandomVP(numberOfDimension,numberOfAlphabet):
    ret = []
    alphabet = list(string.ascii_uppercase)
    for i in xrange(numberOfDimension):
        cur = random.randrange(0,numberOfAlphabet)
        ret.append(alphabet[cur])
    return ret

def generateVantagePoints(options):
    dim          = options['numberOfDimension']
    numberOfVP   = options['numberOfVP']
    cardinality  = options['numberOfAlphabet']

    vps = []
    for i in xrange(numberOfVP):
        vps.append(generateUniformRandomVP(dim,cardinality))
    
    utils.writeDataToFile('vp/vp_%d_%d_%d_random.txt'%(dim,numberOfVP,cardinality),vps)

if __name__ == '__main__':
    options = utils.getOptions()

    utils.createDirectory('vp')
    generateVantagePoints(options)
    
