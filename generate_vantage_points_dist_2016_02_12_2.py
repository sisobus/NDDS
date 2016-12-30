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

def calculateMany(vps):
    s = set()
    if len(vps) == 0:
        return 0,0
    d = [ 0 for i in xrange(len(vps[0])+1) ]
    for i in xrange(len(vps)):
        for j in xrange(len(vps)):
            if i != j :
                now = utils.hammingDistance(vps[i],vps[j])
                d[now] += 1
                s.add(now)
    return len(s),max(d)


def generateVantagePointsWithManyAlgorithm(options):
    numberOfData = options['numberOfData']
    dim          = options['numberOfDimension']
    numberOfVP   = options['numberOfVP']
    cardinality  = options['numberOfAlphabet']
    typeOfVP     = options['typeOfVP']
    datas = utils.getDataInFile(utils.getDataFileName(options))
    threshold = 2

    majorPattern = [] 
    for j in xrange(dim):
        d = {}
        for i in xrange(numberOfData):
            if datas[i][j] in d:
                d[datas[i][j]] += 1
            else :
                d[datas[i][j]] = 0
        d = sorted(d.items(), key=lambda x: x[1], reverse=True)
        majorPattern.append(d[0][0])

    vps = []
    vps.append(majorPattern)
    d = [ 0 for i in xrange(dim+1) ]
    while len(vps) < numberOfVP:
        print len(vps)
        print d
        newPoint = generateUniformRandomVP(dim,cardinality)
        ok = False
        for i in xrange(len(vps)):
            dist = utils.hammingDistance(vps[i],newPoint)
            if d[dist] == 0:
                ok = True
            if d[dist] > 3:
                ok = False
                break
        if ok:
            for i in xrange(len(vps)):
                dist = utils.hammingDistance(vps[i],newPoint)
                d[dist] += 1
            vps.append(newPoint)


#    utils.writeDataToFile('vp/vp_%d_%d_%d_%s.txt'%(dim,numberOfVP,cardinality,typeOfVP),vps[:numberOfVP])

if __name__ == '__main__':
    options = utils.getOptions()

    utils.createDirectory('vp')
    generateVantagePointsWithManyAlgorithm(options)
