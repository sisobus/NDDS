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

def generateVantagePointsWithPattern(options):
    numberOfData = options['numberOfData']
    dim          = options['numberOfDimension']
    numberOfVP   = options['numberOfVP']
    cardinality  = options['numberOfAlphabet']
    typeOfVP     = options['typeOfVP']
    datas = utils.getDataInFile(utils.getDataFileName(options))
    majorPattern = [ [] for i in xrange(numberOfVP) ]
    minorPattern = [ [] for i in xrange(numberOfVP) ]
    for j in xrange(dim):
        d = {}
        for i in xrange(numberOfData):
            if datas[i][j] in d:
                d[datas[i][j]] += 1
            else :
                d[datas[i][j]] = 0
        d = sorted(d.items(), key=lambda x: x[1], reverse=True)
        for k in xrange(numberOfVP):
            majorPattern[k].append(d[k][0])
            minorPattern[k].append(d[-1-k][0])
    utils.writeDataToFile('vp/vp_%d_%d_%d_%s.txt'%(dim,numberOfVP,cardinality,'major'),majorPattern)
    utils.writeDataToFile('vp/vp_%d_%d_%d_%s.txt'%(dim,numberOfVP,cardinality,'minor'),minorPattern)

def calculateMany(vps):
    s = set()
    for i in xrange(len(vps)):
        for j in xrange(len(vps)):
            if i != j :
                now = utils.hammingDistance(vps[i],vps[j])
                s.add(now)
    return s


def ok(dist,dim):
    for i in xrange(1,dim+1):
        if dist[i] == 0:
            return False
    return True

def generateVantagePointsWithManyAlgorithm(options):
    numberOfData = options['numberOfData']
    dim          = options['numberOfDimension']
    numberOfVP   = options['numberOfVP']
    cardinality  = options['numberOfAlphabet']
    typeOfVP     = options['typeOfVP']

    vps = []
    vps.append(generateUniformRandomVP(dim,cardinality))
    dist = [ 0 for i in xrange(dim+1) ]
    while len(vps) < numberOfVP:
        cur_vp = generateUniformRandomVP(dim,cardinality)
        vps.append(cur_vp)
        if ok(dist,dim):
            continue
        ss = calculateMany(vps)
        nnext = False
        print len(vps),
        for sss in ss:
            if dist[sss] != 0:
                nnext = True
        if nnext:
            vps = vps[:-1]
            continue
        for sss in ss:
            dist[sss] += 1
        

    utils.writeDataToFile('vp/vp_%d_%d_%d_%s.txt'%(dim,numberOfVP,cardinality,typeOfVP),vps)




#utils.writeDataToFile('vp/vp_%d_%d_%d_%s.txt'%(dim,numberOfVP,cardinality,typeOfVP),vps)

if __name__ == '__main__':
    options = utils.getOptions()

    utils.createDirectory('vp')
    generateVantagePointsWithManyAlgorithm(options)
