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

def generateAllRandomVantagePoints(options):
    dim          = options['numberOfDimension']
    numberOfVP   = options['numberOfVP']
    cardinality  = options['numberOfAlphabet']
    typeOfVP     = options['typeOfVP']

    vps = []
    for i in xrange(numberOfVP):
        vps.append(generateUniformRandomVP(dim,cardinality))

    utils.writeDataToFile('vp/vp_%d_%d_%d_%s.txt'%(dim,numberOfVP,cardinality,typeOfVP),vps)

def generateHeuristicVantagePoints(options):
    dim          = options['numberOfDimension']
    numberOfVP   = options['numberOfVP']
    cardinality  = options['numberOfAlphabet']
    typeOfVP     = options['typeOfVP']

    threshold = dim*0.4

    vps = []
    for i in xrange(numberOfVP):
        if i == 0:
            vps.append(generateUniformRandomVP(dim,cardinality))
            continue
        while True:
            nvp = generateUniformRandomVP(dim,cardinality)
            ok = True
            for j in xrange(i):
                dist = utils.hammingDistance(vps[j],nvp)
                if dim-dist > threshold:
                    ok = False
                    break
            if ok:
                vps.append(nvp)
                break
    utils.writeDataToFile('vp/vp_%d_%d_%d_%s.txt'%(dim,numberOfVP,cardinality,typeOfVP),vps)

cornerPoints = []
def getCornerPoints(pos,vp,dim,alphabet,cardinality):
    if pos == dim:
        t = []
        for i in xrange(len(vp)):
            t.append(vp[i])
        cornerPoints.append(t)
        return
    for i in xrange(cardinality):
        vp.append(alphabet[i])
        getCornerPoints(pos+1,vp,dim,alphabet,cardinality)
        vp.pop(pos)

def getTotalDistance(vps,vp):
    tvps = []
    for i in xrange(len(vps)):
        t = []
        for j in xrange(len(vps[i])):
            t.append(vps[i][j])
        tvps.append(t)
    t = []
    for i in xrange(len(vp)):
        t.append(vp[i])
    tvps.append(t)

    s = 0.0
    sSquare = 0.0
    n = 0
    for i in xrange(len(tvps)):
        for j in xrange(len(tvps)):
            if i == j:
                continue
            dist = utils.hammingDistance(tvps[i],tvps[j])
            s = s + dist
            sSquare = sSquare + (dist * dist)
            n = n + 1

    m = float(s) / n
    
    ret =  (float(sSquare) / n) - (m*m)
    return float(ret)

def getTotalCostFunction(vps,vp,base):
    tvps = []
    for i in xrange(len(vps)):
        t = []
        for j in xrange(len(vps[i])):
            t.append(vps[i][j])
        tvps.append(t)
    t = []
    for i in xrange(len(vp)):
        t.append(vp[i])
    tvps.append(t)

    s = 0.0
    sSquare = 0.0
    n = 0
    
    for i in xrange(len(tvps)):
        for j in xrange(len(tvps)):
            if i == j:
                continue
            dist = abs(base-float(utils.hammingDistance(tvps[i],tvps[j])))
            s = s + dist
            sSquare = sSquare + (dist * dist)
            n = n + 1

    return s
    m = float(s) / n

    ret =  (float(sSquare) / n) - (m*m)
    return float(ret)

def generateGreedyVantagePoints(options):
    alphabet = list(string.ascii_uppercase)
    dim          = options['numberOfDimension']
    numberOfVP   = options['numberOfVP']
    cardinality  = options['numberOfAlphabet']
    base         = abs(float(dim)-float(dim)/float(cardinality))

    getCornerPoints(0,[],dim,alphabet,cardinality)
    vps = []
    vps.append([ 'A' for i in xrange(dim) ])
    #vps.append([ 'B' for i in xrange(dim) ])
    #vps.append([ 'C' for i in xrange(dim) ])
    #vps.append([ 'D' for i in xrange(dim) ])
    for i in xrange(numberOfVP - 1):
        print i
        #mx, mx_idx = (0,0)
        mn, mn_idx = (987654321.0,0)
        for j in xrange(len(cornerPoints)):
            variation = getTotalCostFunction(vps,cornerPoints[j],base)
            if mn > variation:
                mn, mn_idx = (variation,j)
            #if mx < dist:
            #    mx, mx_idx = (dist,j)
        print mn,mn_idx,cornerPoints[mn_idx]
        vps.append(cornerPoints[mn_idx])
    utils.writeDataToFile('vp/vp_%d_%d_%d_greedy.txt'%(dim,numberOfVP,cardinality),vps)


if __name__ == '__main__':
    options = utils.getOptions()

    utils.createDirectory('vp')
    #generateVantagePoints(options)
    #generateGreedyVantagePoints(options)
    generateAllRandomVantagePoints(options)
#generateHeuristicVantagePoints(options)
