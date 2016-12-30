#!/usr/bin/python
#-*- coding:utf-8 -*-
import utils
import string
import random

options = utils.getOptions()

def getNewAlphabet(cardinality):
    alphabet = list(string.ascii_uppercase)
    return alphabet[random.randint(0,cardinality-1)]

def getMajorPattern():
    numberOfData = options['numberOfData']
    dim          = options['numberOfDimension']
    numberOfVP   = options['numberOfVP']
    cardinality  = options['numberOfAlphabet']
    typeOfVP     = options['typeOfVP']
    datas = utils.getDataInFile(utils.getDataFileName(options))
    rs = [ random.randint(0,numberOfData-1) for i in xrange(100) ]
    majorPattern = ''
    for j in xrange(dim):
        d = {}
        for i in rs:
            if datas[i][j] in d:
                d[datas[i][j]] += 1
            else :
                d[datas[i][j]] = 0
        d = sorted(d.items(), key=lambda x: x[1], reverse=True)
        majorPattern += d[0][0]
    return majorPattern

def generatePoint(vp,N,newDist):
    newPoint = vp
    changePosition = set()
    while len(changePosition) < newDist:
        pos = random.randint(0,N-1)
        changePosition.add(pos)
    for position in changePosition:
        na = vp[position]
        while na == vp[position]:
            na = getNewAlphabet(options['numberOfAlphabet']) 
        newPoint = newPoint[:position]+na+newPoint[position+1:]
    return newPoint

def generateCandidateSet(vp,N,lVl,K):
    newCandidateSet = set()
    dist = set()
    while len(dist) < (K)/lVl:
        print len(dist) 
        newDist = random.randint(1,N)
        if not newDist in dist:
            dist.add(newDist)
            newCandidateSet.add((newDist,generatePoint(vp,N,newDist)))
    return newCandidateSet

def eraseCandidateSet(candidateSet,dist):
    newCandidateSet = set()
    for candidate in candidateSet:
        if not candidate[0] in dist:
            newCandidateSet.add(candidate)
    return newCandidateSet

def _algorithm(N,M,K):
    V = set()
    dist = set()
    candidateSet = set()
    V.add(getMajorPattern())
    while len(V) < M:
        for vp in V:
            candidateSet.add(generateCandidateSet(vp,N,len(V),K))
        maxDistCadinality = 0
        nextPoint = ()
        nextDist = set()
        for candidate in candidateSet:
            newDist = set()
            for vp in V:
                newDist = newDist.add(utils.hammingDistance(vp,candidate))
            if len(newDist) > maxDistCadinality:
                maxDistCadinality = len(newDist)
                nextPoint = candidate
                nextDist = newDist
        V.add(nextPoint[1])
        dist.add(newDist)
        candidateSet = eraseCandidateSet(candidateSet)
    return V

def algorithm():
    dim          = options['numberOfDimension']
    numberOfVP   = options['numberOfVP']
    cardinality  = options['numberOfAlphabet']
    return _algorithm(dim,numberOfVP,dim)

if __name__ == '__main__':
    utils.createDirectory('vp')
    
    print algorithm()
