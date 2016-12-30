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
    d = [ 0 for i in xrange(len(vps[0])+1) ]
    for i in xrange(len(vps)):
        for j in xrange(len(vps)):
            if i != j :
                now = utils.hammingDistance(vps[i],vps[j])
                d[now] += 1
                s.add(now)
    return len(s),d

def generateRandomVP(options):
    numberOfData = options['numberOfData']
    dim          = options['numberOfDimension']
    numberOfVP   = options['numberOfVP']
    cardinality  = options['numberOfAlphabet']
    typeOfVP     = options['typeOfVP']
    datas = utils.getDataInFile(utils.getDataFileName(options))
    vps = []
    for i in xrange(numberOfVP):
        vps.append(datas[random.randint(0,numberOfData)])

    utils.writeDataToFile('vp/vp_%d_%d_%d_%s.txt'%(dim,numberOfVP,cardinality,typeOfVP),vps)



def generateVantagePointsWithManyAlgorithm(options):
    numberOfData = options['numberOfData']
    dim          = options['numberOfDimension']
    numberOfVP   = options['numberOfVP']
    cardinality  = options['numberOfAlphabet']
    typeOfVP     = options['typeOfVP']
    datas = utils.getDataInFile(utils.getDataFileName(options))

    majorPattern = [ [] for i in xrange(numberOfVP+1) ]
    for j in xrange(dim):
        d = {}
        for i in xrange(numberOfData):
            if datas[i][j] in d:
                d[datas[i][j]] += 1
            else :
                d[datas[i][j]] = 0
        d = sorted(d.items(), key=lambda x: x[1], reverse=True)
        for k in xrange(1):
            majorPattern[k].append(d[k][0])

    vps = []
    vps.append(majorPattern[0])
    d = [ 0 for i in xrange(dim+1) ]
    one_pass = True
    while len(vps) < numberOfVP:
        print len(vps)
        ans,ansDataIndex = -1,-1
        if one_pass:
            for i in xrange(len(datas)):
                ok = False
                for j in xrange(len(vps)):
                    dist = utils.hammingDistance(datas[i],vps[j])
                    if d[dist] > 1:
                        ok = True
                    if datas[i] == vps[j]:
                        ok = True
                if ok:
                    continue
                for j in xrange(len(vps)):
                    dist = utils.hammingDistance(datas[i],vps[j])
                    d[dist] += 1
                vps.append(datas[i])
            one_pass = False
        else :
            ans,ansDataIndex = -1,-1
            fx = 987654321
            for i in xrange(len(datas)):
                ok = False
                for j in xrange(len(vps)):
                    if datas[i] == vps[j]:
                        ok = True
                if ok:
                    continue
                vps.append(datas[i])
                cur,dists = calculateMany(vps)
                if cur > ans:
                    ans,ansDataIndex = cur,i
                    fx = max(dists)
                elif cur == ans:
                    if max(dists) < fx:
                        ans,ansDataIndex = cur,i
                        fx = max(dists)
                vps = vps[:-1]
            vps.append(datas[ansDataIndex])

    utils.writeDataToFile('vp/vp_%d_%d_%d_%s.txt'%(dim,numberOfVP,cardinality,typeOfVP),vps)




#utils.writeDataToFile('vp/vp_%d_%d_%d_%s.txt'%(dim,numberOfVP,cardinality,typeOfVP),vps)

if __name__ == '__main__':
    options = utils.getOptions()

    utils.createDirectory('vp')
#generateVantagePointsWithManyAlgorithm(options)
    generateRandomVP(options)
