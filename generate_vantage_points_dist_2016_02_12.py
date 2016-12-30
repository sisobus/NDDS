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

def popDataAtIndex(vps,n):
    return vps[:n]+vps[(n+1):]

def getWorstVP(vps):
    cnt,fx,idx = -1,987654321,-1
    for i in xrange(len(vps)):
        curVPS = popDataAtIndex(vps,i)
        curCnt,curFx = calculateMany(curVPS)
        if cnt == curCnt:
            if fx > curFx:
                cnt,fx,idx = curCnt,curFx,i
        elif cnt < curCnt:
            cnt,fx,idx = curCnt,curFx,i

    return idx

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
    isSelected = [ False for i in xrange(len(datas)) ]
    for i in xrange(len(datas)):
        if datas[i] == majorPattern:
            isSelected[i] = True
    notChangedCount = 0
    while len(vps) < numberOfVP:
        print len(vps)
        changed = False
        for i in xrange(len(datas)):
            if isSelected[i]: 
                continue
            is_pass = False
            for j in xrange(len(vps)):
                if datas[i] == vps[j]:
                    is_pass = True
                dist = utils.hammingDistance(datas[i],vps[j])
                if d[dist] > threshold:
                    is_pass = True
            if is_pass:
                continue
            for j in xrange(len(vps)):
                dist = utils.hammingDistance(datas[i],vps[j])
                d[dist] += 1
            vps.append(datas[i])
            isSelected[i] = True
            changed = True
        if not changed:
            print 'not changed so pop worst (%d)'%notChangedCount
            worstIdx = getWorstVP(vps)
            for j in xrange(len(vps)):
                if j == worstIdx: continue
                dist = utils.hammingDistance(vps[j],vps[worstIdx])
                d[dist] -= 1
            nextVPS = popDataAtIndex(vps,worstIdx)
            vps = nextVPS
            notChangedCount += 1
            if notChangedCount > numberOfVP/2:
                notChangedCount = 0
            for k in xrange(notChangedCount-1):
                worstIdx = getWorstVP(vps)
                for j in xrange(len(vps)):
                    if j == worstIdx: continue
                    dist = utils.hammingDistance(vps[j],vps[worstIdx])
                    d[dist] -= 1
                nextVPS = popDataAtIndex(vps,worstIdx)
                vps = nextVPS

        print len(vps)
        print d


#    utils.writeDataToFile('vp/vp_%d_%d_%d_%s.txt'%(dim,numberOfVP,cardinality,typeOfVP),vps[:numberOfVP])

if __name__ == '__main__':
    options = utils.getOptions()

    utils.createDirectory('vp')
    generateVantagePointsWithManyAlgorithm(options)
