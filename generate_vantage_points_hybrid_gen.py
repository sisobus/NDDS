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
    d = [ 0 for i in xrange(len(vps[0])+1) ]
    for i in xrange(len(vps)):
        for j in xrange(len(vps)):
            if i != j :
                now = utils.hammingDistance(vps[i],vps[j])
                d[now] += 1
                s.add(now)
    return len(s),d

def getDifferentCharacter(c,numberOfAlphabet):
    alphabet = list(string.ascii_uppercase)
    while True:
        cur = random.randrange(0,numberOfAlphabet)
        if alphabet[cur] != c:
            return alphabet[cur]
    return None

def changeStringAtIndex(s,i,c):
    return s[:i]+list(c)+s[i+1:]

def generateVpWithDist_backup(numberOfDimension,numberOfAlphabet,vp,distance):
    ret = ''
    for i in xrange(distance):
        ret = ret + getDifferentCharacter(vp[i],numberOfAlphabet)
    for i in xrange(distance,numberOfDimension):
        ret = ret + vp[i]
    return ret


def generateVpWithDist(numberOfDimension,numberOfAlphabet,vp,distance):
    alphabet = list(string.ascii_uppercase)
    ret = vp
    s = set()
    for i in xrange(distance):
        curIndex = random.randrange(0,numberOfDimension)
        while curIndex in s:
            curIndex = random.randrange(0,numberOfDimension)
        s.add(curIndex)
        ret = changeStringAtIndex(ret,curIndex,getDifferentCharacter(vp[curIndex],numberOfAlphabet))
    return ret

def generateVantagePointsWithHybridAlgorithm(options):
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
    one_pass = False 
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
            change = False
            for i in xrange(dim+1):
                if d[i] == 0:
                    change = True
                    ans,ans_vp = -1, ''
                    fx = 987654321
                    for j in xrange(len(vps)):
                        cur_vp = generateVpWithDist(dim,cardinality,vps[j],i)
                        vps.append(cur_vp)
                        cur,dists = calculateMany(vps)
                        if cur > ans:
                            ans,ans_vp = cur,cur_vp
                            fx = max(dists)
                        elif cur == ans:
                            if max(dists) < fx:
                                ans,ans_vp = cur,cur_vp
                                fx = max(dists)
                        vps = vps[:-1]
                    for j in xrange(len(vps)):
                        dist = utils.hammingDistance(vps[j],ans_vp)
                        d[dist] += 1
                    vps.append(ans_vp)
                    break
            if not change:
                vps.append(datas[random.randrange(0,numberOfData)])


    utils.writeDataToFile('vp/vp_%d_%d_%d_%s.txt'%(dim,numberOfVP,cardinality,typeOfVP),vps)




#utils.writeDataToFile('vp/vp_%d_%d_%d_%s.txt'%(dim,numberOfVP,cardinality,typeOfVP),vps)

if __name__ == '__main__':
    options = utils.getOptions()

    utils.createDirectory('vp')
    generateVantagePointsWithHybridAlgorithm(options)
