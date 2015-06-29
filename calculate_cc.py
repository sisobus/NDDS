#!/usr/bin/python
#-*- coding:utf-8 -*-
import utils
import os
import glob

class CC:
    def __init__(self):
        self.x = 0
        self.id1 = 0
        self.id2 = 0
    def __init__(self,x,id1,id2):
        self.x = x
        self.id1 = id1
        self.id2 = id2
    def __lt__(self,other):
        if self.x == other.x:
            return self.id1 < self.id2
        return self.x < other.x

if __name__ == '__main__':
    options = utils.getOptions()
    dataFileName = 'data/data_%d_%d_%s_%d.txt'%(options['numberOfData'],options['numberOfDimension'],options['distribution'],options['numberOfAlphabet'])
    
    print dataFileName
    datas = utils.getDataInFile(dataFileName)

    filenames = glob.glob('vp/*.txt')
    for filename in filenames:
        if int(filename.split('_')[1]) <> options['numberOfDimension']:
            continue
        print filename
        with open(filename,'r') as fp:
            lines = fp.read().rstrip().split('\n')
        vps = []
        for line in lines:
            vps.append(line.rstrip().split(' '))
        cc = []
        for i in xrange(len(vps)):
            for j in xrange(len(vps)):
                if i == j:
                    continue
                cc.append(CC(abs(utils.calculateCorrelationCoefficient(vps[i],vps[j],datas)),i,j))
        cc.sort()
        for i in xrange(4):
            cur = -(i+1)
            print cc[cur].x, vps[cc[cur].id1], vps[cc[cur].id2]
        print ''
        for i in xrange(4):
            cur = i
            print cc[cur].x, vps[cc[cur].id1], vps[cc[cur].id2]
        print ''
