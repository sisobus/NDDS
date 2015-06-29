#!/usr/bin/python
#-*- coding:utf-8 -*-
import utils
import os
import glob

class A:
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

datas = utils.getDataInFile('data/data_100000_10_u_4.txt')

filenames = glob.glob('vp/*.txt')
for filename in filenames:
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
            cc.append(A(abs(utils.calculateCorrelationCoefficient(vps[i],vps[j],datas)),i,j))
    cc.sort()
    print cc[-1].x, vps[cc[-1].id1], vps[cc[-1].id2]
    print cc[-2].x, vps[cc[-2].id1], vps[cc[-2].id2]
    print cc[-3].x, vps[cc[-3].id1], vps[cc[-3].id2]
    print cc[-4].x, vps[cc[-4].id1], vps[cc[-4].id2]
    print ''
    print cc[0].x, vps[cc[0].id1], vps[cc[0].id2]
    print cc[1].x, vps[cc[1].id1], vps[cc[1].id2]
    print cc[2].x, vps[cc[2].id1], vps[cc[2].id2]
    print cc[3].x, vps[cc[3].id1], vps[cc[3].id2]
    print ''
