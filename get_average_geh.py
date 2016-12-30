#!/usr/bin/python
#-*- coding:utf-8 -*-
# python convert_ndds_to_cds.py -n 100000 -d 10 -m 10 -a 4 -b u -v greedy
import os
import utils

def convertNDDSToCDS(options):
    size            = options['numberOfData']
    dim             = options['numberOfDimension']
    distribution    = options['distribution']
    cardinality     = options['numberOfAlphabet']
    numberOfVP      = options['numberOfVP']
    typeOfVP        = options['typeOfVP']

    dataFileName    = 'data/data_%d_%d_%s_%d.txt'%(size,dim,distribution,cardinality)
    queryFileName   = 'query/query_%d_%d_%s_%d.txt'%(size,dim,distribution,cardinality)
    vpFileName      = 'vp/vp_%d_%d_%d_%s.txt'%(dim,numberOfVP,cardinality,typeOfVP)
#    cdsDataFileName = utils.getCDSDataFileName(options)
#    cdsQueryFileName= utils.getCDSQueryFileName(options)

    datas   = utils.getDataInFile(dataFileName)
    querys  = utils.readDataFromFile(queryFileName)
    vps     = utils.readDataFromFile(vpFileName)
    print len(datas),len(querys),len(vps)

    d = [ {} for i in xrange(dim) ]
    for i in xrange(len(datas)):
        for j in xrange(dim):
            if datas[i][j] in d[j]:
                d[j][datas[i][j]] += 1
            else:
                d[j][datas[i][j]] = 1

    for i in xrange(dim):
        for key in d[i]:
            d[i][key] = 1.0-float(d[i][key])/float(len(datas))
            
    def geh(a,b):
        ret = 0.0
        for i in xrange(len(a)):
            if a[i] != b[i]:
                ret += 1.0
            else:
                ret += d[i][a[i]]/float(dim)
        return ret/float(dim)
    datas = datas[:1000]
    average = 0.0
    for i in xrange(len(datas)):
        cur = 98765432
        for j in xrange(len(datas)):
            if i != j:
                cur = min(cur,geh(datas[i],datas[j]))
        average += cur
    average /= len(datas)
    print average



if __name__ == '__main__':
    utils.createDirectory('cds_data')
    utils.createDirectory('cds_query')

    options = utils.getOptions()
    convertNDDSToCDS(options)
