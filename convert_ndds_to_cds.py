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
    cdsDataFileName = utils.getCDSDataFileName(options)
    cdsQueryFileName= utils.getCDSQueryFileName(options)
#    cdsDataFileName     = 'cds_data/data_%d_%d_%s_%d_%s.txt'%(size,numberOfVP,distribution,cardinality,typeOfVP)
#    cdsQueryFileName    = 'cds_query/query_%d_%d_%s_%d_%s.txt'%(size,numberOfVP,distribution,cardinality,typeOfVP)

    datas   = utils.getDataInFile(dataFileName)
    querys  = utils.readDataFromFile(queryFileName)
    vps     = utils.readDataFromFile(vpFileName)
    print len(datas),len(querys),len(vps)

    cdsDatas = []
    for i in xrange(len(datas)):
        t = []
        for j in xrange(len(vps)):
            t.append(utils.hammingDistance(datas[i],vps[j]))
        cdsDatas.append(t)
    utils.writeDataToFile(cdsDataFileName,cdsDatas)

    cdsQuerys = []
    for i in xrange(len(querys)):
        t = []
        for j in xrange(len(vps)):
            t.append(utils.hammingDistance(querys[i],vps[j]))
        cdsQuerys.append(t)
    utils.writeDataToFile(cdsQueryFileName,cdsQuerys)
    print cdsDataFileName, cdsQueryFileName


if __name__ == '__main__':
    utils.createDirectory('cds_data')
    utils.createDirectory('cds_query')

    options = utils.getOptions()
    convertNDDSToCDS(options)
