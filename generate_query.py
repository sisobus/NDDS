#!/usr/bin/python
#-*- coding:utf-8 -*-
import random
import os
import sys
import glob
import utils

if __name__ == '__main__':
    utils.createDirectory('query')
    dataFileNames = glob.glob('data/*.txt')
    for dataFileName in dataFileNames:
        onlyFileName = dataFileName.split('.')[0].split('/')[1]
        size        = onlyFileName.split('_')[1]
        dim         = onlyFileName.split('_')[2]
        vptype      = onlyFileName.split('_')[3]
        cardinality = onlyFileName.split('_')[4]
        
        queryFileName = 'query/query_%s_%s_%s_%s.txt'%(size,dim,vptype,cardinality)
        print queryFileName
        datas = utils.getDataInFile(dataFileName)
        queryDatas = []
        for i in xrange(100):
            queryDatas.append(datas[random.randrange(0,int(size))])
        utils.writeDataToFile(queryFileName,queryDatas)
