#!/usr/bin/python
#-*- coding:utf-8 -*-
import os
import glob
import sys
import utils
import string
alphabet     = list(string.ascii_uppercase)
def makeDictionaryKeyIsAlphabet():
    ret = {}
    for i in xrange(len(alphabet)):
        ret[alphabet[i]] = i
    return ret

if __name__ == '__main__':
    utils.createDirectory('ndt_data')
    utils.createDirectory('ndt_query')
    dictionary = makeDictionaryKeyIsAlphabet()

    dataFileNames = glob.glob('data/*.txt')
    for dataFileName in dataFileNames:
        print dataFileName
        onlyFileName = dataFileName.split('.')[0].split('/')[1]
        size        = onlyFileName.split('_')[1]
        dim         = onlyFileName.split('_')[2]
        vptype      = onlyFileName.split('_')[3]
        cardinality = onlyFileName.split('_')[4]

        queryFileName = 'query/query_%s_%s_%s_%s.txt'%(size,dim,vptype,cardinality)
        datas = utils.getDataInFile(dataFileName)
        querys = utils.readDataFromFile(queryFileName)

        ndtDataFileName = 'ndt_data/data_%s_%s_%s_%s.txt'%(size,dim,vptype,cardinality)
        ndtQueryFileName = 'ndt_query/query_%s_%s_%s_%s.txt'%(size,dim,vptype,cardinality)

        if os.path.exists(ndtDataFileName):
            print '%s is exists'%ndtDataFileName
        else :
            with open(ndtDataFileName,'w') as fp:
                for i in xrange(len(datas)):
                    for j in xrange(len(datas[i])):
                        fp.write(str(dictionary[datas[i][j]])+' ')
                    fp.write('\n')
        if os.path.exists(ndtQueryFileName):
            print '%s is exsts'%ndtQueryFileName
        else :
            with open(ndtQueryFileName,'w') as fp:
                for i in xrange(len(querys)):
                    for j in xrange(len(querys[i])):
                        fp.write(str(dictionary[datas[i][j]])+' ')
                    fp.write('\n')

