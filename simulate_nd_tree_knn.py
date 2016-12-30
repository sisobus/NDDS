#!/usr/bin/python
#-*- coding:utf-8 -*-
import utils
import os
import commands
import glob

def replaceIntInString(s,x):
    ret = s.split('=')[0]
    ret = ret + '= %d;'%(x)
    return ret

def setConfigHeaderFile(options):
    size = options['numberOfData']
    dim  = options['numberOfDimension']
    numberOfAlphabet = options['numberOfAlphabet']

    originalConfigFileName = 'nd_tree_knn/original_config.h'
    resultConfigFileName = 'nd_tree_knn/config.h'

    with open(originalConfigFileName,'r') as fp:
        lines = fp.read().rstrip().split('\n')
    with open(resultConfigFileName,'w') as fp:
        for line in lines:
            line = line.rstrip()
            if line.find('MUST_CHANGE_FOR_TEST') == -1:
                fp.write(line+'\n')
                continue
            if line.find('MAX_LINE_IN_SOURCE_FILE') <> -1:
                replaceString = replaceIntInString(line,size)
                fp.write(replaceString+'\n')
                continue
            if line.find('DIM') <> -1:
                replaceString = replaceIntInString(line,dim)
                fp.write(replaceString+'\n')
                continue
            if line.find('const int ALPHA') <> -1:
                replaceString = replaceIntInString(line,numberOfAlphabet)
                print replaceString
                fp.write(replaceString+'\n')
                continue
    print 'setConfigHeaderFile complete'

def copyDataFile(options):
    dataFileName    = utils.getNDTDataFileName(options)
    dim             = options['numberOfDimension']

    if not os.path.exists(dataFileName):
        print 'copyDataFile error:%s is not exists'%dataFileName
        exit(-1)
    command = 'cp %s nd_tree_knn/sourceData%d+0.txt'%(dataFileName,dim)
    utils.executeCommand(command)
    print 'copyDataFile complete'

def copyQueryFile(options):
    queryFileName   = utils.getNDTQueryFileName(options)

    if not os.path.exists(queryFileName):
        print 'copyQueryFile error:%s is not exists'%queryFileName
        exit(-1)
    command = 'cp %s nd_tree_knn/rangequeryAll.txt'%queryFileName
    utils.executeCommand(command)
    command = 'cp %s nd_tree_knn/boxqueryAll.txt'%queryFileName
    utils.executeCommand(command)
    print 'copyQueryFile complete'

def executeMake():
    command = 'make -C nd_tree_knn'
    utils.executeCommand(command)
    print 'make complete'

def executeNDT():
    command = 'cd nd_tree_knn ; ./NDT'
    result = utils.executeCommand(command)
    print result


if __name__ == '__main__':
    options = utils.getOptions()

    setConfigHeaderFile(options)
    copyDataFile(options)
    copyQueryFile(options)
    executeMake()
    executeNDT()
