#!/usr/bin/python
#-*- coding:utf-8 -*-
import utils
import os
import commands
import glob

def copyDataFile(options):
    dataFileName    = utils.getCDSDataFileName(options)
    dim             = options['numberOfDimension']

    if not os.path.exists(dataFileName):
        print 'copyDataFile error:%s is not exists'%dataFileName
        exit(-1)
    command = 'cp %s kd_tree/data.txt'
    utils.executeCommand(command)
    print 'copyDataFile complete'

def copyQueryFile(options):
    queryFileName   = utils.getCDSQueryFileName(options)

    if not os.path.exists(queryFileName):
        print 'copyQueryFile error:%s is not exists'%queryFileName
        exit(-1)
    command = 'cp %s kd_tree/query.txt'
    utils.executeCommand(command)
    print 'copyQueryFile complete'

def executeMake():
    command = 'make -C kd_tree'
    utils.executeCommand(command)
    print 'make complete'

def executeNDT(options):
    dim = options['numberOfDimension']
    size = options['numberOfData']
    queryRange = options['queryRange']
    numberOfVP = options['numberOfVP']
    command = 'cd kd_tree ; ./kdtree -load_file data.txt -orig_dim %d -dim %d -rqfile query.txt -range %d -count %d'%(dim,numberOfVP,queryRange,size)
    result = utils.executeCommand(command)
    print result


if __name__ == '__main__':
    options = utils.getOptions()

    copyDataFile(options)
    copyQueryFile(options)
    executeMake()
    executeNDT(options)
