#!/usr/bin/python
#-*- coding:utf-8 -*-
import utils
import os
import commands
import glob

def copyDataFile(options):
    dataFileName    = utils.getNDTDataFileName(options)
    print dataFileName
    dim             = options['numberOfDimension']

    if not os.path.exists(dataFileName):
        print 'copyDataFile error:%s is not exists'%dataFileName
        exit(-1)
    command = 'cp %s mvptree/data.txt'%dataFileName
    utils.executeCommand(command)
    print 'copyDataFile complete'

def copyQueryFile(options):
    queryFileName   = utils.getNDTQueryFileName(options)
    print queryFileName

    if not os.path.exists(queryFileName):
        print 'copyQueryFile error:%s is not exists'%queryFileName
        exit(-1)
    command = 'cp %s mvptree/query.txt'%queryFileName
    utils.executeCommand(command)
    print 'copyQueryFile complete'

def executeMake():
    command = 'make -C mvptree'
    utils.executeCommand(command)
    print 'make complete'

def executeNDT(options):
    dim = options['numberOfDimension']
    size = options['numberOfData']
    numberOfVP = options['numberOfVP']
    for i in xrange(7):
        queryRange = i+1
        command = 'cd mvptree ; ./mvptree -load_file data.txt -dim %d -ref_method mvp -rqfile query.txt -range %d -count %d'%(dim,queryRange,size)
        print command
        result = utils.executeCommand(command)
        print result


if __name__ == '__main__':
    options = utils.getOptions()

    copyDataFile(options)
    copyQueryFile(options)
    executeMake()
    executeNDT(options)
