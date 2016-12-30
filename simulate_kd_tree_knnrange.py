#!/usr/bin/python
#-*- coding:utf-8 -*-
import utils
import os
import commands
import glob

def copyDataFile(options):
    dataFileName    = utils.getCDSDataFileName(options)
    print dataFileName
    dim             = options['numberOfDimension']

    if not os.path.exists(dataFileName):
        print 'copyDataFile error:%s is not exists'%dataFileName
        exit(-1)
    command = 'cp %s kd_tree/data.txt'%dataFileName
    utils.executeCommand(command)
    print 'copyDataFile complete'

def copyQueryFile(options):
    queryFileName   = utils.getCDSQueryFileName(options)
    print queryFileName

    if not os.path.exists(queryFileName):
        print 'copyQueryFile error:%s is not exists'%queryFileName
        exit(-1)
    command = 'cp %s kd_tree/query.txt'%queryFileName
    utils.executeCommand(command)
    print 'copyQueryFile complete'

def executeMake():
    command = 'make -C kd_tree'
    utils.executeCommand(command)
    print 'make complete'

def executeNDT(options):
    dim = options['numberOfDimension']
    size = options['numberOfData']
    knn = options['knn']
    numberOfVP = options['numberOfVP']
    ans = 0.0
    d = []

#    rs = [0.212838/2, 0.249343/2, 0.245334/2, 0.281172/2, 0.301700/2, 0.296108/2]
#    rsb = [0.01070,0.00806,0.00560,0.00538,0.00448,0.00398]
    rs = [4.014654/20, 9.830230/20,14.944167/20,22.749706/20,30.443965/20,35.801397/20]
    rsb = [0.474555/10,1.407575/10,2.077860/10,3.521868/10,4.563981/10,5.814988/10]
    i = 0.01
    dy = 0.01
    if int(dim) == 20:
        i = rs[0] 
        dy = rsb[0]
    elif int(dim) == 40:
        i = rs[1]
        dy = rsb[1]
    elif int(dim) == 60:
        i = rs[2]
        dy = rsb[2]
    elif int(dim) == 80:
        i = rs[3]
        dy = rsb[3]
    elif int(dim) == 100:
        i = rs[4]
        dy = rsb[4]
    elif int(dim) == 120:
        i = rs[5]
        dy = rsb[5]
    while i <= 10.0:
        command = 'cd kd_tree ; ./kdtree -load_file data.txt -orig_dim %d -dim %d -rqfile query.txt -range %f -count %d'%(dim,numberOfVP,i,size)
        result = utils.executeCommand(command)
        print result
        hits = float(result.split(',')[7])
        dist = float(result.split(',')[11])
        ans += dist
        if hits >= 1.0 and len(d) == 0:
            d.append(ans)
        if hits >= 5.0 and len(d) == 1:
            d.append(ans)
        if hits >= 10.0 and len(d) == 2:
            d.append(ans)

        if hits >= knn:
            break
        i += dy
    print ans
    print d

#AVG_VISITED_NODE449.549988
#        REF_PT_TIME,3e-05,INSRT_TIME,0.115137,AVG_QUERY_TIME,5.473e-05,AVG_HITS,0.68,AVG_COMP,2744.7,AVG_DIST,1,END


if __name__ == '__main__':
    options = utils.getOptions()

    copyDataFile(options)
    copyQueryFile(options)
    executeMake()
    executeNDT(options)
