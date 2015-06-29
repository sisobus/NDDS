#/usr/bin/python
#-*- coding:utf-8 -*-
import sys
import math
import string
import os
import glob
import commands
from optparse import OptionParser

def getOptions():
    parser = OptionParser()
    parser.add_option('-n','--size',default='100000',
            help='input number of data',dest='numberOfData')
    parser.add_option('-m','--vpsize',default='25',
            help='input number of vp',dest='numberOfVP')
    parser.add_option('-d','--dimension',default='25',
            help='input number of dimension',dest='numberOfDimension')
    parser.add_option('-a','--cardinality',default='4',
            help='input number of alphabet',dest='numberOfAlphabet')
    parser.add_option('-b','--distribution',default='u',
            help='uniform: u skewed: c1 10-cluster: c10',dest='distribution')
    (options, args) = parser.parse_args(sys.argv[1:])

    numberOfData        = int(options.numberOfData)
    numberOfVP          = int(options.numberOfVP)
    numberOfDimension   = int(options.numberOfDimension)
    numberOfAlphabet    = int(options.numberOfAlphabet)
    distribution        = str(options.distribution)
    ret = {
        'numberOfData'      : numberOfData,
        'numberOfVP'        : numberOfVP,
        'numberOfDimension' : numberOfDimension,
        'numberOfAlphabet'  : numberOfAlphabet,
        'distribution'      : distribution,
            }
    return ret

def createDirectory(directoryName):
    if not os.path.exists(directoryName):
        command = 'mkdir %s'%directoryName
        ret = commands.getoutput(command)

def saveGraph(imageFileName,xp,yp):
    plt.plot(xp,yp,lw=2)
    plt.xlim(0,len(xp)+1)
    plt.savefig(imageFileName,dpi=100)
    plt.clf()

def getDataInFile(filename):
    onlyFileName = filename.split('.')[0]
    numberOfDataSize  = int(onlyFileName.split('_')[1])
    numberOfDimension = int(onlyFileName.split('_')[2])
    numberOfAlphabet  = int(onlyFileName.split('_')[4])

    ret = []
    with open(filename,'r') as fp:
        r = fp.read().rstrip().split('\n')
    for i in xrange(1,numberOfDataSize+1):
        t = r[i].split(':')[1].split(',')
        ret.append(t)
    return ret

def hammingDistance(a,b):
    n = len(a)
    if n <> len(b):
        print 'calculate hamming distance error:length <>'
        exit(-1)
    ret = 0
    for i in xrange(n):
        if a[i] <> b[i]:
            ret = ret + 1
    return ret

def writeDataToFile(filename,datas):
    if os.path.exists(filename):
        print 'utils.writeDataToFile error: %s is exists'%filename
        return
    with open(filename,'w') as fp:
        for i in xrange(len(datas)):
            for j in xrange(len(datas[i])):
                fp.write(str(datas[i][j])+' ')
            fp.write('\n')

def readDataFromFile(filename):
    datas = []
    with open(filename,'r') as fp:
        lines = fp.read().rstrip().split('\n')
    for line in lines:
        datas.append(line.aplit(' '))
    return datas

def calculateCorrelationCoefficient(vp1,vp2,datas):
    sumOfSubD1iAndD2i = 0
    sumOfD1i          = 0
    sumOfD2i          = 0
    sumOfD1iSquare    = 0
    sumOfD2iSquare    = 0
    
    size = len(datas)
    dim  = len(vp1)
    for i in xrange(size):
        d1i = hammingDistance(vp1,datas[i])
        d2i = hammingDistance(vp2,datas[i])

        sumOfSubD1iAndD2i = sumOfSubD1iAndD2i + (d1i*d2i)
        sumOfD1i          = sumOfD1i + d1i
        sumOfD2i          = sumOfD2i + d2i
        sumOfD1iSquare    = sumOfD1iSquare + (d1i*d1i)
        sumOfD2iSquare    = sumOfD2iSquare + (d2i*d2i)

    ret = (size*sumOfSubD1iAndD2i - sumOfD1i*sumOfD2i) / (math.sqrt(size*sumOfD1iSquare - (sumOfD1i*sumOfD1i)) * math.sqrt(size*sumOfD2iSquare - (sumOfD2i*sumOfD2i)) )
    return ret
