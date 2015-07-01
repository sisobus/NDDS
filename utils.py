#/usr/bin/python
#-*- coding:utf-8 -*-
import sys
import math
import string
import os
import glob
import commands
from optparse import OptionParser
import matplotlib.pyplot as plt

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
    parser.add_option('-v','--vptype',default='greedy',
            help='type of vantage point [greedy,random,AA,AB]',dest='vptype')
    (options, args) = parser.parse_args(sys.argv[1:])

    numberOfData        = int(options.numberOfData)
    numberOfVP          = int(options.numberOfVP)
    numberOfDimension   = int(options.numberOfDimension)
    numberOfAlphabet    = int(options.numberOfAlphabet)
    distribution        = str(options.distribution)
    typeOfVP            = str(options.vptype)
    ret = {
        'numberOfData'      : numberOfData,
        'numberOfVP'        : numberOfVP,
        'numberOfDimension' : numberOfDimension,
        'numberOfAlphabet'  : numberOfAlphabet,
        'distribution'      : distribution,
        'typeOfVP'          : typeOfVP,
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

def saveGraphWithHighValue(imageFileName,xp,yp,highValue):
    plt.plot(xp,yp,lw=2)
    plt.text(len(xp),0,str(highValue))
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
        datas.append(line.rstrip().split(' '))
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

def getDataFileName(options):
    size            = options['numberOfData']
    dim             = options['numberOfDimension']
    distribution    = options['distribution']
    cardinality     = options['numberOfAlphabet']
    dataFileName    = 'data/data_%d_%d_%s_%d.txt'%(size,dim,distribution,cardinality)
    return dataFileName

def getQueryFileName(options):
    size            = options['numberOfData']
    dim             = options['numberOfDimension']
    distribution    = options['distribution']
    cardinality     = options['numberOfAlphabet']
    queryFileName   = 'query/query_%d_%d_%s_%d.txt'%(size,dim,distribution,cardinality)
    return queryFileName

def getVPFileName(options):
    dim             = options['numberOfDimension']
    numberOfVP      = options['numberOfVP']
    cardinality     = options['numberOfAlphabet']
    typeOfVP        = options['typeOfVP']
    vpFileName      = 'vp/vp_%d_%d_%d_%s.txt'%(dim,numberOfVP,cardinality,typeOfVP)
    return vpFileName

def getCDSDataFileName(options):
    size            = options['numberOfData']
    numberOfVP      = options['numberOfVP']
    distribution    = options['distribution']
    cardinality     = options['numberOfAlphabet']
    typeOfVP        = options['typeOfVP']
    cdsDataFileName = 'cds_data/data_%d_%d_%s_%d_%s.txt'%(size,numberOfVP,distribution,cardinality,typeOfVP)
    return cdsDataFileName

def getCDSQueryFileName(options):
    size                = options['numberOfData']
    distribution        = options['distribution']
    cardinality         = options['numberOfAlphabet']
    numberOfVP          = options['numberOfVP']
    typeOfVP            = options['typeOfVP']
    cdsQueryFileName    = 'cds_query/query_%d_%d_%s_%d_%s.txt'%(size,numberOfVP,distribution,cardinality,typeOfVP)
    return cdsQueryFileName

def getImageFileName(options,tag):
    size            = options['numberOfData']
    dim             = options['numberOfDimension']
    distribution    = options['distribution']
    cardinality     = options['numberOfAlphabet']
    numberOfVP      = options['numberOfVP']
    typeOfVP        = options['typeOfVP']
    imageFileName   = 'figure/figure_%d_%d_%s_%d_%d_%s_%d.png'%(size,dim,distribution,cardinality,numberOfVP,typeOfVP,tag)
    return imageFileName

