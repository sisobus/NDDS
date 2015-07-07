#!/usr/bin/python
#-*- coding:utf-8 -*-
import utils
import os
import string
import numpy as np

if __name__ == '__main__':
    utils.createDirectory('figure_pair')
    options = utils.getOptions()
    dim     = options['numberOfDimension']

    dataFileName    = utils.getDataFileName(options)
    vpFileName      = utils.getVPFileName(options)
    datas           = utils.getDataInFile(dataFileName)
    vps             = utils.readDataFromFile(vpFileName)

    for i in xrange(len(vps)):
        for j in xrange(i+1,len(vps)):
            imageFileName = utils.getFigurePairName(options,i,j)
            print imageFileName
            if os.path.exists(imageFileName):
                print '%s is exists'%imageFileName
                continue
            xp = []
            yp = []
            zp = []
            zcnt = [ [ 0 for ii in xrange(dim+1) ] for jj in xrange(dim+1) ]
            cc = utils.calculateCorrelationCoefficient(vps[i],vps[j],datas)
            for k in xrange(len(datas)):
                x = utils.hammingDistance(vps[i],datas[k])
                y = utils.hammingDistance(vps[j],datas[k])
                zcnt[y][x] = zcnt[y][x] + 1
                xp.append(x)
                yp.append(y)
            for k in xrange(len(xp)):
                zp.append(zcnt[yp[k]][xp[k]])
            #utils.saveGraphUsingPointWithCC(imageFileName,xp,yp,cc,dim)
            utils.saveGraphUsing3DSurfaceWithCC(imageFileName,xp,yp,zp,cc,dim)


