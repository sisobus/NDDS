#!/usr/bin/python
#-*- coding:utf-8 -*-
import utils
import os
import string

if __name__ == '__main__':
    utils.createDirectory('figure')
    options = utils.getOptions()

    dataFileName    = utils.getDataFileName(options)
    vpFileName      = utils.getVPFileName(options)
    datas           = utils.getDataInFile(dataFileName)
    vps             = utils.readDataFromFile(vpFileName)

    curDatas = datas
    for i in xrange(len(vps)):
        print i
        n = len(curDatas)
        x = [ [] for j in xrange(len(vps[i])+1) ]
        for j in xrange(n):
            nextPosition = utils.hammingDistance(vps[i],curDatas[j])
            x[nextPosition].append(j)
        mx, position = (0,0)
        xp = []
        yp = []
        for j in xrange(len(vps[i])+1):
            xp.append(j)
            yp.append(len(x[j]))
            if mx < len(x[j]):
                mx, position = (len(x[j]),j)
        imageFileName = utils.getImageFileName(options,i)
        if os.path.exists(imageFileName):
            print '%s is exists'%(imageFileName)
        elif not os.path.exists(imageFileName):
            utils.saveGraphWithHighValue(imageFileName,xp,yp,mx)
        
        curDatas = []
        for dataId in x[position]:
            curDatas.append(datas[dataId])

