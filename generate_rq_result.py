#!/usr/bin/python
#-*- coding:utf-8 -*-
import utils

if __name__ == '__main__':
    utils.createDirectory('rq_result')
    options = utils.getOptions()
    queryRange = options['queryRange']

    dataFileName    = utils.getDataFileName(options)
    queryFileName   = utils.getQueryFileName(options)
    rqResultFileName= utils.getRQResultFileName(options)
    datas           = utils.getDataInFile(dataFileName)
    querys          = utils.readDataFromFile(queryFileName)

    with open(rqResultFileName,'w') as fp:
        for i in xrange(len(querys)):
            print '#%d'%i
            fp.write(('#%d'%(i))+'\n')
            ans = []
            for j in xrange(len(datas)):
                dist = utils.hammingDistance(querys[i],datas[j])
                if dist <= queryRange:
                    ans.append(datas[j])
            for data in ans:
                for j in xrange(len(data)):
                    fp.write('%c '%data[j])
                fp.write('\n')
