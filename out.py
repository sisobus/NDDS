#!/usr/bin/python
#-*- coding:utf-8 -*-
import utils
import os
import glob
# nd_100000_20_u.out
def get_out_filenames(options):
    ret = []
    b = options['distribution']
    for i in xrange(20,121,20):
        s = 'result/nd_100000_%d_%s.out'%(i,b)
        ret.append(s)
    return ret

if __name__ == '__main__':
    options = utils.getOptions()

    filenames = get_out_filenames(options)
    heuristic = []
    rq_io = []
    number_of_nodes = []
    for filename in filenames:
        with open(filename,'r') as fp:
            lines = fp.read().rstrip().split('\n')
        bc = []
        rq = []
        for line in lines:
            if line.find('BestChild') != -1:
                bc.append(line.split(':')[-1].lstrip().rstrip())
            if line.find('I/O') != -1:
                rq.append(line.split(':')[-1].lstrip().rstrip())
            if line.find('Number of disk blocks') != -1:
                number_of_nodes.append(line.split(':')[-1].lstrip().rstrip())
        heuristic.append(bc)
        rq_io.append(rq)

    print 'range query(=2) I/O'
    d = [ i for i in xrange(20,121,20) ]
    for i in xrange(len(rq_io)):
        print 'dim',d[i],':',rq_io[i][1]
    print 'Bestchild heuristic [ covered_area,notcovered_overlap_enlarge,notcovered_area_enlarge,notcovered_area ]'
    for i in xrange(len(heuristic)):
        print 'dim',d[i],':',
        for j in xrange(len(heuristic[i])):
            print heuristic[i][j],
        print ''
