import re
import csv
import sys   
import time 
import difflib    
import urllib2
import StringIO
import xml.etree.ElementTree
from operator import itemgetter 
from xml.dom.minidom import parse, parseString    
from Bio import SeqIO
from Bio.Seq import Seq    
from Bio import pairwise2
from Bio.SeqRecord import SeqRecord

import numpy
import pysam
from scipy.stats import mode
from scipy import stats

hg19_offset = {1:[0,[0,1]],
2:[249250621,[1,0]],
3:[492449994,[0,1]],
4:[690472424,[1,0]],
5:[881626700,[0,1]],
6:[1062541960,[1,0]],
7:[1233657027,[0,1]],
8:[1392795690,[1,0]],
9:[1539159712,[0,1]],
10:[1680373143,[1,0]],
11:[1815907890,[0,1]],
12:[1950914406,[1,0]],
13:[2084766301,[0,1]],
14:[2199936179,[1,0]],
15:[2307285719,[0,1]],
16:[2409817111,[1,0]],
17:[2500171864,[0,1]],
18:[2581367074,[1,0]],
19:[2659444322,[0,1]],
20:[2718573305,[1,0]],
21:[2781598825,[0,1]],
22:[2829728720,[1,0]],
'X':[2881033286,[0,1]],
'Y':[3036303846,[1,0]],
'M':[3095677412,[0,1]]}

def compile_manhat(rawdata_file,out_file):

    rawdata = [] 
    try:
        with open(rawdata_file, 'rb') as csvfile:
            print 'opened csv'
            csvreader = csv.reader(csvfile, delimiter=',')
            rawdata = [row for row in csvreader]
    except:
        print "Could not open file " + rawdata_file
        sys.exit(1)         
  
    rawdata_n = map(itemgetter(1), rawdata)
    rawdata_n = map(int, rawdata_n)
    rawdata_loc = map(itemgetter(12), rawdata)
    rawdata_loc = map(int, rawdata_loc)
    rawdata_names = map(itemgetter(13), rawdata)
    rawdata_ct = [len(el.split(',')) for el in rawdata_names]
    rawdata_info = map(itemgetter(6), rawdata)
    rawdata_chrid = [el.split('_')[1] for el in rawdata_info]
    rawdata_offset = [el.split('_')[2].split('.')[0] for el in rawdata_info]
    rawdata_offset = map(int, rawdata_offset)

    chrid = 0
    chrmed = 0

    with open(out_file,'w') as csvfile:
        mywriter = csv.writer(csvfile, delimiter=',')
        for i in range(0,len(rawdata_chrid)):
            chrid = rawdata_chrid[i]
            chrmed = rawdata_loc[i]+rawdata_offset[i]+hg19_offset[chrid][0]

            print 'rawdata_names[i]: '+rawdata_names[i]
            print 'rawdata_info[i]: '+rawdata_info[i]
            print 'chrid: '+str(chrid)
            print 'chrmed: '+str(chrid)
            print 'rawdata_ct[i]: '+str(rawdata_ct[i])

            mywriter.writerow([chrmed,
                float(hg19_offset[chrid][1][0]),
                float(hg19_offset[chrid][1][1])])

    return True


def main():

    fin = sys.argv[1]
    fout = sys.argv[2]
    # todo_done = sys.argv[3]

    success = compile_manhat(fin,fout)         

    # line = fin+' '+fout
    # if success:
    #     with open(todo_done,'a') as f:
    #         f.write(line+'\n')

main()


