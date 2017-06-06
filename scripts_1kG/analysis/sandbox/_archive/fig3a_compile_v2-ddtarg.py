import os
import csv
import sys
import numpy

def main(fin,fout):

    print 'importing data'

    ids = (['transcript','targID','PAM','PAMori','str','targ','searchID',
        'sMM','sTarg','sPAM','sPAMori','sstr','sloc','snames',
        'rMM','rPAMflag','rTarg','rPAM','rPAMori','rstr','namect','targct','not_ontarg'])

    ct = 0
    rawdata_n = 0
    rawdata_loc = 0
    rawdata_names = ''
    rawdata_ct = 0
    rawdata_info = ''
    rawdata_chrid = ''
    rawdata_offset = 0
    out = [numpy.zeros((2,6)), 0]
    with open(fin, 'rb') as csvin:
        csvreader = csv.reader(csvin, delimiter=',')
        for row in csvreader:
            if ct == 0:
                feat_id = row[0]
                targct = int(row[len(row)-2])
                namect = int(row[len(row)-3]) 
                out[1]+=targct           

            ct+=1
            if (ct % 1000) == 0:
                print ct 

            rawdata_info = row[6]
            rawdata_loc = int(row[12])
            rawdata_chrid = rawdata_info.split('_')[1]
            rawdata_offset = int(rawdata_info.split('_')[2].split('.')[0])

            #remove on-targets
            not_ontarg = int(row[len(row)-1])

            #increment counts for each level of MM
            if not_ontarg:
                rawdata_n = int(row[7])
                rawdata_names = row[13].replace('[','').replace(']','')
                rawdata_ct = len(rawdata_names.split(','))
                if rawdata_names == '\'0\'':                 
                    out[0][0,rawdata_n]+=1
                else:
                    #only increment array if:
                    #off-target is present from drop in distance of alt allele to targ
                    #compared to distance between ref allele and targ
                    if int(row[14]) > int(row[7]):
                        out[0][1,rawdata_n]+=1
            else:
                print 'on-target: continuing'
    
    #write data structure for input into R
    with open(fout,'wb') as csvout:
        mywriter = csv.writer(csvout, delimiter=',')
        out[0] = out[0].astype(float)/float(out[1])            
        for i in range(0,2):
            mywriter.writerow(list(out[0][i,:])+[feat_id,i])                                      

fin = sys.argv[1]
fout = sys.argv[2]

main(fin,fout)


