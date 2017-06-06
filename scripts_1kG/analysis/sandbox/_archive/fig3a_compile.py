import os
import csv
import sys
import numpy

def main(fin,fin_annot,namect,targct,fout):

    print 'importing data'

    ids = (['transcript','targID','PAM','PAMori','str','targ','searchID',
        'sMM','sTarg','sPAM','sPAMori','sstr','sloc','snames',
        'rMM','rPAMflag','rTarg','rPAM','rPAMori','rstr','namect','targct','not_ontarg'])
    feat_id = fin_annot.split('_')[2]+'_'+fin_annot.split('_')[3]

    ct = 0
    bounds = []
    with open(fin_annot, 'rb') as csvin:
        csvreader = csv.reader(csvin, delimiter=',')
        for row in csvreader:
            #has header
            if ct == 0:
                ct+=1
                continue
            el3 = int(row[3].replace(',',''))
            el4 = int(row[4].replace(',',''))
            if el3 <= el4:
                bounds+=[[row[2],el3,el4]]
            else:
                bounds+=[[row[2],el4,el3]]

    ct = 0
    rawdata_n = 0
    rawdata_loc = 0
    rawdata_names = ''
    rawdata_ct = 0
    rawdata_info = ''
    rawdata_chrid = ''
    rawdata_offset = 0
    out = [numpy.zeros((2,6)), targct]
    with open(fin, 'rb') as csvin:
        csvreader = csv.reader(csvin, delimiter=',')
        for row in csvreader:
            ct+=1
            if (ct % 1000) == 0:
                print ct 
            
            rawdata_info = row[6]
            rawdata_loc = int(row[12])
            rawdata_chrid = rawdata_info.split('_')[1]
            rawdata_offset = int(rawdata_info.split('_')[2].split('.')[0])

            #remove on-targets
            not_ontarg = 1
            #change from index space to 1-rooted
            chrmed = rawdata_loc+rawdata_offset+1
            for el in bounds:
                if (rawdata_chrid == el[0] and (chrmed >= el[1]) 
                    and (chrmed <= el[2])):
                    print row
                    not_ontarg = 0
                    break

            #increment counts for each level of MM
            if not_ontarg:
                rawdata_n = int(row[7])
                rawdata_names = row[13].replace('[','').replace(']','')
                rawdata_ct = len(rawdata_names.split(','))
                if rawdata_names == '\'0\'':                 
                    out[0][0,rawdata_n]+=1
                else:
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
fin_annot = sys.argv[2]
namect = sys.argv[3]
targct = sys.argv[4]
fout = sys.argv[5]

main(fin,fin_annot,namect,targct,fout)


