import os
import csv
import sys
import numpy

def main(fin,fout):

    print 'importing data'

    ids = (['chrID','src','ftype','lb','ub','dot1','str','dot2','info',
        'af1e-5','af1e-4','af1e-3','af1e-2','af1e-1','af1',
        'het1e-5','het1e-4','het1e-3','het1e-2','het1e-1','het1',
        'ac','an','ptype'])
    #Import the sample data and concatenate
    offset = 9
    offset_het = 15
    out = [numpy.zeros((6,2)),0,0]
    out[0] = out[0].astype(float)

    ct = 0
    with open(fin, 'rb') as csvin:
        csvreader = csv.reader(csvin, delimiter=',')
        for row in csvreader:
            ct+=1
            if (ct % 1000) == 0:
                print out
                print ct 

            out[1]+=int(row[21])
            out[2]+=int(row[22])
            if int(row[21]) > 0:                  
                for i in range(0,6):
                    el = row[i+offset].replace('[','').replace(']','').replace(',','').strip().split()
                    el_het = row[i+offset_het].replace('[','').replace(']','').replace(',','').strip().split()

                    # print el
                    # print el_het

                    for j in range(len(el)):
                        out[0][i,0]+=float(el[j])
                        out[0][i,1]+=float(el_het[j])
    
    print out

    with open(fout,'wb') as csvout:
        mywriter = csv.writer(csvout, delimiter=',')           

        for i in range(0,6):
            # act_sum = out[0][i,0]
            # het_sum = out[0][i,1]
            # hom_ratio = (float(1)-(float(het_sum)/float(act_sum)))
            mywriter.writerow(list(out[0][i,:])+[out[1],out[2]]+[ids[i+offset]])                                     

fin = sys.argv[1]
fout = sys.argv[2]

main(fin,fout)


