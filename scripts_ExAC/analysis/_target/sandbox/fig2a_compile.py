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
    #26 is longer than most guides
    #equivalent to including extra N's
    #makes lengths of output for all prots equivalent
    out = [numpy.zeros((6,26)),0]

    ct = 0
    with open(fin, 'rb') as csvin:
        csvreader = csv.reader(csvin, delimiter=',')
        for row in csvreader:
            ct+=1
            if (ct % 1000) == 0:
                print ct 

            # print ''
            # print row

            temp_sum = 0
            out[1]+=int(row[22]) 
            if int(row[21]) > 0:                  
                for i in range(0,6):
                    el = row[i+offset].replace('[','').replace(']','').replace(',','').strip().split()
                    el2 = numpy.zeros(26)

                    for j in range(len(el)):
                        if int(el[j]) > 0:
                            # el2+=[float(x)/float(row[22])]
                            el2[j] = int(el[j])
                            temp_sum+=int(el[j])

                    out[0][i,:]+=el2
    
                if temp_sum > int(row[22]):
                    print 'error: sum is greater than whole'
                    print temp_sum
                    print row[22]
                    print 'exiting'
                    print exit(1)

    #PAM currenly on left
    #fliplr to orient PAM on right 
    # out[0] = numpy.fliplr(out[0])
    
    with open(fout,'wb') as csvout:
        mywriter = csv.writer(csvout, delimiter=',')
        var_total = numpy.matrix(out[0]).sum()
        if out[1] > 0:
            out[0] = out[0].astype(float)/float(out[1])            

        for i in range(0,6):
            mywriter.writerow(list(out[0][i,:])+[var_total,out[1]]+[ids[i+offset]])                                      

fin = sys.argv[1]
fout = sys.argv[2]

main(fin,fout)


