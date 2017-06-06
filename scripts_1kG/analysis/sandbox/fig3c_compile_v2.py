import os
import csv
import sys
import numpy

def main(fin,fout):

    print 'importing data'

    abin = (['1e-1','1e-2','1e-3','1e-4'])

    #fields in scr-hom rows
    ids = (['tid','tchr','tloc','tPAM','tPAMori','tstr','ttarg',
        'ochr','oMM','otarg','oPAM','oPAMori','ostr','oloc',
        'ac','ac_het','ac_hom','an','names',
        'rMM','rPAMflag','rtarg','rPAM','rPAMori','rstr',
        'targct','not_ontarg','src'])

    id2ind = {}
    for i in range(len(ids)):
        id2ind[ids[i]] = i

    fin_base = os.path.basename(fin)
    prot_id = fin_base.split('.')[1]

    ct = 0
    alt_ct = 0
    ref_ct = 0
    nmm = 0
    af = 0
    out = [numpy.zeros((4,5)), 0]
    with open(fin, 'rb') as csvin:
        csvreader = csv.reader(csvin, delimiter=',')
        for row in csvreader:
            if ct == 0:
                gene_id = row[id2ind['tid']].split('-')[0]
                targct = int(row[id2ind['targct']])
                out[1] = targct          

            ct+=1
            if (ct % 100) == 0:
                print ct 

            #remove on-targets
            not_ontarg = int(row[id2ind['not_ontarg']])

            #increment counts for each level of MM
            if not_ontarg:
                if row[id2ind['src']] == 'alt':
                    alt_ct+=1
                elif row[id2ind['src']] == 'ref':
                    ref_ct+=1

                nmm = int(row[id2ind['oMM']])  

                af = float(row[id2ind['ac']])/float(row[id2ind['an']])
                if af >= float(0.1):
                    out[0][0,nmm]+=1
                    out[0][0,4]+=1
                elif af >= float(0.01):
                    out[0][1,nmm]+=1
                    out[0][1,4]+=1
                elif af >= float(0.001):
                    out[0][2,nmm]+=1
                    out[0][2,4]+=1
                elif af >= float(0.0001):
                    out[0][3,nmm]+=1
                    out[0][3,4]+=1 
                else:
                    print 'error: freq lower than 1/namect'
                    print exit(1)                                       

    #write data structure for input into R
    with open(fout,'wb') as csvout:
        mywriter = csv.writer(csvout, delimiter=',')
        out[0] = out[0].astype(float)/float(out[1])            
        for i in range(0,4):
            mywriter.writerow(list(out[0][i,:])+[gene_id,prot_id,abin[i]])                                      

fin = sys.argv[1]
fout = sys.argv[2]

main(fin,fout)


