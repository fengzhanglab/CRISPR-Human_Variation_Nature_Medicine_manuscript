import os
import csv
import sys
import numpy

def main(fin,namesin,fout):

    print 'importing data'

    mmbin = (['0-MM','1-MM','2-MM','3-MM','all'])

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

    out = {}
    with open(namesin, 'rb') as txtin:
        for line in txtin:
            print line.rstrip()
            out[line.rstrip()] = numpy.zeros(5)

    ct = 0
    alt_ct = 0
    ref_ct = 0
    nmm = 0
    targct = 0
    names = []
    name_strip = ''
    with open(fin, 'rb') as csvin:
        csvreader = csv.reader(csvin, delimiter=',')
        for row in csvreader:
            if ct == 0:
                gene_id = row[id2ind['tid']].split('-')[0]
                targct = int(row[id2ind['targct']])        

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
                names = row[id2ind['names']].replace('[','').replace(']','').split(',')

                for name in names:
                    name_strip = name.replace('\'','').replace(' ','')
                    # print out[name_strip]
                    # print out[name_strip].astype(float)
                    out[name_strip][nmm]+=1 
                    out[name_strip][4]+=1                                      

    #write data structure for input into R
    with open(fout,'wb') as csvout:
        mywriter = csv.writer(csvout, delimiter=',')
        out2 = []
        for el in out:
            print el
            print out[el]
            out2+=[[el]+list(out[el].astype(float)/float(targct))]  

        out2.sort(key=lambda x: x[0]) 
        out2 = map(list, zip(*out2))         
        ct = -1
        for row in out2:
            if ct >= 0:
                mywriter.writerow(row+[gene_id,prot_id,mmbin[ct]])
            ct+=1                                      

fin = sys.argv[1]
namesin = sys.argv[2]
fout = sys.argv[3]

main(fin,namesin,fout)


