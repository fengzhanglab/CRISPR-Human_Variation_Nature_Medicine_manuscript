import os
import csv
import sys
import numpy

def main(fin,fout):

    print 'importing data'

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
    out = {}
    idstr = ''
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
                idstr = (row[id2ind['tchr']]+"|"+row[id2ind['tloc']]+"|"+
                    row[id2ind['tPAMori']]+"|"+row[id2ind['tstr']])

                if idstr not in out:
                    out[idstr] = ([numpy.zeros(5,dtype=float),row[id2ind['tchr']],
                        row[id2ind['tloc']],row[id2ind['tPAMori']],row[id2ind['tstr']]])

                out[idstr][0][nmm]+=1 
                out[idstr][0][4]+=1                                       

    #write data structure for input into R
    with open(fout,'wb') as csvout:
        mywriter = csv.writer(csvout, delimiter=',')           
        for el in out:
            mywriter.writerow(list(out[el][0])+out[el][1:5]+[gene_id,prot_id])

fin = sys.argv[1]
fout = sys.argv[2]

main(fin,fout)


