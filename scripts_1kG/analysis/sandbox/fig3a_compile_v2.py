import os
import csv
import sys
import numpy

def main(fin,fout):

    print 'importing data'

    #allele type
    atype = (['all','ref','adtarg','cpam','ddtarg','dtarg'])

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
    out = [numpy.zeros((6,5)), 0]
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

                #all: presuming less than or equal to 3MM
                #ind 4 for each row is total for <= 3MM
                nmm = int(row[id2ind['oMM']])  
                out[0][0,nmm]+=1
                out[0][0,4]+=1

                #adtarg: only increment array if:
                #off-target is present increase in distance of alt allele to targ
                #compared to distance between ref allele and targ
                if row[id2ind['src']] == 'ref':
                    out[0][1,nmm]+=1
                    out[0][1,4]+=1

                #adtarg: only increment array if:
                #off-target is present increase in distance of alt allele to targ
                #compared to distance between ref allele and targ
                if int(row[id2ind['rMM']]) < int(row[id2ind['oMM']]):
                    out[0][2,nmm]+=1
                    out[0][2,4]+=1

                #cpam: only increment array if:
                #off-target is present from pam creation in alt allele
                #where none existed in ref_allele                    
                if int(row[id2ind['rPAMflag']]) == 0:
                    out[0][3,nmm]+=1
                    out[0][3,4]+=1

                #ddtarg: only increment array if:
                #off-target is present from drop in distance of alt allele to targ
                #compared to distance between ref allele and targ
                if int(row[id2ind['rMM']]) > int(row[id2ind['oMM']]):
                    out[0][4,nmm]+=1
                    out[0][4,4]+=1

                #dtarg: only increment array if:
                #off-target is present from same distance between alt allele and targ
                #change in mismatch locations between ref allele and targ
                if ((not (row[id2ind['rtarg']] == row[id2ind['otarg']])) and
                    (int(row[id2ind['rMM']]) == int(row[id2ind['oMM']]))):
                    out[0][5,nmm]+=1
                    out[0][5,4]+=1

            else:
                print 'on-target: continuing'
    
    print 'alt_ct: '+str(alt_ct)
    print 'ref_ct: '+str(ref_ct)

    #write data structure for input into R
    with open(fout,'wb') as csvout:
        mywriter = csv.writer(csvout, delimiter=',')
        #get average number of OT at each mm level per target
        out[0] = out[0].astype(float)/float(out[1])  
        for i in range(6):           
            mywriter.writerow(list(out[0][i,:])+[gene_id,prot_id,atype[i]])                                      

fin = sys.argv[1]
fout = sys.argv[2]

main(fin,fout)


