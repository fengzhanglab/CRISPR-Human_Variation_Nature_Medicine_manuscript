import os
import csv
import sys
import numpy

def main(fin,genes,fout,prot_id):

    print 'importing data'

    #ids = (['chrID','src','ftype','lb','ub','dot1','str','dot2','info',...
    #Import the sample data and concatenate
    g_param = {}
    with open(genes,'rb') as csvgene:
        greader = csv.reader(csvgene, delimiter=',')
        for g in greader:
            g_param[g[0]] = g

    ct = 0
    prows = {}
    prowsl = []
    g_dict = {}
    a_dict = {}
    with open(fin, 'rb') as csvin:
        csvreader = csv.reader(csvin, delimiter=',')
        for row in csvreader:
            rowc = numpy.asarray(row[9:35],dtype=float)
            rowmax = numpy.amax(rowc)

            ct+=1
            if (ct % 100) == 0:
                print ct

            if not (row[40] == 'ALL'):
                continue

            #platinum threshold is 1e-4
            if rowmax < 1e-4:
                #translate target root coord to cut site
                #Cas9, TS offset 17, BS offset 6
                if row[36] == 'TS':
                    #+18 to make cut sites 5' to the break on target strand
                    row[35] = int(row[35])+18
                elif row[36] == 'BS':
                    row[35] = int(row[35])+9
                
                el = row[len(row)-1]
                if el not in prows:
                    prows[el] = [row]
                else:
                    prows[el]+=[row]

        #sort by cut sites in direction of transcription
        for el in prows:
            if prows[el][0][5] == '+':
                prows[el].sort(key=lambda x: x[35])
            elif prows[el][0][5] == '-':
                prows[el].sort(key=lambda x: x[35], reverse=True)
            prowsl+=prows[el]

        for row in prowsl:
            rowc = numpy.asarray(row[9:35],dtype=float)

            el = row[len(row)-1]
            if g_param[el][2] == row[7]:
                if el not in g_dict:
                    g_dict[el] = [open(fout.replace(
                        '.'+prot_id+'.gq.csv','_'+el+'.'+prot_id+'.gq.csv'),'wb')]
                    g_dict[el]+=[open(fout.replace(
                        '.'+prot_id+'.gq.csv','_'+el+'.'+prot_id+'.gq.fa'),'wb')]
                    mywriter = csv.writer(g_dict[el][0], delimiter=',')
                    temp = [numpy.amax(rowc[0:26])]
                    mywriter.writerow(list(temp)+[el]+row[35:37]+[row[40]]) 

                    trans_str = str(row[42])+'-'+str(row[41])
                    loc_str = str(row[35])
                    str_str = str(row[36])
                    pam_str = str(row[37][20:26])
                    targ_str = str(row[37][0:20])
                    chr_str = str(row[0].replace('chr',''))
                    head_str = '>'+trans_str+'|'+chr_str+'|'+loc_str+'|'+pam_str+'|'+'R|'+str_str  
                    g_dict[el][1].write(head_str+'\n')
                    g_dict[el][1].write(targ_str+'\n')   

                    a_dict[el] = {}
                    annot = str(row[0:9])
                    if annot not in a_dict[el]:
                        print row[0:9]
                        a_dict[el][annot] = row[0:9]                                                         
                else:
                    mywriter = csv.writer(g_dict[el][0], delimiter=',')
                    temp = [numpy.amax(rowc[0:26])]
                    mywriter.writerow(list(temp)+[el]+row[35:37]+[row[40]]) 

                    trans_str = str(row[42])+'-'+str(row[41])
                    loc_str = str(row[35])
                    str_str = str(row[36])
                    pam_str = str(row[37][20:26])
                    targ_str = str(row[37][0:20])
                    chr_str = str(row[0].replace('chr',''))
                    head_str = '>'+trans_str+'|'+chr_str+'|'+loc_str+'|'+pam_str+'|'+'R|'+str_str  
                    g_dict[el][1].write(head_str+'\n')
                    g_dict[el][1].write(targ_str+'\n')   

                    annot = str(row[0:9])
                    if annot not in a_dict[el]:
                        print row[0:9]
                        a_dict[el][annot] = row[0:9]                                                    
    
    for el in g_dict:
        g_dict[el][0].close() 
        g_dict[el][1].close()  

    for el in a_dict:
        annots = a_dict[el]
        with open(fout.replace(
            '.'+prot_id+'.gq.csv','_'+el+'.'+prot_id+'.gq.annot.csv'),'wb') as csvout:
            csvwriter = csv.writer(csvout, delimiter=',')
            for annot in annots:
                csvwriter.writerow(annots[annot])                                                           

fin = sys.argv[1]
genes = sys.argv[2]
fout = sys.argv[3]
prot_id = sys.argv[4]

main(fin,genes,fout,prot_id)


