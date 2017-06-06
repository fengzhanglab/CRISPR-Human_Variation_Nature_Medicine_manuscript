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
    rows = {}
    rowsl = []
    g_dict = {}
    with open(fin, 'rb') as csvin:
        csvreader = csv.reader(csvin, delimiter=',')
        for row in csvreader:

            ct+=1
            if (ct % 100) == 0:
                print ct

            #translate target root coord to cut site
            #Cas9, TS offset 17, BS offset 6
            if row[36] == 'TS':
                #+18 to make cut sites 5' to the break on target strand
                row[35] = int(row[35])+18
            elif row[36] == 'BS':
                row[35] = int(row[35])+9
            
            el = row[len(row)-1]
            if el not in rows:
                rows[el] = [row]
            else:
                rows[el]+=[row]

        #sort by cut sites in direction of transcription
        for el in rows:
            if rows[el][0][5] == '+':
                rows[el].sort(key=lambda x: x[35])
            elif rows[el][0][5] == '-':
                rows[el].sort(key=lambda x: x[35], reverse=True)
            rowsl+=rows[el]

        for row in rowsl:
            rowc = numpy.asarray(row[9:36],dtype=float)

            el = row[len(row)-1]
            if g_param[el][2] == row[7]:
                if el not in g_dict:
                    g_dict[el] = open(fout.replace(
                        '.'+prot_id+'.gq.csv','_'+el+'.'+prot_id+'.gq.csv'),'wb')
                    mywriter = csv.writer(g_dict[el], delimiter=',')
                    temp = [numpy.amax(rowc[0:26])]
                    mywriter.writerow( temp+[el]+row[35:37]+[row[40]]) 
                                                       
                else:
                    mywriter = csv.writer(g_dict[el], delimiter=',')
                    temp = [numpy.amax(rowc[0:26])]
                    mywriter.writerow( temp+[el]+row[35:37]+[row[40]])                                                    
    
    for el in g_dict:
        g_dict[el].close()                                                          

fin = sys.argv[1]
genes = sys.argv[2]
fout = sys.argv[3]
prot_id = sys.argv[4]

main(fin,genes,fout,prot_id)


