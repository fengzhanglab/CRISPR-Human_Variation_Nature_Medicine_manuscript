import os
import csv
import sys
import numpy

def main(fin,fout):

    print 'importing data'

    #ids = (['chrID','src','ftype','lb','ub','dot1','str','dot2','info',...
    #Import the sample data and concatenate
    ct = 0
    g_dict = {}
    with open(fin, 'rb') as csvin:
        csvreader = csv.reader(csvin, delimiter=',')
        for row in csvreader:
            rowc = row[9:32]

            # #demonstrates that by excluding regions with
            # #no variant call in exac, rates increase 2-3X
            # #most of these null rows are due to no coverage by exac
            # #thf, motivates filtering gencode exons for exac coverage
            # #see 160915 todo
            # rowsum = float(0)
            # for el in rowc:
            #     rowsum+=float(el)
            # if not (rowsum > 0):
            #     # print 'no var for row'
            #     continue 

            ct+=1
            if (ct % 100) == 0:
                print ct 
            if (row[37] == 'ALL') and ('transcript_type=protein_coding' in row[8]):
                el = row[len(row)-1]
                # if el not in g_dict:
                #     g_dict[el] = [el,numpy.zeros((6,23)),0]
                #     g_dict[el][2]+=1
                #     for i in range(0,len(rowc)):
                #         if float(rowc[i]) >= 0.1:
                #             g_dict[el][1][0,i]+=float(rowc[i])
                #         elif float(rowc[i]) >= 0.01:
                #             g_dict[el][1][1,i]+=float(rowc[i])
                #         elif float(rowc[i]) >= 0.001:
                #             g_dict[el][1][2,i]+=float(rowc[i])
                #         elif float(rowc[i]) >= 0.0001:
                #             g_dict[el][1][3,i]+=float(rowc[i])
                #         elif float(rowc[i]) >= 0.00001:
                #             g_dict[el][1][4,i]+=float(rowc[i])
                #         elif float(rowc[i]) > 0.000001:
                #             g_dict[el][1][5,i]+=float(rowc[i])
                # else:
                #     g_dict[el][2]+=1
                #     for i in range(0,len(rowc)):
                #         if float(rowc[i]) >= 0.1:
                #             g_dict[el][1][0,i]+=float(rowc[i])
                #         elif float(rowc[i]) >= 0.01:
                #             g_dict[el][1][1,i]+=float(rowc[i])
                #         elif float(rowc[i]) >= 0.001:
                #             g_dict[el][1][2,i]+=float(rowc[i])
                #         elif float(rowc[i]) >= 0.0001:
                #             g_dict[el][1][3,i]+=float(rowc[i])
                #         elif float(rowc[i]) >= 0.00001:
                #             g_dict[el][1][4,i]+=float(rowc[i])
                #         elif float(rowc[i]) > 0.000001:
                #             g_dict[el][1][5,i]+=float(rowc[i])

                if el not in g_dict:
                    g_dict[el] = [el,numpy.zeros((6,23)),0]
                    g_dict[el][2]+=1
                    for i in range(0,len(rowc)):
                        if float(rowc[i]) >= 0.1:
                            g_dict[el][1][0,i]+=1
                        if float(rowc[i]) >= 0.01:
                            g_dict[el][1][1,i]+=1
                        if float(rowc[i]) >= 0.001:
                            g_dict[el][1][2,i]+=1
                        if float(rowc[i]) >= 0.0001:
                            g_dict[el][1][3,i]+=1
                        if float(rowc[i]) >= 0.00001:
                            g_dict[el][1][4,i]+=1
                        if float(rowc[i]) >= 0.000001:
                            g_dict[el][1][5,i]+=1
                else:
                    g_dict[el][2]+=1
                    for i in range(0,len(rowc)):
                        if float(rowc[i]) >= 0.1:
                            g_dict[el][1][0,i]+=1
                        if float(rowc[i]) >= 0.01:
                            g_dict[el][1][1,i]+=1
                        if float(rowc[i]) >= 0.001:
                            g_dict[el][1][2,i]+=1
                        if float(rowc[i]) >= 0.0001:
                            g_dict[el][1][3,i]+=1
                        if float(rowc[i]) >= 0.00001:
                            g_dict[el][1][4,i]+=1
                        if float(rowc[i]) >= 0.000001:
                            g_dict[el][1][5,i]+=1               

    #query_gene automatically orients
    #no need for reorientation below
    # for g in g_dict:
    #      numpy.fliplr(g_dict[g][1])           

    for el in g_dict:
        with open(fout.replace('.spcas9NGG.gq.csv','_'+el+'.spcas9NGG.gq.csv'),'wb') as csvout:
            mywriter = csv.writer(csvout, delimiter=',')        
            if g_dict[el][2] > 0:
                g_dict[el][1] = g_dict[el][1].astype(float)/float(g_dict[el][2])            

            for i in range(0,6):
                mywriter.writerow(list(g_dict[el][1][i,:])+[el, i]) 

    bins = [10,100,1000,10000,100000,1000000]
    with open(fout.replace('.spcas9NGG.gq.csv','c.spcas9NGG.gq.csv'),'wb') as csvout:
        mywriter = csv.writer(csvout, delimiter=',')  
        for el in g_dict:
            mywriter.writerow([0,0,0, el, 0])                  
            for i in range(0,6):
                temp = ([numpy.sum(g_dict[el][1][i,0:10],dtype=float),
                    numpy.sum(g_dict[el][1][i,10:20],dtype=float),
                    numpy.sum(g_dict[el][1][i,20:23],dtype=float)])
                mywriter.writerow(temp+[el, bins[i]])                                                      

fin = sys.argv[1]
fout = sys.argv[2]

main(fin,fout)


