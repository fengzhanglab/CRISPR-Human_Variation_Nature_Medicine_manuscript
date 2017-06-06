#!/bin/bash

f=Homo_sapiens_assembly19_21_27972000.gg
ff=gencode.v19.annotation_21_27972000.gff3
f_base=`basename $f`
fout=${f_base%.gg*}".ntq.csv"

call="../exac_genome-nt_query_v0_2_1.py $f $ff $fout"

qsub -P varicas -q short -o nt_query_test.o -e nt_query_test.e -b y -cwd -l h_vmem=3g python2.7 $call


