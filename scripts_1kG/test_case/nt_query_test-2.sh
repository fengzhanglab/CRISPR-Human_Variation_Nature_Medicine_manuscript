#!/bin/bash

f=2_1089000.gg
f_base=`basename $f`
fout=${f_base%.gg*}".ntq.csv"

qsub -q short -o nt_query_test.o -e nt_query_test.e -b y -cwd -l h_vmem=3G python2.7 ../geno_genome-nt_query_v0_2_1.py $f $fout


