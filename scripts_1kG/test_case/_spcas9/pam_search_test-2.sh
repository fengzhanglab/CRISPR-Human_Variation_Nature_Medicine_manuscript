#!/bin/bash

threshold=5
f=2_1089000.spcas9.pam
targ=Homo_sapiens_EMX1_001_sequence_exoc.spcas9.fa
fout=$f".csv"

qsub -q short -o pam_search_test.o -e pam_search_test.e -b y -cwd -l m_mem_free=3G python2.7 ../geno_PAM-search_v0_2_1.py $f $targ $threshold $fout

# python2.7 ../geno_PAM-search_v0_2_1.py $f $targ $threshold $fout
