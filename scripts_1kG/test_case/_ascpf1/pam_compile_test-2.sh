#!/bin/bash

PAM_param="20 L TTTA,TTTC,TTTG,TTTT"
f=../2_1089000.gg
f_base=`basename $f`
fout=../${f_base%.gg*}".ascpf1.pam"

qsub -q short -o pam_compile_test.o -e pam_compile_test.e -b y -cwd -l h_vmem=3G python2.7 ../../geno_genome-PAM_compile_v0_2_1.py $f $PAM_param $fout

# python2.7 ../geno_genome-PAM_compile_v0_2_1.py $f $PAM_param $fout 
