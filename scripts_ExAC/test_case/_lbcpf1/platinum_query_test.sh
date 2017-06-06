#!/bin/bash

PAM_param="20 L TTTA,TTTC,TTTG,TTTT"
f=../Homo_sapiens_assembly19_21_27972000.gg
ff=../gencode.v19.annotation_21_27972000.gff3
f_base=`basename $f`
fout=${f_base%.gg*}".lbcpf1.ptq.csv"

call="../../exac_genome-platinum_query_v0_2_1.py $f $ff $PAM_param $fout"

qsub -q short -o platinum_query_test.o -e platinum_query_test.e -b y -cwd -l h_vmem=3g python2.7 $call

# python2.7 ../exac_genome-target_query_v0_2_1.py $f $PAM_param $fout
