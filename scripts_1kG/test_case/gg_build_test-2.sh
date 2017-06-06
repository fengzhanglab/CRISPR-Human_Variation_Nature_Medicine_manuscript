#!/bin/bash

ref=2_1089000.fa
vcf=ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_1089000.vcf

qsub -q short -o gg_build_test.o -e gg_build_test.e -b y -cwd -l m_mem_free=3G python2.7 ../geno_genome-build_v0_2_1.py build_gg $ref $vcf

# python2.7 ../geno_genome-PAM_compile_v0_2_1.py $f $PAM_param $fout 
