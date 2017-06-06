#!/bin/bash

ref=Homo_sapiens_assembly19_21.fa
vcf=ExAC.r0.3.1.sites_21.vep.vcf

# qsub -q short -o gg_build_test.o -e gg_build_test.e -b y -cwd -l m_mem_free=3G python2.7 ../geno_genome-build_v0_2_1.py build_gg $ref $vcf

python2.7 ../exac_genome-build_v0_2_1.py split_vcf $ref $vcf 5000000 1000
