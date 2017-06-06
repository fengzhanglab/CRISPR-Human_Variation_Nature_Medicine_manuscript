#!/bin/bash

ref=Homo_sapiens_assembly19_21_27972000.fa
vcf=ExAC.r0.3.1.sites_21_27972000.vep.vcf
cov=Panel.chr21.coverage_27972000.txt

qsub -q short -o gg_build_test.o -e gg_build_test.e -b y -cwd -l h_vmem=3g python2.7 ../exac_genome-build_v0_2_1.py build_gg $ref $vcf $cov

# python2.7 ../exac_genome-build_v0_2_1.py build_gg $ref $vcf
