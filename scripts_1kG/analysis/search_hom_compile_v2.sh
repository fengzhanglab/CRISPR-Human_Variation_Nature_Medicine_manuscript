#!/bin/bash

for fin in `ls *.pam_src.csv`; do
    names=info/ALL.phase3.names.txt
    males=info/ALL.phase3.males.txt

    fout=${fin%.csv*}"-hom.csv"
    echo $fout

    py_call="search_hom_compile_v2.py $fin $names $males $fout"
    echo $py_call

    qsub -P varicas -q short -o ${fin%.csv*}".o" -e ${fin%.csv*}".e" -b y -cwd -l h_vmem=15g python $py_call
    # python $py_call

    # break
done


