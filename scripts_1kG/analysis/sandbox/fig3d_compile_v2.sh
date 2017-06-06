#!/bin/bash

for fin in `ls ../**_src-hom.csv`; do
    echo $fin
    fout=${fin%.csv*}"_fig3d.csv"
    echo $fout

    python fig3d_compile_v2.py $fin ../info/ALL.phase3.names.txt $fout
    # break
done

cat ../**_search**_fig3d.csv > ../_fig3d.csv