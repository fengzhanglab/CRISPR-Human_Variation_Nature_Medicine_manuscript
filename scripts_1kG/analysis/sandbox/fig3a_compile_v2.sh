#!/bin/bash

for fin in `ls ../**_src-hom.csv`; do
    echo $fin
    fout=${fin%.csv*}"_fig3a.csv"
    echo $fout

    python fig3a_compile_v2.py $fin $fout
    # break
done

cat ../**_search**_fig3a.csv > ../_fig3a.csv