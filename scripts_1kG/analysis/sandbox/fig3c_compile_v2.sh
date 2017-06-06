#!/bin/bash

for fin in `ls ../**_src-hom.csv`; do
    echo $fin
    fout=${fin%.csv*}"_fig3c.csv"
    echo $fout

    python fig3c_compile_v2.py $fin $fout
    # break
done

cat ../**_search**_fig3c.csv > ../_fig3c.csv