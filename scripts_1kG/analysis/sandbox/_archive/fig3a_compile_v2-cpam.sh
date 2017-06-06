#!/bin/bash

for fin in `ls ../**_src.csv`; do
    echo $fin
    fout=${fin%.csv*}"_fig3a-cpam.csv"
    echo $fout

    python fig3a_compile_v2-cpam.py $fin $fout
    # break
done

cat ../**_search**_fig3a-cpam.csv > ../_fig3a-cpam.csv