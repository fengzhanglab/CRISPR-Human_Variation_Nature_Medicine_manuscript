#!/bin/bash

for fin in `ls ../**_src.csv`; do
    echo $fin
    fout=${fin%.csv*}"_fig3a-adtarg.csv"
    echo $fout

    python fig3a_compile_v2-adtarg.py $fin $fout
    # break
done

cat ../**_search**_fig3a-adtarg.csv > ../_fig3a-adtarg.csv