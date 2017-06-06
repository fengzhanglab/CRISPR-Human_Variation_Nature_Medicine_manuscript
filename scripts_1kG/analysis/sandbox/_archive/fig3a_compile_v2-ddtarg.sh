#!/bin/bash

for fin in `ls ../**_src.csv`; do
    echo $fin
    fout=${fin%.csv*}"_fig3a-ddtarg.csv"
    echo $fout

    python fig3a_compile_v2-ddtarg.py $fin $fout
    # break
done

cat ../**_search**_fig3a-ddtarg.csv > ../_fig3a-ddtarg.csv