#!/bin/bash

for fin in `ls ../**_src.csv`; do
    echo $fin
    fout=${fin%.csv*}"_fig3a-dtarg.csv"
    echo $fout

    python fig3a_compile_v2-dtarg.py $fin $fout
    # break
done

cat ../**_search**_fig3a-dtarg.csv > ../_fig3a-dtarg.csv