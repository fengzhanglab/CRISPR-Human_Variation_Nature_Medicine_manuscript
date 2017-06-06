#!/bin/bash

# for fin in `ls ../**_src-hom.csv`; do
#     echo $fin
#     fout=${fin%.csv*}"_fig3x.csv"
#     echo $fout

#     python fig3x_compile_v2.py $fin $fout
#     # break
# done

for fin in `ls ../**_src-hom.csv`; do
    echo $fin
    fout=${fin%.csv*}"_fig3x.csv"
    fout_base=`basename $fout`

    Rscript fig3x_v2.R $fout ${fout_base%.csv}".pdf"
    # break
done

# cat ../**_search**_fig3x-ddtarg-dtarg-cpam.csv > ../_fig3x-ddtarg-dtarg-cpam.csv