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

    gene_id=${fin%.spcas9.pam_src-hom.csv*}
    gene_id=${gene_id#*_search-}
    echo $gene_id
    annot=../info/"_search-fig2xp_"$gene_id".spcas9NGG.gq.annot.csv"

    Rscript fig3xx_v4.R $fout $annot ${fout_base%.csv}".pdf"
    # break
done


