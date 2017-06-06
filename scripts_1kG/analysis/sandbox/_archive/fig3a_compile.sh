#!/bin/bash

PAM_param="AGG,CGG,GGG,TGG,AAG,CAG,GAG,TAG"

for f in `ls ../**_exoc.spcas9.fa`; do
    echo $f
    targct=`wc -l $f | cut -f1 -d' '`
    echo $targct
    targct=$((targct/2))
    echo $targct
    namect=`wc -l ../ALL.phase3.names.txt | cut -f1 -d' '`
    echo $namect

    IFS='_' read -r -a els <<< "$f"
    echo ${els[*]}
    gene_id=${els[2]}
    echo $gene_id
    fin="../_search-"$gene_id".csv"
    echo $fin
    annot=${f%.spcas9.fa*}".csv"
    echo $annot
    fout=${fin%.csv*}"_src.csv"
    echo $fout

    python fig3a_compile.py $fin $annot $namect $targct $fout
    
done

cat ../**_search**_fig3a.csv > ../_fig3a.csv