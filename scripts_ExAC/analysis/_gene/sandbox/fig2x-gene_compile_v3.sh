#!/bin/bash

prot_id="sacas9"
fin='../_search.'$prot_id'.gq.csv'
targ=../../../ensembl/GRCh37/annot/therapeutic_targets_v2.csv
fout='../_search-fig2x.'$prot_id'.gq.csv'
fout_root=${fout%'.'$prot_id'.gq.csv'*}

echo $fout_root
rm ../_search-fig2x_*
python fig2x-gene_compile_v3.py $fin $targ $fout $prot_id

rm _search-fig2x_*
for f in `ls $fout_root"_"**`; do
	f_base=`basename $f`
	Rscript fig2x-gene_v3.R $f ${f_base%'.'$prot_id'.gq.csv'*}".pdf"
done


