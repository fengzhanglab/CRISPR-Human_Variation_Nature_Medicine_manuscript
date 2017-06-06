#!/bin/bash

prot_id="sacas9"
fin='../_search.'$prot_id'.gq.csv'
targ=../../../ensembl/GRCh37/annot/therapeutic_targets_v2.csv
fout='../_search-fig2xp.'$prot_id'.gq.csv'
fout_root=${fout%'.'$prot_id'.gq.csv'*}

echo $fout_root
rm ../_search-fig2xp_*
python fig2xp-gene_compile_v3.py $fin $targ $fout $prot_id

for f in `ls $fout_root"_"**.gq.fa`; do
	head -n 200 $f > ${f%'.'$prot_id'.gq.fa'*}".gqc.fa"
done

rm _search-fig2xp_*
for f in `ls $fout_root"_"**.gq.csv`; do
	f_base=`basename $f`
	Rscript fig2xp-gene_v3.R $f ${f_base%'.'$prot_id'.gq.csv'*}".pdf"
done


