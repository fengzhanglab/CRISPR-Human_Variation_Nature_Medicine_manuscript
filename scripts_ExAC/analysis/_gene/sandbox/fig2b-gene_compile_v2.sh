#!/bin/bash

fin=../_search.spcas9NGG.gq.csv
fout=../_search-fig2b.spcas9NGG.gq.csv
fout_root=${fout%.spcas9NGG.gq.csv*}

# python fig2b-gene_compile_v2.py $fin $fout

for f in `ls $fout_root"_"**`; do
	f_base=`basename $f`
	Rscript fig2b-gene_v2.R $f ${f_base%.spcas9NGG.gq.csv*}".pdf"
done

fin=$fout_root"c.spcas9NGG.gq.csv"
fout=`basename $fout`
Rscript fig2b-gene_v2c.R $fin ${fout%.spcas9NGG.gq.csv*}"c.pdf"


