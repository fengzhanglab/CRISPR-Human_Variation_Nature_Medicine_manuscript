#!/bin/bash

fin="../_search-fig2xs.sacas9.gq.csv"
fout="_search-fig2xs.sacas9.gq.pdf"

cat ../_search-fig2x_**.gq.csv > $fin

Rscript fig2xs-gene_v3.R $fin $fout


