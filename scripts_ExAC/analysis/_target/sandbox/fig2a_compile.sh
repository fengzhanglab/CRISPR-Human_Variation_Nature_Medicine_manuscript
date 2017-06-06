#!/bin/bash

fin=../_search.vrer.tq.csv
fout=../_search.vrer.tq.af.csv

python fig2a_compile.py $fin $fout

Rscript fig2a.R $fout