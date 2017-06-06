#!/bin/bash

fin=../_search.vrer.tq.csv
fout=../_search.vrer.tq.het.csv

python fig2a_compile-het.py $fin $fout

Rscript fig2a-2.R $fout