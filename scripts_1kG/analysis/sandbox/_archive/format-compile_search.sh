#!/bin/bash

for f in `ls _search/**.pam.csv`; do
	echo $f
	f_base=`basename $f`
	tr -d $'\r' < $f > .tmp.csv
	sed -e "s/$/,$f_base/" -i .tmp.csv
	cat .tmp.csv >> _search.csv		
done
