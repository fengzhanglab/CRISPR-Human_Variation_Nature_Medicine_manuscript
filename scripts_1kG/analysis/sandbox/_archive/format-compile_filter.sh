#!/bin/bash

ct=0
while IFS='' read -r line || [[ -n "$line" ]]; do
	IFS=',' read -r -a els <<< "$line"
	nmm=${els[7]}
	if ((nmm <= 2)); then
		echo $line >> _search_le2mm.csv
	fi
	echo $ct
	((ct+=1))
done < "_search.csv"
