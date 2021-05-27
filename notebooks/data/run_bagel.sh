#!/usr/bin/env

mkdir data/BAGEL_output 

for i in data/CFGs/*
do 
	b=$(echo $i | sed 's/\//\t/g' | rev | cut -f1 | rev | sed 's/\.txt//')
	mkdir data/BAGEL_output/$b

	for j in data/CRISPRcleaned_logFC/*
	do
		a=$(echo $j | sed 's/\//\t/g' | rev | cut -f1 | rev | sed 's/\.gz//')

		python3 data/BAGEL.py bf -i $j -o data/BAGEL_output/$b/$a -e $i -n data/curated_NEG.txt -c 1 -s 1234 -r
	done
done
