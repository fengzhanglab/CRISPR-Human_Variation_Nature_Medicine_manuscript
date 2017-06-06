#!/bin/bash

prot_id=ascpf1
PAM_param="TTTA,TTTC,TTTG,TTTT"

root_path=/broad/zhanglabdata/DAS_sandbox/_genoGen
ftp_path=$root_path"/Reference/ftp.1000genomes.ebi.ac.uk/vol1/ftp"
ref_path=$ftp_path"/technical/reference/phase2_reference_assembly_sequence/_source"
ref_genome=$ref_path'/hs37d5.fa'

# ref_genome='hs37d5.fa'

for f in `ls info/*.gqc.fa`; do
    echo $f
    temp=`wc -l $f`
    echo $temp
    IFS=' ' read -r -a els <<< "$temp"
    targct=${els[0]}
    echo $targct
    targct=$((targct/2))
    echo $targct
    names=info/ALL.phase3.names.txt
    males=info/ALL.phase3.males.txt

    IFS='_' read -r -a els <<< "$f"
    IFS='.' read -r -a els2 <<< "${els[2]}"
    echo ${els[*]}
    echo ${els2[*]}
    gene_id=${els2[0]}
    echo $gene_id
    fin="_search-"$gene_id"."$prot_id".pam.csv"
    echo $fin
    fout=${fin%.csv*}"_src.csv"
    echo $fout

    py_call="search_ref_compile_v2.py $fin $names $males $targct $ref_genome $PAM_param $fout"
    echo $py_call

    qsub -P varicas -q short -o ${fin%.csv*}".o" -e ${fin%.csv*}".e" -b y -cwd -l h_vmem=15g python $py_call
    # python $py_call

    # break
done


