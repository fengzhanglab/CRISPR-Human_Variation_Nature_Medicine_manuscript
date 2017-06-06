#!/bin/bash

root_path=/broad/zhanglabdata/DAS_sandbox/_genoGen
ftp_path=$root_path"/Reference/ftp.1000genomes.ebi.ac.uk/vol1/ftp"
ref_path=$ftp_path"/technical/reference/phase2_reference_assembly_sequence/_source/_source"
ref_genome=$ref_path'/hs37d5.fa'

PAM_param="AGG,CGG,GGG,TGG,AAG,CAG,GAG,TAG"

for f in `ls *_exoc.spcas9.fa`; do
    echo $f
    targct=`wc -l $f | cut -f1 -d' '`
    echo $targct
    targct=$((targct/2))
    echo $targct
    namect=`wc -l ALL.phase3.names.txt | cut -f1 -d' '`
    echo $namect

    IFS='_' read -r -a els <<< "$f"
    echo ${els[*]}
    gene_id=${els[2]}
    echo $gene_id
    fin="_search-"$gene_id".csv"
    echo $fin
    annot=${f%.spcas9.fa*}".csv"
    echo $annot
    fout=${fin%.csv*}"_src.csv"
    echo $fout

    py_call="search_ref_compile.py $fin $annot $namect $targct $ref_genome $PAM_param $fout"

    qsub -P varicas -q short -o ${fin%.csv*}".o" -e ${fin%.csv*}".e" -b y -cwd -l h_vmem=15g python $py_call
    # break
done