#!/bin/bash
#David Scott, MIT, 2016

# WD=`pwd`
# THIS=$WD/$0
THIS="/broad/zhanglabdata/DAS_sandbox/_genoGen/Release/scripts_1kG/bin/genoGen.sh"
bin=${THIS%/*}
scripts=${bin%/*}
echo $THIS
echo $bin
echo $scripts

USEUSE=/broad/software/scripts/useuse
source $USEUSE
use UGER
use Samtools
use Python-2.7

geno_gen_build="geno_genome-build_v0_2_1.py"
geno_gen_pam="geno_genome-PAM_compile_v0_2_1.py"
geno_gen_nt="geno_genome-nt_query_v0_2_1.py"
geno_pam_search="geno_PAM-search_v0_2_1.py"
geno_ref_search="geno_REF-search_v0_2_1.py"

root_path=/broad/zhanglabdata/DAS_sandbox/_genoGen
ftp_path=$root_path"/Reference/ftp.1000genomes.ebi.ac.uk/vol1/ftp"
vcf_path=$ftp_path"/release/20130502"
ref_path=$ftp_path"/technical/reference/phase2_reference_assembly_sequence"

ref_genome='hs37d5.fa'
vcf_prefix='ALL.chr'

bin_size=1000000
bin_overlap=1000

##################################################
#options
##################################################
# Use > 1 to consume two arguments per pass in the loop (e.g. each
# argument has a corresponding value to go with it).
# Use > 0 to consume one or more arguments per pass in the loop (e.g.
# some arguments don't have a corresponding value to go with it such
# as in the --default example).
# note: if this is set to > 0 the /etc/hosts part is not recognized ( may be a bug )
while [[ $# > 1 ]]
do
key="$1"

case $key in
    -m|--method)
    method="$2"
    echo "method specified:"$method
    shift # past argument
    ;;
    -t|--todo)
    todo="$2"
    echo "todo specified:"$todo
    shift # past argument
    ;;
    --targ)
    targ="$2"
    echo "target_fasta specified:"$targ
    shift # past argument
    ;;
    --spec)
    spec="$2"
    echo "spec specified:"$spec
    shift # past argument
    ;;        
    *)
	echo "unknown option:"$1
	echo "exiting"
	exit 0
    ;;
esac
shift # past argument or value
done

if [ -z "$method" ]; then
	echo "method unset, exiting"
	exit 0
elif [ -z "$todo" ]; then
	echo "todo unset, exiting"
	exit 0
fi	

echo $method
##################################################
#parallel_ops wrapper
#split reference genome file to chromosomes
##################################################
if [ $method == "split_ref" ]; then
	todo_root="${todo%.*}"	

	rm -f $todo
	echo $ref_path/$ref_genome > $todo

	echo "initializing reference split"
	array_params="-s $THIS -k split_ref_job -t $todo"
	cmd="$bin/parallel.ops.sh -m run_job_array $array_params --mem 2g --max_jobs 1"
	qsub -P varicas -q long -o $todo_root".o" -e $todo_root".e" -cwd $cmd

##################################################
#split vcf genome file to chromosomes
##################################################
elif [ $method == "split_ref_job" ]; then
	todo_root="${todo%.*}"
	todo_done=$todo_root"_done.txt"	

	while IFS='' read -r line || [[ -n "$line" ]]; do
		echo $line
		success=true
		{ 
			python2.7 $scripts/$geno_gen_build split_ref $line
			#gzip $line
		} || {
			echo "split ref job unsuccessful"
			success=false
			exit 0
		}	
		if $success; then
			echo $line >> $todo_done
		fi			
	done < "$todo"

##################################################
#parallel_ops wrapper
#split vcf files
##################################################
elif [ $method == "split_vcf" ]; then
	todo_root="${todo%.*}"

	if [ ! -f $todo ]; then
		for f in `ls $ref_path/*.fa`; do
			f_base=${f%.*}
			chrID=${f_base##*_}
			echo $chrID
			for ff in `ls $vcf_path"/"$vcf_prefix$chrID"."**".vcf.gz"`; do
				echo $f
				echo $ff			
				echo $f $ff $bin_size $bin_overlap >> $todo
			done
		done
	fi

	echo "initializing reference vcf split"
	array_params="-s $THIS -k split_vcf_job -t $todo"
	cmd="$bin/parallel.ops.sh -m run_long_job_array $array_params --mem 10g --max_jobs 50"
	qsub -P varicas -q long -o $todo_root".o" -e $todo_root".e" -cwd $cmd

##################################################
#split vcf files
##################################################
elif [ $method == "split_vcf_job" ]; then
	todo_root="${todo%.*}"
	todo_done=$todo_root"_done.txt"	

	while IFS='' read -r line || [[ -n "$line" ]]; do
		echo $line
		success=true
		{ 
			python2.7 $scripts/$geno_gen_build split_vcf $line
		} || {
			echo "split vcf job unsuccessful"
			success=false
			exit 0
		}	
		if $success; then
			echo $line >> $todo_done
		fi			
	done < "$todo"

##################################################
#parallel_ops wrapper
#build gg objects
##################################################
elif [ $method == "build_gg" ]; then
	todo_root="${todo%.*}"

	if [ ! -f $todo ]; then
		for f in `ls $ref_path/*_**_**.fa`; do
			f_base=${f%.*}
			chrBIN=${f_base##*_}
			f_base=${f%_*}
			chrID=${f_base##*_}
			echo $chrID
			if [ -e $vcf_path"/"$vcf_prefix$chrID"."**"_"$chrBIN".vcf" ]; then
				for ff in `ls $vcf_path"/"$vcf_prefix$chrID"."**"_"$chrBIN".vcf"`; do
					echo $f
					echo $ff			
					echo $f $ff >> $todo
				done
			else
				echo "NO VCF FOR FASTA: "
				echo $f
				echo $f >> $todo
			fi
		done
	fi

	echo "initializing reference build gg"
	array_params="-s $THIS -k build_gg_job -t $todo"
	# cmd="$bin/parallel.ops.sh -m run_long_job_array $array_params --mem 5g --max_jobs 50"
	cmd="$bin/parallel.ops.sh -m run_job_array $array_params --mem 2g --max_jobs 2000"
	qsub -P varicas -q long -o $todo_root".o" -e $todo_root".e" -cwd $cmd

##################################################
#build gg objects
##################################################
elif [ $method == "build_gg_job" ]; then
	todo_root="${todo%.*}"
	todo_done=$todo_root"_done.txt"	

	while IFS='' read -r line || [[ -n "$line" ]]; do
		echo $line
		success=true
		{ 
			python2.7 $scripts/$geno_gen_build build_gg $line
		} || {
			echo "build gg job unsuccessful"
			success=false
			exit 0
		}	
		if $success; then
			echo $line >> $todo_done
		fi			
	done < "$todo"

##################################################
#parallel_ops wrapper
#compile cas9 pam objects
##################################################
elif [ $method == "compile_ascpf1_pam" ]; then
	todo_root="${todo%.*}"

	PAM_param="20 L TTTA,TTTC,TTTG,TTTT"
	if [ ! -f $todo ]; then
		for f in `ls $ref_path/*_**_**.gg`; do
			echo $f $PAM_param ${f%.gg*}".ascpf1.pam"
			echo $f $PAM_param ${f%.gg*}".ascpf1.pam" >> $todo
		done
	fi

	echo "initializing pam compile"
	array_params="-s $THIS -k compile_pam_job -t $todo"
	# cmd="$bin/parallel.ops.sh -m run_long_job_array $array_params --mem 30g --max_jobs 50"
	cmd="$bin/parallel.ops.sh -m run_job_array $array_params --mem 10g --max_jobs 2000"
	echo $cmd
	qsub -P varicas -q long -o $todo_root".o" -e $todo_root".e" -cwd $cmd

##################################################
#parallel_ops wrapper
#compile cas9 pam objects
##################################################
elif [ $method == "compile_spcas9_pam" ]; then
	todo_root="${todo%.*}"

	PAM_param="20 R AGG,CGG,GGG,TGG,AAG,CAG,GAG,TAG"
	if [ ! -f $todo ]; then
		for f in `ls $ref_path/*_**_**.gg`; do
			echo $f $PAM_param ${f%.gg*}".spcas9.pam"
			echo $f $PAM_param ${f%.gg*}".spcas9.pam" >> $todo
		done
	fi

	echo "initializing pam compile"
	array_params="-s $THIS -k compile_pam_job -t $todo"
	# cmd="$bin/parallel.ops.sh -m run_long_job_array $array_params --mem 30g --max_jobs 50"
	cmd="$bin/parallel.ops.sh -m run_job_array $array_params --mem 10g --max_jobs 2000"
	echo $cmd
	qsub -P varicas -q long -o $todo_root".o" -e $todo_root".e" -cwd $cmd

##################################################
#parallel_ops wrapper
#compile cas9 pam objects
##################################################
elif [ $method == "compile_vqr_pam" ]; then
	todo_root="${todo%.*}"

	#append NGAN PAMs
	PAM_param="20 R AGAA,CGAA,GGAA,TGAA"
	PAM_param=$PAM_param",AGAC,CGAC,GGAC,TGAC"
	PAM_param=$PAM_param",AGAG,CGAG,GGAG,TGAG"
	PAM_param=$PAM_param",AGAT,CGAT,GGAT,TGAT"
	#append rest of NGNG PAMs
	PAM_param=$PAM_param",AGCG,CGCG,GGCG,TGCG"
	PAM_param=$PAM_param",AGGG,CGGG,GGGG,TGGG"
	PAM_param=$PAM_param",AGTG,CGTG,GGTG,TGTG"
	if [ ! -f $todo ]; then
		for f in `ls $ref_path/*_**_**.gg`; do
			echo $f $PAM_param ${f%.gg*}".vqr.pam"
			echo $f $PAM_param ${f%.gg*}".vqr.pam" >> $todo
		done
	fi

	echo "initializing pam compile"
	array_params="-s $THIS -k compile_pam_job -t $todo"
	# cmd="$bin/parallel.ops.sh -m run_long_job_array $array_params --mem 30g --max_jobs 50"
	cmd="$bin/parallel.ops.sh -m run_job_array $array_params --mem 10g --max_jobs 2000"
	echo $cmd
	qsub -P varicas -q long -o $todo_root".o" -e $todo_root".e" -cwd $cmd

##################################################
#parallel_ops wrapper
#compile cas9 pam objects
##################################################
elif [ $method == "compile_sacas9_pam" ]; then
	todo_root="${todo%.*}"

	PAM_NNGRRT="AAGAAT,ACGAAT,AGGAAT,ATGAAT,AAGAGT,ACGAGT,AGGAGT,ATGAGT"
	PAM_NNGRRT=$PAM_NNGRRT",AAGGAT,ACGGAT,AGGGAT,ATGGAT,AAGGGT,ACGGGT,AGGGGT,ATGGGT"
	PAM_NNGRRT=$PAM_NNGRRT",CAGAAT,CCGAAT,CGGAAT,CTGAAT,CAGAGT,CCGAGT,CGGAGT,CTGAGT"
	PAM_NNGRRT=$PAM_NNGRRT",CAGGAT,CCGGAT,CGGGAT,CTGGAT,CAGGGT,CCGGGT,CGGGGT,CTGGGT"
	PAM_NNGRRT=$PAM_NNGRRT",GAGAAT,GCGAAT,GGGAAT,GTGAAT,GAGAGT,GCGAGT,GGGAGT,GTGAGT"
	PAM_NNGRRT=$PAM_NNGRRT",GAGGAT,GCGGAT,GGGGAT,GTGGAT,GAGGGT,GCGGGT,GGGGGT,GTGGGT"
	PAM_NNGRRT=$PAM_NNGRRT",TAGAAT,TCGAAT,TGGAAT,TTGAAT,TAGAGT,TCGAGT,TGGAGT,TTGAGT"
	PAM_NNGRRT=$PAM_NNGRRT",TAGGAT,TCGGAT,TGGGAT,TTGGAT,TAGGGT,TCGGGT,TGGGGT,TTGGGT"			
	PAM_param="20 R $PAM_NNGRRT"
	echo $PAM_param
	if [ ! -f $todo ]; then
		for f in `ls $ref_path/*_**_**.gg`; do
			echo $f $PAM_param ${f%.gg*}".sacas9.pam"
			echo $f $PAM_param ${f%.gg*}".sacas9.pam" >> $todo
		done
	fi

	echo "initializing pam compile"
	array_params="-s $THIS -k compile_pam_job -t $todo"
	# cmd="$bin/parallel.ops.sh -m run_long_job_array $array_params --mem 30g --max_jobs 50"
	cmd="$bin/parallel.ops.sh -m run_job_array $array_params --mem 10g --max_jobs 2000"
	echo $cmd
	qsub -P varicas -q long -o $todo_root".o" -e $todo_root".e" -cwd $cmd

##################################################
#compile cas9 pam object
##################################################
elif [ $method == "compile_pam_job" ]; then
	todo_root="${todo%.*}"
	todo_done=$todo_root"_done.txt"	

	while IFS='' read -r line || [[ -n "$line" ]]; do
		echo $line
		success=true
		{ 
			python2.7 $scripts/$geno_gen_pam $line
		} || {
			echo "compile pam job unsuccessful"
			success=false
			exit 0
		}	
		if $success; then
			echo $line >> $todo_done
		fi			
	done < "$todo"

##################################################
#parallel_ops wrapper
#compile cas9 pam objects
##################################################
elif [ $method == "query_nt" ]; then
	todo_root="${todo%.*}"

	if [ ! -f $todo ]; then
		for f in `ls $ref_path/*_**_**.gg`; do
			echo $f $PAM_param ${f%.gg*}".ntq.csv"
			echo $f ${f%.gg*}".ntq.csv" >> $todo
		done
	fi

	echo "initializing nt query"
	array_params="-s $THIS -k query_nt_job -t $todo"
	# cmd="$bin/parallel.ops.sh -m run_long_job_array $array_params --mem 30g --max_jobs 50"
	cmd="$bin/parallel.ops.sh -m run_job_array $array_params --mem 10g --max_jobs 2000"
	echo $cmd
	qsub -P varicas -q long -o $todo_root".o" -e $todo_root".e" -cwd $cmd

##################################################
#compile cas9 pam object
##################################################
elif [ $method == "query_nt_job" ]; then
	todo_root="${todo%.*}"
	todo_done=$todo_root"_done.txt"	

	while IFS='' read -r line || [[ -n "$line" ]]; do
		echo $line
		success=true
		{ 
			python2.7 $scripts/$geno_gen_nt $line
		} || {
			echo "compile pam job unsuccessful"
			success=false
			exit 0
		}	
		if $success; then
			echo $line >> $todo_done
		fi			
	done < "$todo"

##################################################
#parallel_ops wrapper
#compile cas9 pam objects
##################################################
elif [ $method == "search_ascpf1_pam" ]; then
	if [ -z "$targ" ]; then
		echo "targ fasta unset, exiting"
		exit 0
	elif [ -z "$spec" ]; then
		echo "spec unset, exiting"
		exit 0
	fi

	todo_root="${todo%.*}"
	threshold="5"

	mkdir -p _search
	if [ ! -f $todo ]; then
		for f in `ls $ref_path/*_**_**.ascpf1.pam`; do
			f_base=`basename $f`
			echo $f $targ $threshold "_search/"${f_base%.ascpf1.pam.csv*}"-"$spec".ascpf1.pam.csv"
			echo $f $targ $threshold "_search/"${f_base%.ascpf1.pam.csv*}"-"$spec".ascpf1.pam.csv" >> $todo
		done
	fi

	echo "initializing pam search"
	array_params="-s $THIS -k search_pam_job -t $todo"
	# cmd="$bin/parallel.ops.sh -m run_long_job_array $array_params --mem 5g --max_jobs 50"
	cmd="$bin/parallel.ops.sh -m run_job_array $array_params --mem 5G --max_jobs 2000"
	echo $cmd
	qsub -P varicas -q long -o $todo_root".o" -e $todo_root".e" -cwd $cmd

##################################################
#parallel_ops wrapper
#compile cas9 pam objects
##################################################
elif [ $method == "search_spcas9_pam" ]; then
	if [ -z "$targ" ]; then
		echo "targ fasta unset, exiting"
		exit 0
	elif [ -z "$spec" ]; then
		echo "spec unset, exiting"
		exit 0
	fi

	todo_root="${todo%.*}"
	threshold="5"

	mkdir -p _search
	if [ ! -f $todo ]; then
		for f in `ls $ref_path/*_**_**.spcas9.pam`; do
			f_base=`basename $f`
			echo $f $targ $threshold "_search/"${f_base%.spcas9.pam.csv*}"-"$spec".spcas9.pam.csv"
			echo $f $targ $threshold "_search/"${f_base%.spcas9.pam.csv*}"-"$spec".spcas9.pam.csv" >> $todo
		done
	fi

	echo "initializing pam search"
	array_params="-s $THIS -k search_pam_job -t $todo"
	# cmd="$bin/parallel.ops.sh -m run_long_job_array $array_params --mem 5g --max_jobs 50"
	cmd="$bin/parallel.ops.sh -m run_job_array $array_params --mem 5G --max_jobs 2000"
	echo $cmd
	qsub -P varicas -q long -o $todo_root".o" -e $todo_root".e" -cwd $cmd

##################################################
#parallel_ops wrapper
#compile cas9 pam objects
##################################################
elif [ $method == "search_vqr_pam" ]; then
	if [ -z "$targ" ]; then
		echo "targ fasta unset, exiting"
		exit 0
	elif [ -z "$spec" ]; then
		echo "spec unset, exiting"
		exit 0
	fi

	todo_root="${todo%.*}"
	threshold="5"

	mkdir -p _search
	if [ ! -f $todo ]; then
		for f in `ls $ref_path/*_**_**.vqr.pam`; do
			f_base=`basename $f`
			echo $f $targ $threshold "_search/"${f_base%.vqr.pam.csv*}"-"$spec".vqr.pam.csv"
			echo $f $targ $threshold "_search/"${f_base%.vqr.pam.csv*}"-"$spec".vqr.pam.csv" >> $todo
		done
	fi

	echo "initializing pam search"
	array_params="-s $THIS -k search_pam_job -t $todo"
	# cmd="$bin/parallel.ops.sh -m run_long_job_array $array_params --mem 5g --max_jobs 50"
	cmd="$bin/parallel.ops.sh -m run_job_array $array_params --mem 5G --max_jobs 2000"
	echo $cmd
	qsub -P varicas -q long -o $todo_root".o" -e $todo_root".e" -cwd $cmd

##################################################
#parallel_ops wrapper
#compile cas9 pam objects
##################################################
elif [ $method == "search_sacas9_pam" ]; then
	if [ -z "$targ" ]; then
		echo "targ fasta unset, exiting"
		exit 0
	elif [ -z "$spec" ]; then
		echo "spec unset, exiting"
		exit 0
	fi

	todo_root="${todo%.*}"
	threshold="5"

	mkdir -p _search
	if [ ! -f $todo ]; then
		for f in `ls $ref_path/*_**_**.sacas9.pam`; do
			f_base=`basename $f`
			echo $f $targ $threshold "_search/"${f_base%.sacas9.pam.csv*}"-"$spec".sacas9.pam.csv"
			echo $f $targ $threshold "_search/"${f_base%.sacas9.pam.csv*}"-"$spec".sacas9.pam.csv" >> $todo
		done
	fi

	echo "initializing pam search"
	array_params="-s $THIS -k search_pam_job -t $todo"
	# cmd="$bin/parallel.ops.sh -m run_long_job_array $array_params --mem 5g --max_jobs 50"
	cmd="$bin/parallel.ops.sh -m run_job_array $array_params --mem 5G --max_jobs 2000"
	echo $cmd
	qsub -P varicas -q long -o $todo_root".o" -e $todo_root".e" -cwd $cmd

##################################################
#compile cas9 pam object
##################################################
elif [ $method == "search_pam_job" ]; then
	todo_root="${todo%.*}"
	todo_done="$todo_root""_done.txt"	

	while IFS='' read -r line || [[ -n "$line" ]]; do
		echo "$line"
		success=true
		{ 
			python2.7 $scripts/$geno_pam_search $line
		} || {
			echo "compile pam job unsuccessful"
			success=false
			exit 0
		}	
		if $success; then
			echo $line >> $todo_done
		fi			
	done < "$todo"

##################################################
#parallel_ops wrapper
#compile cas9 pam objects
##################################################
elif [ $method == "search_spcas9_ref" ]; then
	if [ -z "$targ" ]; then
		echo "targ fasta unset, exiting"
		exit 0
	fi

	todo_root="${todo%.*}"
	threshold="5"

	PAM_param="$targ R AGG,CGG,GGG,TGG,AAG,CAG,GAG,TAG $threshold"
	if [ ! -f $todo ]; then
		for f in `ls $ref_path/*_**_**.fa`; do
			echo $f $PAM_param ${f%.fa*}".spcas9.ref.csv"
			echo $f $PAM_param ${f%.fa*}".spcas9.ref.csv" >> $todo
		done
	fi

	echo "initializing ref search"
	array_params="-s $THIS -k search_ref_job -t $todo"
	# cmd="$bin/parallel.ops.sh -m run_long_job_array $array_params --mem 30g --max_jobs 50"
	cmd="$bin/parallel.ops.sh -m run_job_array $array_params --mem 5g --max_jobs 2000"
	echo $cmd
	qsub -P varicas -q long -o $todo_root".o" -e $todo_root".e" -cwd $cmd

##################################################
#compile cas9 pam object
##################################################
elif [ $method == "search_ref_job" ]; then
	todo_root="${todo%.*}"
	todo_done=$todo_root"_done.txt"	

	while IFS='' read -r line || [[ -n "$line" ]]; do
		echo $line
		success=true
		{ 
			python2.7 $scripts/$geno_ref_search $line
		} || {
			echo "compile pam job unsuccessful"
			success=false
			exit 0
		}	
		if $success; then
			echo $line >> $todo_done
		fi			
	done < "$todo"

fi
echo "Done!"


