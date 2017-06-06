#!/bin/bash
#David Scott, MIT, 2016

# WD=`pwd`
# THIS=$WD/$0
THIS="/broad/zhanglabdata/DAS_sandbox/_genoGen/Release/scripts_ExAC/bin/exacGen.sh"
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

geno_gen_build="exac_genome-build_v0_2_1.py"
geno_gen_pam="exac_genome-PAM_query_v0_2_1.py"
geno_gen_target="exac_genome-target_query_v0_2_1.py"
geno_gen_platinum="exac_genome-platinum_query_v0_2_1.py"
geno_gen_gene="exac_genome-gene_query_v0_2_1.py"
geno_gen_nt="exac_genome-nt_query_v0_2_1.py"
geno_gen_search="exac_genome-search_v0_2.py"

root_path=/broad/zhanglabdata/DAS_sandbox/_genoGen
ftp_path=$root_path"/Reference/ftp.broadinstitute.org/pub"
ref_path=$ftp_path"/seq/references"
vcf_path=$ftp_path"/ExAC_release/release0.3.1"
cov_path=$ftp_path"/ExAC_release/release0.3.1/coverage"
gff3_path=$root_path"/Reference/gencode19"

ref_genome='Homo_sapiens_assembly19.fa'
ref_genome_base=${ref_genome%.*}
vcf_prefix='ExAC.r0.3.1.sites'
gff3_prefix='gencode.v19.annotation'
cov_prefix='Panel.chr'

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
    -g|--gene)
    gene_list="$2"
    echo "gene_list specified:"$gene_list
    shift # past argument
    ;;
    -p|--pam)
    pam_list="$2"
    echo "pam_list specified:"$pam_list
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
if [ $method == "split_chrom" ]; then
	todo_root="${todo%.*}"	

	rm -f $todo
	f=$ref_path/$ref_genome
	ff=$vcf_path/$vcf_prefix".vep.vcf.gz"
	fff=$gff3_path/$gff3_prefix".gff3.gz"
	echo $f $ff $fff > $todo

	echo "initializing chrom split"
	array_params="-s $THIS -k split_chrom_job -t $todo"
	cmd="$bin/parallel.ops.sh -m run_job_array $array_params --mem 10g --max_jobs 1"
	qsub -q long -o $todo_root".o" -e $todo_root".e" -cwd $cmd

##################################################
#split vcf genome file to chromosomes
##################################################
elif [ $method == "split_chrom_job" ]; then
	todo_root="${todo%.*}"
	todo_done=$todo_root"_done.txt"	

	while IFS='' read -r line || [[ -n "$line" ]]; do
		echo $line
		success=true
		{ 
			python2.7 $scripts/$geno_gen_build split_chrom $line
			#gzip $line
		} || {
			echo "split chrom job unsuccessful"
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
		# for f in `ls $ref_path/*.fa`; do
		for f in `ls $cov_path/*.txt.gz`; do		
			f_base=${f%.*}
			chrID=${f_base##*_}
			echo $chrID
			ff=$vcf_path"/"$vcf_prefix"_"$chrID".vep.vcf"
			fff=$gff3_path"/"$gff3_prefix"_"$chrID".gff3"		
			echo $f $ff $fff $f $bin_size $bin_overlap >> $todo
		done
	fi

	echo "initializing reference vcf split"
	array_params="-s $THIS -k split_vcf_job -t $todo"
	cmd="$bin/parallel.ops.sh -m run_long_job_array $array_params --mem 10g --max_jobs 50"
	qsub -q long -o $todo_root".o" -e $todo_root".e" -cwd $cmd

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
		for f in `ls $ref_path/$ref_genome_base**_**_**.fa`; do
			f_base=${f%.*}
			chrBIN=${f_base##*_}
			f_base=${f%_*}
			chrID=${f_base##*_}
			echo $chrID
			if [ -e $vcf_path"/"$vcf_prefix"_"$chrID"_"$chrBIN".vep.vcf" ]; then
				ff=$vcf_path"/"$vcf_prefix"_"$chrID"_"$chrBIN".vep.vcf"
				fff=$cov_path"/"$cov_prefix$chrID".coverage_"$chrBIN".txt"
				echo $f
				echo $ff			
				echo $f $ff $fff >> $todo
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
	cmd="$bin/parallel.ops.sh -m run_job_array $array_params --mem 10g --max_jobs 2000"
	qsub -q long -o $todo_root".o" -e $todo_root".e" -cwd $cmd

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
elif [ $method == "query_spcas9_pam" ]; then
	todo_root="${todo%.*}"

	PAM_param="AGG,CGG,GGG,TGG,AAG,CAG,GAG,TAG"
	if [ ! -f $todo ]; then
		for f in `ls $ref_path/$ref_genome_base**_**_**.gg`; do
			f_base=${f%.gg*}
			chrBIN=${f_base##*_}
			f_base=${f%_*}
			chrID=${f_base##*_}
			echo $chrID
			if [ -e $gff3_path"/"$gff3_prefix"_"$chrID"_"$chrBIN".gff3" ]; then
				ff=$gff3_path"/"$gff3_prefix"_"$chrID"_"$chrBIN".gff3"
				echo $f
				echo $ff			
				echo $f $ff $PAM_param ${f%.gg*}".spcas9.pq.csv" >> $todo
			else
				echo "NO GFF3 FOR GG: "
			fi
		done
	fi

	echo "initializing pam query"
	array_params="-s $THIS -k query_pam_job -t $todo"
	# cmd="$bin/parallel.ops.sh -m run_long_job_array $array_params --mem 30g --max_jobs 50"
	cmd="$bin/parallel.ops.sh -m run_job_array $array_params --mem 6g --max_jobs 2000"
	echo $cmd
	qsub -q long -o $todo_root".o" -e $todo_root".e" -cwd $cmd

##################################################
#parallel_ops wrapper
#compile cas9 pam objects
##################################################
elif [ $method == "query_spcas9NGG_pam" ]; then
	todo_root="${todo%.*}"

	PAM_param="AGG,CGG,GGG,TGG"
	if [ ! -f $todo ]; then
		for f in `ls $ref_path/$ref_genome_base**_**_**.gg`; do
			f_base=${f%.gg*}
			chrBIN=${f_base##*_}
			f_base=${f%_*}
			chrID=${f_base##*_}
			echo $chrID
			if [ -e $gff3_path"/"$gff3_prefix"_"$chrID"_"$chrBIN".gff3" ]; then
				ff=$gff3_path"/"$gff3_prefix"_"$chrID"_"$chrBIN".gff3"
				echo $f
				echo $ff			
				echo $f $ff $PAM_param ${f%.gg*}".spcas9NGG.pq.csv" >> $todo
			else
				echo "NO GFF3 FOR GG: "
			fi
		done
	fi

	echo "initializing pam query"
	array_params="-s $THIS -k query_pam_job -t $todo"
	# cmd="$bin/parallel.ops.sh -m run_long_job_array $array_params --mem 30g --max_jobs 50"
	cmd="$bin/parallel.ops.sh -m run_job_array $array_params --mem 6g --max_jobs 2000"
	echo $cmd
	qsub -q long -o $todo_root".o" -e $todo_root".e" -cwd $cmd

##################################################
#parallel_ops wrapper
#compile cas9 pam objects
##################################################
elif [ $method == "query_vrer_pam" ]; then
	todo_root="${todo%.*}"

	PAM_param="AGCG,CGCG,GGCG,TGCG"
	if [ ! -f $todo ]; then
		for f in `ls $ref_path/$ref_genome_base**_**_**.gg`; do
			f_base=${f%.gg*}
			chrBIN=${f_base##*_}
			f_base=${f%_*}
			chrID=${f_base##*_}
			echo $chrID
			if [ -e $gff3_path"/"$gff3_prefix"_"$chrID"_"$chrBIN".gff3" ]; then
				ff=$gff3_path"/"$gff3_prefix"_"$chrID"_"$chrBIN".gff3"
				echo $f
				echo $ff			
				echo $f $ff $PAM_param ${f%.gg*}".vrer.pq.csv" >> $todo
			else
				echo "NO GFF3 FOR GG: "
			fi
		done
	fi

	echo "initializing pam query"
	array_params="-s $THIS -k query_pam_job -t $todo"
	# cmd="$bin/parallel.ops.sh -m run_long_job_array $array_params --mem 30g --max_jobs 50"
	cmd="$bin/parallel.ops.sh -m run_job_array $array_params --mem 6g --max_jobs 2000"
	echo $cmd
	qsub -q long -o $todo_root".o" -e $todo_root".e" -cwd $cmd	

##################################################
#parallel_ops wrapper
#compile cas9 pam objects
##################################################
elif [ $method == "query_vqr_pam" ]; then
	todo_root="${todo%.*}"

	PAM_param="AGA,CGA,GGA,TGA"
	if [ ! -f $todo ]; then
		for f in `ls $ref_path/$ref_genome_base**_**_**.gg`; do
			f_base=${f%.gg*}
			chrBIN=${f_base##*_}
			f_base=${f%_*}
			chrID=${f_base##*_}
			echo $chrID
			if [ -e $gff3_path"/"$gff3_prefix"_"$chrID"_"$chrBIN".gff3" ]; then
				ff=$gff3_path"/"$gff3_prefix"_"$chrID"_"$chrBIN".gff3"
				echo $f
				echo $ff			
				echo $f $ff $PAM_param ${f%.gg*}".vqr.pq.csv" >> $todo
			else
				echo "NO GFF3 FOR GG: "
			fi
		done
	fi

	echo "initializing pam query"
	array_params="-s $THIS -k query_pam_job -t $todo"
	# cmd="$bin/parallel.ops.sh -m run_long_job_array $array_params --mem 30g --max_jobs 50"
	cmd="$bin/parallel.ops.sh -m run_job_array $array_params --mem 6g --max_jobs 2000"
	echo $cmd
	qsub -q long -o $todo_root".o" -e $todo_root".e" -cwd $cmd	

##################################################
#parallel_ops wrapper
#compile cas9 pam objects
##################################################
elif [ $method == "query_sacas9_pam" ]; then
	todo_root="${todo%.*}"

	PAM_NNGRRT="AAGAAT,ACGAAT,AGGAAT,ATGAAT,AAGAGT,ACGAGT,AGGAGT,ATGAGT"
	PAM_NNGRRT=$PAM_NNGRRT",AAGGAT,ACGGAT,AGGGAT,ATGGAT,AAGGGT,ACGGGT,AGGGGT,ATGGGT"
	PAM_NNGRRT=$PAM_NNGRRT",CAGAAT,CCGAAT,CGGAAT,CTGAAT,CAGAGT,CCGAGT,CGGAGT,CTGAGT"
	PAM_NNGRRT=$PAM_NNGRRT",CAGGAT,CCGGAT,CGGGAT,CTGGAT,CAGGGT,CCGGGT,CGGGGT,CTGGGT"
	PAM_NNGRRT=$PAM_NNGRRT",GAGAAT,GCGAAT,GGGAAT,GTGAAT,GAGAGT,GCGAGT,GGGAGT,GTGAGT"
	PAM_NNGRRT=$PAM_NNGRRT",GAGGAT,GCGGAT,GGGGAT,GTGGAT,GAGGGT,GCGGGT,GGGGGT,GTGGGT"
	PAM_NNGRRT=$PAM_NNGRRT",TAGAAT,TCGAAT,TGGAAT,TTGAAT,TAGAGT,TCGAGT,TGGAGT,TTGAGT"
	PAM_param=$PAM_NNGRRT",TAGGAT,TCGGAT,TGGGAT,TTGGAT,TAGGGT,TCGGGT,TGGGGT,TTGGGT"
	if [ ! -f $todo ]; then
		for f in `ls $ref_path/$ref_genome_base**_**_**.gg`; do
			f_base=${f%.gg*}
			chrBIN=${f_base##*_}
			f_base=${f%_*}
			chrID=${f_base##*_}
			echo $chrID
			if [ -e $gff3_path"/"$gff3_prefix"_"$chrID"_"$chrBIN".gff3" ]; then
				ff=$gff3_path"/"$gff3_prefix"_"$chrID"_"$chrBIN".gff3"
				echo $f
				echo $ff			
				echo $f $ff $PAM_param ${f%.gg*}".sacas9.pq.csv" >> $todo
			else
				echo "NO GFF3 FOR GG: "
			fi
		done
	fi

	echo "initializing pam query"
	array_params="-s $THIS -k query_pam_job -t $todo"
	# cmd="$bin/parallel.ops.sh -m run_long_job_array $array_params --mem 30g --max_jobs 50"
	cmd="$bin/parallel.ops.sh -m run_job_array $array_params --mem 6g --max_jobs 2000"
	echo $cmd
	qsub -q long -o $todo_root".o" -e $todo_root".e" -cwd $cmd	

##################################################
#parallel_ops wrapper
#compile cas9 pam objects
##################################################
elif [ $method == "query_lbcpf1_pam" ]; then
	todo_root="${todo%.*}"

	PAM_param="TTTA,TTTC,TTTG,TTTT"
	if [ ! -f $todo ]; then
		for f in `ls $ref_path/$ref_genome_base**_**_**.gg`; do
			f_base=${f%.gg*}
			chrBIN=${f_base##*_}
			f_base=${f%_*}
			chrID=${f_base##*_}
			echo $chrID
			if [ -e $gff3_path"/"$gff3_prefix"_"$chrID"_"$chrBIN".gff3" ]; then
				ff=$gff3_path"/"$gff3_prefix"_"$chrID"_"$chrBIN".gff3"
				echo $f
				echo $ff			
				echo $f $ff $PAM_param ${f%.gg*}".lbcpf1.pq.csv" >> $todo
			else
				echo "NO GFF3 FOR GG: "
			fi
		done
	fi

	echo "initializing pam query"
	array_params="-s $THIS -k query_pam_job -t $todo"
	# cmd="$bin/parallel.ops.sh -m run_long_job_array $array_params --mem 30g --max_jobs 50"
	cmd="$bin/parallel.ops.sh -m run_job_array $array_params --mem 6g --max_jobs 2000"
	echo $cmd
	qsub -q long -o $todo_root".o" -e $todo_root".e" -cwd $cmd

##################################################
#parallel_ops wrapper
#compile cas9 pam objects
##################################################
elif [ $method == "query_2nt_pam" ]; then
	if [ -z "$pam_list" ]; then
		echo "pam_list unset, exiting"
		exit 0
	fi

	todo_root="${todo%.*}"

	PAM_param="A${pam_list},C${pam_list},G${pam_list},T${pam_list}"
	echo $PAM_param
	if [ ! -f $todo ]; then
		mkdir _search
		for f in `ls $ref_path/$ref_genome_base**_**_**.gg`; do
			f_base=${f%.gg*}
			chrBIN=${f_base##*_}
			f_base=${f%_*}
			chrID=${f_base##*_}
			echo $chrID
			if [ -e $gff3_path"/"$gff3_prefix"_"$chrID"_"$chrBIN".gff3" ]; then
				ff=$gff3_path"/"$gff3_prefix"_"$chrID"_"$chrBIN".gff3"
				echo $f
				echo $ff
				f_base=`basename $f`
				echo "_search/"${f_base%.gg*}".N${pam_list}.pq.csv"			
				echo $f $ff $PAM_param "_search/"${f_base%.gg*}".N${pam_list}.pq.csv" >> $todo
			else
				echo "NO GFF3 FOR GG: "
			fi
		done
	fi

	echo "initializing pam query"
	array_params="-s $THIS -k query_pam_job -t $todo"
	# cmd="$bin/parallel.ops.sh -m run_long_job_array $array_params --mem 30g --max_jobs 50"
	cmd="$bin/parallel.ops.sh -m run_job_array $array_params --mem 6g --max_jobs 2000"
	echo $cmd
	# qsub -q long -o $todo_root".o" -e $todo_root".e" -cwd $cmd
	qsub -q long -o $todo_root".o" -e $todo_root".e" -cwd $cmd

##################################################
#compile cas9 pam object
##################################################
elif [ $method == "query_pam_job" ]; then
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
elif [ $method == "query_spcas9_target" ]; then
	todo_root="${todo%.*}"

	PAM_param="20 R AGG,CGG,GGG,TGG,AAG,CAG,GAG,TAG"
	if [ ! -f $todo ]; then
		for f in `ls $ref_path/$ref_genome_base**_**_**.gg`; do
			f_base=${f%.gg*}
			chrBIN=${f_base##*_}
			f_base=${f%_*}
			chrID=${f_base##*_}
			echo $chrID
			if [ -e $gff3_path"/"$gff3_prefix"_"$chrID"_"$chrBIN".gff3" ]; then
				ff=$gff3_path"/"$gff3_prefix"_"$chrID"_"$chrBIN".gff3"
				echo $f
				echo $ff			
				echo $f $ff $PAM_param ${f%.gg*}".spcas9.tq.csv" >> $todo
			else
				echo "NO GFF3 FOR GG: "
			fi
		done
	fi

	echo "initializing target query"
	array_params="-s $THIS -k query_target_job -t $todo"
	# cmd="$bin/parallel.ops.sh -m run_long_job_array $array_params --mem 30g --max_jobs 50"
	cmd="$bin/parallel.ops.sh -m run_job_array $array_params --mem 6g --max_jobs 2000"
	echo $cmd
	qsub -q long -o $todo_root".o" -e $todo_root".e" -cwd $cmd

##################################################
#parallel_ops wrapper
#compile cas9 pam objects
##################################################
elif [ $method == "query_spcas9NGG_target" ]; then
	todo_root="${todo%.*}"

	PAM_param="20 R AGG,CGG,GGG,TGG"
	if [ ! -f $todo ]; then
		for f in `ls $ref_path/$ref_genome_base**_**_**.gg`; do
			f_base=${f%.gg*}
			chrBIN=${f_base##*_}
			f_base=${f%_*}
			chrID=${f_base##*_}
			echo $chrID
			if [ -e $gff3_path"/"$gff3_prefix"_"$chrID"_"$chrBIN".gff3" ]; then
				ff=$gff3_path"/"$gff3_prefix"_"$chrID"_"$chrBIN".gff3"
				echo $f
				echo $ff			
				echo $f $ff $PAM_param ${f%.gg*}".spcas9NGG.tq.csv" >> $todo
			else
				echo "NO GFF3 FOR GG: "
			fi
		done
	fi

	echo "initializing target query"
	array_params="-s $THIS -k query_target_job -t $todo"
	# cmd="$bin/parallel.ops.sh -m run_long_job_array $array_params --mem 30g --max_jobs 50"
	cmd="$bin/parallel.ops.sh -m run_job_array $array_params --mem 6g --max_jobs 2000"
	echo $cmd
	qsub -q long -o $todo_root".o" -e $todo_root".e" -cwd $cmd

##################################################
#parallel_ops wrapper
#compile cas9 pam objects
##################################################
elif [ $method == "query_vrer_target" ]; then
	todo_root="${todo%.*}"

	PAM_param="20 R AGCG,CGCG,GGCG,TGCG"
	if [ ! -f $todo ]; then
		for f in `ls $ref_path/$ref_genome_base**_**_**.gg`; do
			f_base=${f%.gg*}
			chrBIN=${f_base##*_}
			f_base=${f%_*}
			chrID=${f_base##*_}
			echo $chrID
			if [ -e $gff3_path"/"$gff3_prefix"_"$chrID"_"$chrBIN".gff3" ]; then
				ff=$gff3_path"/"$gff3_prefix"_"$chrID"_"$chrBIN".gff3"
				echo $f
				echo $ff			
				echo $f $ff $PAM_param ${f%.gg*}".vrer.tq.csv" >> $todo
			else
				echo "NO GFF3 FOR GG: "
			fi
		done
	fi

	echo "initializing target query"
	array_params="-s $THIS -k query_target_job -t $todo"
	# cmd="$bin/parallel.ops.sh -m run_long_job_array $array_params --mem 30g --max_jobs 50"
	cmd="$bin/parallel.ops.sh -m run_job_array $array_params --mem 6g --max_jobs 2000"
	echo $cmd
	qsub -q long -o $todo_root".o" -e $todo_root".e" -cwd $cmd

##################################################
#parallel_ops wrapper
#compile cas9 pam objects
##################################################
elif [ $method == "query_vqr_target" ]; then
	todo_root="${todo%.*}"

	PAM_param="20 R AGA,CGA,GGA,TGA"
	if [ ! -f $todo ]; then
		for f in `ls $ref_path/$ref_genome_base**_**_**.gg`; do
			f_base=${f%.gg*}
			chrBIN=${f_base##*_}
			f_base=${f%_*}
			chrID=${f_base##*_}
			echo $chrID
			if [ -e $gff3_path"/"$gff3_prefix"_"$chrID"_"$chrBIN".gff3" ]; then
				ff=$gff3_path"/"$gff3_prefix"_"$chrID"_"$chrBIN".gff3"
				echo $f
				echo $ff			
				echo $f $ff $PAM_param ${f%.gg*}".vqr.tq.csv" >> $todo
			else
				echo "NO GFF3 FOR GG: "
			fi
		done
	fi

	echo "initializing target query"
	array_params="-s $THIS -k query_target_job -t $todo"
	# cmd="$bin/parallel.ops.sh -m run_long_job_array $array_params --mem 30g --max_jobs 50"
	cmd="$bin/parallel.ops.sh -m run_job_array $array_params --mem 6g --max_jobs 2000"
	echo $cmd
	qsub -q long -o $todo_root".o" -e $todo_root".e" -cwd $cmd

##################################################
#parallel_ops wrapper
#compile cas9 pam objects
##################################################
elif [ $method == "query_sacas9_target" ]; then
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
		for f in `ls $ref_path/$ref_genome_base**_**_**.gg`; do
			f_base=${f%.gg*}
			chrBIN=${f_base##*_}
			f_base=${f%_*}
			chrID=${f_base##*_}
			echo $chrID
			if [ -e $gff3_path"/"$gff3_prefix"_"$chrID"_"$chrBIN".gff3" ]; then
				ff=$gff3_path"/"$gff3_prefix"_"$chrID"_"$chrBIN".gff3"
				echo $f
				echo $ff			
				echo $f $ff $PAM_param ${f%.gg*}".sacas9.tq.csv" >> $todo
			else
				echo "NO GFF3 FOR GG: "
			fi
		done
	fi

	echo "initializing target query"
	array_params="-s $THIS -k query_target_job -t $todo"
	# cmd="$bin/parallel.ops.sh -m run_long_job_array $array_params --mem 30g --max_jobs 50"
	cmd="$bin/parallel.ops.sh -m run_job_array $array_params --mem 6g --max_jobs 2000"
	echo $cmd
	qsub -q long -o $todo_root".o" -e $todo_root".e" -cwd $cmd

##################################################
#parallel_ops wrapper
#compile cas9 pam objects
##################################################
elif [ $method == "query_lbcpf1_target" ]; then
	todo_root="${todo%.*}"

	PAM_param="20 L TTTA,TTTC,TTTG,TTTT"
	if [ ! -f $todo ]; then
		for f in `ls $ref_path/$ref_genome_base**_**_**.gg`; do
			f_base=${f%.gg*}
			chrBIN=${f_base##*_}
			f_base=${f%_*}
			chrID=${f_base##*_}
			echo $chrID
			if [ -e $gff3_path"/"$gff3_prefix"_"$chrID"_"$chrBIN".gff3" ]; then
				ff=$gff3_path"/"$gff3_prefix"_"$chrID"_"$chrBIN".gff3"
				echo $f
				echo $ff			
				echo $f $ff $PAM_param ${f%.gg*}".lbcpf1.tq.csv" >> $todo
			else
				echo "NO GFF3 FOR GG: "
			fi
		done
	fi

	echo "initializing target query"
	array_params="-s $THIS -k query_target_job -t $todo"
	# cmd="$bin/parallel.ops.sh -m run_long_job_array $array_params --mem 30g --max_jobs 50"
	cmd="$bin/parallel.ops.sh -m run_job_array $array_params --mem 6g --max_jobs 2000"
	echo $cmd
	qsub -q long -o $todo_root".o" -e $todo_root".e" -cwd $cmd

##################################################
#compile cas9 pam object
##################################################
elif [ $method == "query_target_job" ]; then
	todo_root="${todo%.*}"
	todo_done=$todo_root"_done.txt"	

	while IFS='' read -r line || [[ -n "$line" ]]; do
		echo $line
		success=true
		{ 
			python2.7 $scripts/$geno_gen_target $line
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
elif [ $method == "query_spcas9NGG_platinum" ]; then
	todo_root="${todo%.*}"

	PAM_param="20 R AGG,CGG,GGG,TGG"
	if [ ! -f $todo ]; then
		mkdir _search
		for f in `ls $ref_path/$ref_genome_base**_**_**.gg`; do
			f_base=${f%.gg*}
			chrBIN=${f_base##*_}
			f_base=${f%_*}
			chrID=${f_base##*_}
			echo $chrID
			if [ -e $gff3_path"/"$gff3_prefix"_"$chrID"_"$chrBIN".gff3" ]; then
				ff=$gff3_path"/"$gff3_prefix"_"$chrID"_"$chrBIN".gff3"
				echo $f
				echo $ff
				f_base=`basename $f`			
				echo $f $ff $PAM_param "_search/"${f_base%.gg*}".spcas9NGG.tq.csv" >> $todo
			else
				echo "NO GFF3 FOR GG: "
			fi
		done
	fi

	echo "initializing platinum query"
	array_params="-s $THIS -k query_platinum_job -t $todo"
	# cmd="$bin/parallel.ops.sh -m run_long_job_array $array_params --mem 30g --max_jobs 50"
	cmd="$bin/parallel.ops.sh -m run_job_array $array_params --mem 6g --max_jobs 2000"
	echo $cmd
	qsub -q long -o $todo_root".o" -e $todo_root".e" -cwd $cmd

##################################################
#parallel_ops wrapper
#compile cas9 pam objects
##################################################
elif [ $method == "query_vqr_platinum" ]; then
	todo_root="${todo%.*}"

	PAM_param="20 R AGA,CGA,GGA,TGA"
	if [ ! -f $todo ]; then
		mkdir _search
		for f in `ls $ref_path/$ref_genome_base**_**_**.gg`; do
			f_base=${f%.gg*}
			chrBIN=${f_base##*_}
			f_base=${f%_*}
			chrID=${f_base##*_}
			echo $chrID
			if [ -e $gff3_path"/"$gff3_prefix"_"$chrID"_"$chrBIN".gff3" ]; then
				ff=$gff3_path"/"$gff3_prefix"_"$chrID"_"$chrBIN".gff3"
				echo $f
				echo $ff			
				f_base=`basename $f`			
				echo $f $ff $PAM_param "_search/"${f_base%.gg*}".vqr.tq.csv" >> $todo
			else
				echo "NO GFF3 FOR GG: "
			fi
		done
	fi

	echo "initializing platinum query"
	array_params="-s $THIS -k query_platinum_job -t $todo"
	# cmd="$bin/parallel.ops.sh -m run_long_job_array $array_params --mem 30g --max_jobs 50"
	cmd="$bin/parallel.ops.sh -m run_job_array $array_params --mem 6g --max_jobs 2000"
	echo $cmd
	qsub -q long -o $todo_root".o" -e $todo_root".e" -cwd $cmd

##################################################
#parallel_ops wrapper
#compile cas9 pam objects
##################################################
elif [ $method == "query_sacas9_platinum" ]; then
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
		mkdir _search
		for f in `ls $ref_path/$ref_genome_base**_**_**.gg`; do
			f_base=${f%.gg*}
			chrBIN=${f_base##*_}
			f_base=${f%_*}
			chrID=${f_base##*_}
			echo $chrID
			if [ -e $gff3_path"/"$gff3_prefix"_"$chrID"_"$chrBIN".gff3" ]; then
				ff=$gff3_path"/"$gff3_prefix"_"$chrID"_"$chrBIN".gff3"
				echo $f
				echo $ff			
				f_base=`basename $f`			
				echo $f $ff $PAM_param "_search/"${f_base%.gg*}".sacas9.tq.csv" >> $todo
			else
				echo "NO GFF3 FOR GG: "
			fi
		done
	fi

	echo "initializing platinum query"
	array_params="-s $THIS -k query_platinum_job -t $todo"
	# cmd="$bin/parallel.ops.sh -m run_long_job_array $array_params --mem 30g --max_jobs 50"
	cmd="$bin/parallel.ops.sh -m run_job_array $array_params --mem 6g --max_jobs 2000"
	echo $cmd
	qsub -q long -o $todo_root".o" -e $todo_root".e" -cwd $cmd

##################################################
#parallel_ops wrapper
#compile cas9 pam objects
##################################################
elif [ $method == "query_lbcpf1_platinum" ]; then
	todo_root="${todo%.*}"

	PAM_param="20 L TTTA,TTTC,TTTG,TTTT"
	if [ ! -f $todo ]; then
		mkdir _search
		for f in `ls $ref_path/$ref_genome_base**_**_**.gg`; do
			f_base=${f%.gg*}
			chrBIN=${f_base##*_}
			f_base=${f%_*}
			chrID=${f_base##*_}
			echo $chrID
			if [ -e $gff3_path"/"$gff3_prefix"_"$chrID"_"$chrBIN".gff3" ]; then
				ff=$gff3_path"/"$gff3_prefix"_"$chrID"_"$chrBIN".gff3"
				echo $f
				echo $ff			
				f_base=`basename $f`			
				echo $f $ff $PAM_param "_search/"${f_base%.gg*}".lbcpf1.tq.csv" >> $todo
			else
				echo "NO GFF3 FOR GG: "
			fi
		done
	fi

	echo "initializing platinum query"
	array_params="-s $THIS -k query_platinum_job -t $todo"
	# cmd="$bin/parallel.ops.sh -m run_long_job_array $array_params --mem 30g --max_jobs 50"
	cmd="$bin/parallel.ops.sh -m run_job_array $array_params --mem 6g --max_jobs 2000"
	echo $cmd
	qsub -q long -o $todo_root".o" -e $todo_root".e" -cwd $cmd

##################################################
#compile cas9 pam object
##################################################
elif [ $method == "query_platinum_job" ]; then
	todo_root="${todo%.*}"
	todo_done=$todo_root"_done.txt"	

	while IFS='' read -r line || [[ -n "$line" ]]; do
		echo $line
		success=true
		{ 
			python2.7 $scripts/$geno_gen_platinum $line
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
elif [ $method == "query_spcas9_gene" ]; then
	if [ -z "$gene_list" ]; then
		echo "gene_list unset, exiting"
		exit 0
	fi	

	todo_root="${todo%.*}"

	PAM_param="20 R AGG,CGG,GGG,TGG,AAG,CAG,GAG,TAG"
	if [ ! -f $todo ]; then
		for f in `ls $ref_path/$ref_genome_base**_**_**.gg`; do
			f_base=${f%.gg*}
			chrBIN=${f_base##*_}
			f_base=${f%_*}
			chrID=${f_base##*_}
			echo $chrID
			if [ -e $gff3_path"/"$gff3_prefix"_"$chrID"_"$chrBIN".gff3" ]; then
				ff=$gff3_path"/"$gff3_prefix"_"$chrID"_"$chrBIN".gff3"
				echo $f
				echo $ff
				echo $gene_list			
				echo $f $ff $gene_list $PAM_param ${f%.gg*}".spcas9.gq.csv" >> $todo
			else
				echo "NO GFF3 FOR GG: "
			fi
		done
	fi

	echo "initializing target query"
	array_params="-s $THIS -k query_gene_job -t $todo"
	# cmd="$bin/parallel.ops.sh -m run_long_job_array $array_params --mem 30g --max_jobs 50"
	cmd="$bin/parallel.ops.sh -m run_job_array $array_params --mem 6g --max_jobs 2000"
	echo $cmd
	qsub -q long -o $todo_root".o" -e $todo_root".e" -cwd $cmd

##################################################
#parallel_ops wrapper
#compile cas9 pam objects
##################################################
elif [ $method == "query_spcas9NGG_gene" ]; then
	if [ -z "$gene_list" ]; then
		echo "gene_list unset, exiting"
		exit 0
	fi	

	todo_root="${todo%.*}"

	PAM_param="20 R AGG,CGG,GGG,TGG"
	if [ ! -f $todo ]; then
		for f in `ls $ref_path/$ref_genome_base**_**_**.gg`; do
			f_base=${f%.gg*}
			chrBIN=${f_base##*_}
			f_base=${f%_*}
			chrID=${f_base##*_}
			echo $chrID
			if [ -e $gff3_path"/"$gff3_prefix"_"$chrID"_"$chrBIN".gff3" ]; then
				ff=$gff3_path"/"$gff3_prefix"_"$chrID"_"$chrBIN".gff3"
				echo $f
				echo $ff
				echo $gene_list			
				echo $f $ff $gene_list $PAM_param ${f%.gg*}".spcas9NGG.gq.csv" >> $todo
			else
				echo "NO GFF3 FOR GG: "
			fi
		done
	fi

	echo "initializing target query"
	array_params="-s $THIS -k query_gene_job -t $todo"
	# cmd="$bin/parallel.ops.sh -m run_long_job_array $array_params --mem 30g --max_jobs 50"
	cmd="$bin/parallel.ops.sh -m run_job_array $array_params --mem 6g --max_jobs 2000"
	echo $cmd
	qsub -q long -o $todo_root".o" -e $todo_root".e" -cwd $cmd

##################################################
#parallel_ops wrapper
#compile cas9 pam objects
##################################################
elif [ $method == "query_vqr_gene" ]; then
	if [ -z "$gene_list" ]; then
		echo "gene_list unset, exiting"
		exit 0
	fi	

	todo_root="${todo%.*}"

	PAM_param="20 R AGA,CGA,GGA,TGA"
	if [ ! -f $todo ]; then
		for f in `ls $ref_path/$ref_genome_base**_**_**.gg`; do
			f_base=${f%.gg*}
			chrBIN=${f_base##*_}
			f_base=${f%_*}
			chrID=${f_base##*_}
			echo $chrID
			if [ -e $gff3_path"/"$gff3_prefix"_"$chrID"_"$chrBIN".gff3" ]; then
				ff=$gff3_path"/"$gff3_prefix"_"$chrID"_"$chrBIN".gff3"
				echo $f
				echo $ff
				echo $gene_list			
				echo $f $ff $gene_list $PAM_param ${f%.gg*}".vqr.gq.csv" >> $todo
			else
				echo "NO GFF3 FOR GG: "
			fi
		done
	fi

	echo "initializing target query"
	array_params="-s $THIS -k query_gene_job -t $todo"
	# cmd="$bin/parallel.ops.sh -m run_long_job_array $array_params --mem 30g --max_jobs 50"
	cmd="$bin/parallel.ops.sh -m run_job_array $array_params --mem 6g --max_jobs 2000"
	echo $cmd
	qsub -q long -o $todo_root".o" -e $todo_root".e" -cwd $cmd

##################################################
#parallel_ops wrapper
#compile cas9 pam objects
##################################################
elif [ $method == "query_sacas9_gene" ]; then
	if [ -z "$gene_list" ]; then
		echo "gene_list unset, exiting"
		exit 0
	fi	

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
		for f in `ls $ref_path/$ref_genome_base**_**_**.gg`; do
			f_base=${f%.gg*}
			chrBIN=${f_base##*_}
			f_base=${f%_*}
			chrID=${f_base##*_}
			echo $chrID
			if [ -e $gff3_path"/"$gff3_prefix"_"$chrID"_"$chrBIN".gff3" ]; then
				ff=$gff3_path"/"$gff3_prefix"_"$chrID"_"$chrBIN".gff3"
				echo $f
				echo $ff
				echo $gene_list			
				echo $f $ff $gene_list $PAM_param ${f%.gg*}".sacas9.gq.csv" >> $todo
			else
				echo "NO GFF3 FOR GG: "
			fi
		done
	fi

	echo "initializing target query"
	array_params="-s $THIS -k query_gene_job -t $todo"
	# cmd="$bin/parallel.ops.sh -m run_long_job_array $array_params --mem 30g --max_jobs 50"
	cmd="$bin/parallel.ops.sh -m run_job_array $array_params --mem 6g --max_jobs 2000"
	echo $cmd
	qsub -q long -o $todo_root".o" -e $todo_root".e" -cwd $cmd

##################################################
#parallel_ops wrapper
#compile cas9 pam objects
##################################################
elif [ $method == "query_lbcpf1_gene" ]; then
	if [ -z "$gene_list" ]; then
		echo "gene_list unset, exiting"
		exit 0
	fi	

	todo_root="${todo%.*}"

	PAM_param="20 L TTTA,TTTC,TTTG,TTTT"
	if [ ! -f $todo ]; then
		for f in `ls $ref_path/$ref_genome_base**_**_**.gg`; do
			f_base=${f%.gg*}
			chrBIN=${f_base##*_}
			f_base=${f%_*}
			chrID=${f_base##*_}
			echo $chrID
			if [ -e $gff3_path"/"$gff3_prefix"_"$chrID"_"$chrBIN".gff3" ]; then
				ff=$gff3_path"/"$gff3_prefix"_"$chrID"_"$chrBIN".gff3"
				echo $f
				echo $ff
				echo $gene_list			
				echo $f $ff $gene_list $PAM_param ${f%.gg*}".lbcpf1.gq.csv" >> $todo
			else
				echo "NO GFF3 FOR GG: "
			fi
		done
	fi

	echo "initializing target query"
	array_params="-s $THIS -k query_gene_job -t $todo"
	# cmd="$bin/parallel.ops.sh -m run_long_job_array $array_params --mem 30g --max_jobs 50"
	cmd="$bin/parallel.ops.sh -m run_job_array $array_params --mem 6g --max_jobs 2000"
	echo $cmd
	qsub -q long -o $todo_root".o" -e $todo_root".e" -cwd $cmd

##################################################
#compile cas9 pam object
##################################################
elif [ $method == "query_gene_job" ]; then
	todo_root="${todo%.*}"
	todo_done=$todo_root"_done.txt"	

	while IFS='' read -r line || [[ -n "$line" ]]; do
		echo $line
		success=true
		{ 
			python2.7 $scripts/$geno_gen_gene $line
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
#compile nt objects
##################################################
elif [ $method == "query_nt" ]; then
	todo_root="${todo%.*}"

	if [ ! -f $todo ]; then
		for f in `ls $ref_path/$ref_genome_base**_**_**.gg`; do
			f_base=${f%.gg*}
			chrBIN=${f_base##*_}
			f_base=${f%_*}
			chrID=${f_base##*_}
			echo $chrID
			if [ -e $gff3_path"/"$gff3_prefix"_"$chrID"_"$chrBIN".gff3" ]; then
				ff=$gff3_path"/"$gff3_prefix"_"$chrID"_"$chrBIN".gff3"
				echo $f
				echo $ff			
				echo $f $ff ${f%.gg*}".ntq.csv" >> $todo
			else
				echo "NO GFF3 FOR GG: "
			fi
		done
	fi

	echo "initializing nt query"
	array_params="-s $THIS -k query_nt_job -t $todo"
	# cmd="$bin/parallel.ops.sh -m run_long_job_array $array_params --mem 30g --max_jobs 50"
	cmd="$bin/parallel.ops.sh -m run_job_array $array_params --mem 6g --max_jobs 2000"
	echo $cmd
	# qsub -q long -o $todo_root".o" -e $todo_root".e" -cwd $cmd
	qsub -q long -o $todo_root".o" -e $todo_root".e" -cwd $cmd

##################################################
#compile nt object
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
			echo "compile nt job unsuccessful"
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
elif [ $method == "search_cas9_gg" ]; then
	todo_root="${todo%.*}"

	targ="GAGTCCGAGCAGAAGAAGAA"
	PAM_param="$targ R AGG,CGG,GGG,TGG,AAG,CAG,GAG,TAG 3"
	if [ ! -f $todo ]; then
		for f in `ls $ref_path/*_**_**.gg`; do
			f_base=${f##*/}
			echo $f $PAM_param `pwd`/${f_base%.gg*}".EMX13.fa"
			echo $f $PAM_param `pwd`/${f_base%.gg*}".EMX13.fa" >> $todo
		done
	fi

	echo "initializing gg search"
	array_params="-s $THIS -k search_gg_job -t $todo"
	# cmd="$bin/parallel.ops.sh -m run_long_job_array $array_params --mem 5g --max_jobs 50"
	cmd="$bin/parallel.ops.sh -m run_job_array $array_params --mem 1G --max_jobs 2000"
	echo $cmd
	qsub -q long -o $todo_root".o" -e $todo_root".e" -cwd $cmd

##################################################
#compile cas9 pam object
##################################################
elif [ $method == "search_gg_job" ]; then
	todo_root="${todo%.*}"
	todo_done=$todo_root"_done.txt"	

	while IFS='' read -r line || [[ -n "$line" ]]; do
		echo $line
		success=true
		{ 
			python2.7 $scripts/$geno_gen_search $line
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


