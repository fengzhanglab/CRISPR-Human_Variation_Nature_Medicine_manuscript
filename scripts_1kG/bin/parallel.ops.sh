#!/bin/sh
##################################################
##################################################
#David Scott, MIT, 2016
#rules and behavior for array parallelization
#INPUT:
#all submissions must consist of a list of:
#1. a shell SCRIPT for individual job handling
#2. a parameter specifying the SCRIPT TASK to be run
#3. a TODO LIST of single files containing SCRIPT commands 
#	OR a list of files for manipulation by the SCRIPT
#4. a string of parameters parameterizing the TASK 
#OUTPUT:
#1. all wrapper scripts must take a INPUT FILE containing 
#	single item from TODO LIST
#2. must output a DONE FILE named "[INPUT FILE]_done.txt"
#	containing the finished lines from INPUT_FILE
#3. upon successful completion, job must output file
#	"[INPUT FILE].success" containing md5sum of DONE FILE
#BEHAVIOR:
#given satisfaction of above parameters, the general mode
#of execution will consist of a master thread that will:
#	1. dispatch job
#	2. monitor job for completion
#	3. clean up after job following completion
#REQUEUING
#jobs with incomplete execution or unexpected termination
#can be requeued by simply resubmitting the initial 
#job submission, and the array execution will continue
#where the master thread was last working
#ONLY REQUEUE FOR JOBS WITH LINEAR EXECUTION THAT WILL
#OVERWRITE ANY PREVIOUSLY GENERATED OUTPUT FOR INCOMPLETE JOBS
##################################################
##################################################

THIS="$0"
echo $THIS

USEUSE=/broad/software/scripts/useuse
source $USEUSE
use UGER

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

echo $1
echo $2

case $key in
    -m|--method)
    method="$2"
    echo "method specified:"$method
    shift # past argument
    ;;	
    -s|--script)
    script="$2"
    echo "script specified:"$script
    shift # past argument
    ;;
   	-t|--todo)
    todo="$2"
    echo "todo specified:"$todo
    shift # past argument
    ;;
   	-k|--task)
    task="$2"
    echo "task specified:"$task
    shift # past argument
    ;;    
    -p|--params)
    params="$2"
    echo "params specified:"$params
    shift # past argument
    ;;
    --mem)
    memory="$2"
    echo "memory specified:"$memory
    shift # past argument
    ;;   
    --max_jobs)
    max_jobs="$2"
    echo "max_jobs specified:"$max_jobs
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
#parallel job handler
##################################################
if [ $method == "run_job_array" ]; then
	if [ -z "$script" ]; then
		echo "script unset, exiting"
		exit 0
	elif [ -z "$task" ]; then
		echo "task unset, exiting"
		exit 0
	elif [ -z "$memory" ]; then
		echo "memory unset, exiting"
		exit 0			
	elif [ -z "$max_jobs" ]; then
		echo "max_jobs unset, exiting"
		exit 0				
	fi

	if [ -z "$params" ]; then
		echo "dispatching job array"
		job_array_param="-s $script -k $task -t $todo --mem $memory --max_jobs $max_jobs"
		bash $THIS -m dispatch_job_array $job_array_param
	else
		echo "dispatching job array with params"
		job_array_param="-s $script -k $task -t $todo --mem $memory --max_jobs $max_jobs"
		bash $THIS -m dispatch_job_array $job_array_param -p "$params"
	fi

	bash $THIS -m monitor_job_array -t $todo
	bash $THIS -m clean_job_array2 -t $todo

##################################################
#parallel long job handler
##################################################
elif [ $method == "run_long_job_array" ]; then
	if [ -z "$script" ]; then
		echo "script unset, exiting"
		exit 0
	elif [ -z "$task" ]; then
		echo "task unset, exiting"
		exit 0	
	elif [ -z "$memory" ]; then
		echo "memory unset, exiting"
		exit 0							
	elif [ -z "$max_jobs" ]; then
		echo "max_jobs unset, exiting"
		exit 0				
	fi	

	if [ -z "$params" ]; then
		echo "dispatching long job array"
		job_array_param="-s $script -k $task -t $todo --mem $memory --max_jobs $max_jobs"
		bash $THIS -m dispatch_long_job_array $job_array_param
	else
		echo "dispatching long job array with params"
		job_array_param="-s $script -k $task -t $todo --mem $memory --max_jobs $max_jobs"
		bash $THIS -m dispatch_long_job_array $job_array_param -p "$params"
	fi

	echo "monitoring long job array"
	bash $THIS -m monitor_job_array -t $todo
	bash $THIS -m clean_job_array2 -t $todo

##################################################
#parallel job dispatcher
##################################################
elif [ $method == "dispatch_job_array" ]; then
	if [ -z "$script" ]; then
		echo "script unset, exiting"
		exit 0
	elif [ -z "$task" ]; then
		echo "task unset, exiting"
		exit 0	
	elif [ -z "$memory" ]; then
		echo "memory unset, exiting"
		exit 0				
	elif [ -z "$max_jobs" ]; then
		echo "max_jobs unset, exiting"
		exit 0								
	fi	

	todo_root="${todo%.*}"
	todo_input=$todo_root"_input.txt"
	todo_array=$todo_root"_array.txt"
	todo_script=$todo_root"_array.sh"
	todo_array_root="${todo_array%.*}"
	todo_array_pend=$todo_array_root"_pend.txt"
	todo_array_done=$todo_array_root"_done.txt"

	cp $todo $todo_input
	if [ -e $todo_array_pend ]; then
		echo "restarting job for todo_array:"
		echo "path:$todo_array"
		bash $THIS -m snap_job_array -t $todo
		bash $THIS -m clean_job_array2 -t $todo

		#setup truncated todo list to requeue
		rm -f $todo
		while IFS='' read -r line || [[ -n "$line" ]]; do
			line_root="${line%.*}"
			line_done=$line_root"_done.txt"			
			line_temp=$line_root"_temp.txt"

			sort $line > $line_temp
			cat $line_temp > $line
			sort $line_done > $line_temp
			cat $line_temp > $line_done
			comm -3 $line $line_done > $line_temp					

			rm -f $line
			while IFS='' read -r line_el || [[ -n "$line_el" ]]; do
				echo $line_el | sed -e 's/^[[:space:]]*//' >> $line
			done < "$line_temp"
			cat $line >> $todo
			rm $line
			rm $line_done
			rm $line_temp
		done < "$todo_array_pend"
		rm -f $todo_array_temp
		rm -f $todo_array_done
	else
		echo "starting new job for todo:"
		echo "path:$todo"
	fi
	nline=$(wc -l < $todo)

	if (($max_jobs >= $nline)); then
		todo_job_size=1
		awk_param="-v job_root=$todo_root -v job_size=$todo_job_size"
		awk $awk_param '{filename = job_root"_job_"int(((NR)/job_size))".txt"; print >> filename}' $todo		
	else
		todo_job_size=$(((nline/max_jobs)+1))
		awk_param="-v job_root=$todo_root -v job_size=$todo_job_size"
		awk $awk_param '{filename = job_root"_job_"int(((NR)/job_size)+1)".txt"; print >> filename}' $todo		
	fi

	#split current searches for current iteration

	njob=0
	for f in `ls $todo_root"_job_"**".txt"`; do
		echo $f >> $todo_array
		((njob+=1))
	done

	#parameterize job_array
	echo '#!/bin/bash' >> $todo_script
	echo '#$'" -P varicas" >> $todo_script
	echo '#$'" -l h_vmem=$memory" >> $todo_script
	echo '#$ -N job_array' >> $todo_script
	echo '#$ -t'" 1-$njob:1" >> $todo_script
	echo '#$ -cwd' >> $todo_script
	if [ -z "$params" ]; then
		echo "$script -m $task -t $todo_root"'_job_${SGE_TASK_ID}.txt' >> $todo_script
	else
		echo "$script -m $task -t $todo_root"'_job_${SGE_TASK_ID}.txt'" $params" >> $todo_script
	fi

	#submit job array and monitor array completion
	qsub $todo_script

	cp $todo_input $todo
	rm -f $todo_input

##################################################
#parallel job long dispatcher
##################################################
elif [ $method == "dispatch_long_job_array" ]; then
	if [ -z "$script" ]; then
		echo "script unset, exiting"
		exit 0
	elif [ -z "$task" ]; then
		echo "task unset, exiting"
		exit 0	
	elif [ -z "$memory" ]; then
		echo "memory unset, exiting"
		exit 0		
	elif [ -z "$max_jobs" ]; then
		echo "max_jobs unset, exiting"
		exit 0								
	fi

	todo_root="${todo%.*}"
	todo_input=$todo_root"_input.txt"
	todo_array=$todo_root"_array.txt"
	todo_script=$todo_root"_array.sh"
	todo_array_root="${todo_array%.*}"
	todo_array_pend=$todo_array_root"_pend.txt"
	todo_array_done=$todo_array_root"_done.txt"

	cp $todo $todo_input
	if [ -e $todo_array_pend ]; then
		echo "restarting job for todo_array:"
		echo "path:$todo_array"
		bash $THIS -m snap_job_array -t $todo
		bash $THIS -m clean_job_array2 -t $todo

		#setup truncated todo list to requeue
		rm -f $todo
		while IFS='' read -r line || [[ -n "$line" ]]; do
			line_root="${line%.*}"
			line_done=$line_root"_done.txt"			
			line_temp=$line_root"_temp.txt"

			sort $line > $line_temp
			cat $line_temp > $line
			sort $line_done > $line_temp
			cat $line_temp > $line_done
			comm -3 $line $line_done > $line_temp					

			rm -f $line
			while IFS='' read -r line_el || [[ -n "$line_el" ]]; do
				echo $line_el | sed -e 's/^[[:space:]]*//' >> $line
			done < "$line_temp"
			cat $line >> $todo
			rm $line
			rm $line_done
			rm $line_temp
		done < "$todo_array_pend"
		rm -f $todo_array_temp
		rm -f $todo_array_done	
	else
		echo "starting new job for todo:"
		echo "path:$todo"
	fi
	nline=$(wc -l < $todo)

	if (($max_jobs >= $nline)); then
		todo_job_size=1
		awk_param="-v job_root=$todo_root -v job_size=$todo_job_size"
		awk $awk_param '{filename = job_root"_job_"int(((NR)/job_size))".txt"; print >> filename}' $todo		
	else
		todo_job_size=$(((nline/max_jobs)+1))
		awk_param="-v job_root=$todo_root -v job_size=$todo_job_size"
		awk $awk_param '{filename = job_root"_job_"int(((NR)/job_size)+1)".txt"; print >> filename}' $todo		
	fi

	for f in `ls $todo_root"_job_"**".txt"`; do
		fe=${f%.*}".e"
		fo=${f%.*}".o"
		echo $f >> $todo_array

		qsub_param="-q long -P varicas -e $fe -o $fo -cwd -l h_vmem=$memory"
		if [ -z "$params" ]; then
			qsub $qsub_param $script -m $task -t $f
		else
			qsub $qsub_param $script -m $task -t $f $params
		fi
	done

	cp $todo_input $todo
	rm -f $todo_input

##################################################
#parallel job monitor
##################################################
elif [ $method == "monitor_job_array" ]; then

	todo_root="${todo%.*}"
	todo_array=$todo_root"_array.txt"
	todo_array_root="${todo_array%.*}"
	todo_array_pend=$todo_array_root"_pend.txt"
	todo_array_done=$todo_array_root"_done.txt"
	todo_array_temp=$todo_array_root"_temp.txt"

	cp $todo_array $todo_array_pend		
	while [ -s $todo_array_pend ]; do
		bash $THIS -m snap_job_array -t $todo
		sleep 30
	done
	rm -f $todo_array_temp
	rm -f $todo_array_done

##################################################
#parallel job state snapshot
##################################################
elif [ $method == "snap_job_array" ]; then

	todo_root="${todo%.*}"
	todo_array=$todo_root"_array.txt"
	todo_array_root="${todo_array%.*}"
	todo_array_pend=$todo_array_root"_pend.txt"
	todo_array_done=$todo_array_root"_done.txt"
	todo_array_temp=$todo_array_root"_temp.txt"	

	cat $todo_array > $todo_array_temp
	sort -u $todo_array_temp > $todo_array

	echo "uger queue status at:"`date`
	qstat
	while IFS='' read -r line || [[ -n "$line" ]]; do
		line_root="${line%.*}"
		if [ -s $line_root"_done.txt" ]; then
			md5_str=`md5sum $line`
			md5_line=${md5_str%$line*}
			md5_str=`md5sum $line_root"_done.txt"`
			md5_done=${md5_str%$line_root"_done.txt"*}
			if [ "$md5_line" == "$md5_done" ]; then
				echo $md5_done > $line_root".success"
				echo $line >> $todo_array_done
			fi	
		fi	
	done < "$todo_array_pend"
	cat $todo_array_done > $todo_array_temp
	sort -u $todo_array_temp > $todo_array_done
	echo "jobs finished at:"`date`
	cat $todo_array_done

	rm -f $todo_array_pend
	comm -3 $todo_array $todo_array_done > $todo_array_temp
	while IFS='' read -r line || [[ -n "$line" ]]; do
		echo $line | sed -e 's/^[[:space:]]*//' >> $todo_array_pend
	done < "$todo_array_temp"
	echo "jobs remianing at:"`date`
	cat $todo_array_pend

##################################################
#parallel job cleaner
##################################################
elif [ $method == "clean_job_array" ]; then

	todo_root="${todo%.*}"
	todo_array=$todo_root"_array.txt"
	todo_array_root="${todo_array%.*}"
	echo "removing:"$todo_array_root".sh"
	rm -f $todo_array_root".sh"
	while IFS='' read -r line || [[ -n "$line" ]]; do
		line_root="${line%.*}"
		if [ -s $line_root".success" ]; then
			echo "removing:"$line
			echo "removing:"$line_root"_done.txt"
			echo "removing:"$line_root".success"				
			rm -f $line
			rm -f $line_root"_done.txt"
			rm -f $line_root".success"
		else
			echo "removing:"$line_root"_done.txt"
			rm -f $line_root"_done.txt"
		fi	
	done < $todo_array
	rm -f $todo_array

##################################################
#parallel job cleaner2
##################################################
elif [ $method == "clean_job_array2" ]; then

	todo_root="${todo%.*}"
	todo_array=$todo_root"_array.txt"
	todo_array_root="${todo_array%.*}"
	echo "removing:"$todo_array_root".sh"
	rm -f $todo_array_root".sh"
	while IFS='' read -r line || [[ -n "$line" ]]; do
		line_root="${line%.*}"
		if [ -s $line_root".success" ]; then
			echo "removing:"$line
			echo "removing:"$line_root"_done.txt"
			echo "removing:"$line_root".success"				
			rm -f $line
			rm -f $line_root"_done.txt"
			rm -f $line_root".success"
		fi	
	done < $todo_array
	rm -f $todo_array

fi

echo "Done!"


