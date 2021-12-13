#!/bin/bash
while getopts s:D:f:S:c:w:o:M: flag
do
	case "${flag}" in
		s) sample_sec=${OPTARG};;   # for wss.pl
		# d) d_memrec=${OPTARG};;
		D) dataset_name=${OPTARG};;
		f) num_fixels=${OPTARG};;   # 0 as full
		S) num_subj=${OPTARG};;
		c) num_cores=${OPTARG};;
		w) run_where=${OPTARG};;    # "sge" or "interactive" or "vmware" or "dopamine"
		o) output_folder=${OPTARG};;  # generated in wrapper.sh
		M) run_memoryProfiler=${OPTARG};;   # TRUE or FALSE
	esac
done

date

echo "sampling every __ sec when memory profiling: $sample_sec"
echo "dataset name: $dataset_name"
echo "num_fixels: $num_fixels"
echo "num_subj: $num_subj"
echo "num_cores: $num_cores"
echo "run_where: $run_where"
echo "output_folder: $output_folder"


#singularity instance start -e -B $TMPDIR:/var -B $HOME:/root    /cbica/projects/fixel_db/myr_r4.1.0forFixelArray.sif myr_r4.1.0forFixelArray
#singularity instance list
#singularity exec instance://myr_r4.1.0forFixelArray Rscript /cbica/projects/fixel_db/FixelArray/notebooks/memoryProfiling_FixelArray.lm.R

if [[ "$run_where" == "vmware" ]]
then
	cmd_memrec="perl memrec"
	# cmd_memrec="perl /home/chenying/Desktop/Apps/memrec/pull_from_github/memrec/bin/memrec.pl"
	# cmd_memrec="/home/chenying/Desktop/Apps/wss-master/wss.pl"
else
	cmd_memrec="memrec"
fi

filename_memrec_output="memprofile.inMB.every${d_memrec}sec"
fn_memrec_output="${output_folder}/${filename_memrec_output}"
fn_R_output="${output_folder}/Routput.txt"
fn_myMemProf="${output_folder}/output_myMemoryProfiler.txt"

# memrec -d $d_memrec -M -o ${fn_memrec_output} 

echo ""
# if [[ "$run_where" == "sge"  ]] || [[ "$run_where" == "interactive"   ]]; then
	# cmd="${cmd_memrec} -d $d_memrec -M -o ${fn_memrec_output} singularity run --cleanenv /cbica/projects/fixel_db/myr_r4.1.0forFixelArray.sif Rscript /cbica/projects/fixel_db/FixelArray/notebooks/memoryProfiling_FixelArray.lm.R $dataset_name $num_fixels $num_subj $num_cores"
# elif [[  "$run_where" == "vmware"  ]]; then
	# which R
	# cmd="${cmd_memrec} -d $d_memrec -M -o ${fn_memrec_output} Rscript ./memoryProfiling_FixelArray.lm.R $dataset_name $num_fixels $num_subj $num_cores"  # /home/chenying/Desktop/fixel_project/FixelArray/notebooks
# fi
# echo $cmd
# $cmd

cmd="Rscript ./memoryProfiling_FixelArray.lm.R $dataset_name $num_fixels $num_subj $num_cores > ${fn_R_output}  2>&1 &"
echo $cmd
Rscript ./memoryProfiling_FixelArray.lm.R $dataset_name $num_fixels $num_subj $num_cores > ${fn_R_output}  2>&1 &     # cannot run at background if using $cmd to execuate..

parent_id=$!  # get the pid of last execuated command

echo "parent id = ${parent_id}"

#echo "sleep for 10sec to let multicore-computing starts...."
#sleep 10
#if [[ "${num_cores}" -gt 1  ]]; then     # more than one cores requested
#    child_id_list=`pgrep -P ${parent_id}`
#fi
#echo "child id(s): ${child_id_list}"

if [[ "$run_memoryProfiler" == "TRUE"  ]]; then
	bash myMemoryProfiler.sh -P ${parent_id} -c ${num_cores} -s ${sample_sec} -o ${output_folder} > ${fn_myMemProf} 2>&1
else
	echo "not to call myMemoryProfiler.sh to profile memory"

fi



#start_time=$(date +%s.%3N) # however, memrec will run in the background --> the elapsed time is not full time of running memrec

# $cmd

#end_time=$(date +%s.%3N)
#echo ""

#elapsed=$(echo "scale=3; $end_time - $start_time" | bc)   # millisec, with precision of three digits after floating point
#echo "elapsed time of memrec: ${elapsed}"

#echo ""

# /cbica/projects/fixel_db/FixelArray/notebooks/

#date
