#!/bin/bash

# example command:
# bash wrapper_benchmark_FixelArray.lm.sh -d 0.01 -D test_n50 -f 100 -s 100 -c 2 -w interactive -O TRUE
# qsub -l h_vmem=30G wrapper_benchmark_FixelArray.lm.sh -d 0.01 -D test_n50 -f 100 -s 100 -c 2 -w sge    # this will add ${JOB_ID} to foldername; run at interactive node to determine the memory requirements... # tried 20G, did not run..

while getopts d:D:f:s:c:w:O: flag
do
        case "${flag}" in
                d) d_memrec=${OPTARG};;
                D) dataset_name=${OPTARG};;
                f) num_fixels=${OPTARG};;   # 0 as full
                s) num_subj=${OPTARG};;
                c) num_cores=${OPTARG};;
                w) run_where=${OPTARG};;    # "sge" or "interactive" or "vmware"
                O) overwrite=${OPTARG};;   # "TRUE"
        esac
done

echo "JOB_ID = ${JOB_ID}"

printf -v date '%(%Y%m%d-%H%M%S)T' -1   # $date, in YYYYmmdd-HHMMSS
echo "date variable: ${date}"

foldername_jobid="lm.${dataset_name}.nfixel-${num_fixels}.nsubj-${num_subj}.ncore-${num_cores}.${run_where}"

if [  "$run_where" = "sge" ]; then
        folder_benchmark="/cbica/projects/fixel_db/FixelArray_benchmark"
        
        echo "adding JOB_ID to foldername"
        foldername_jobid="${foldername_jobid}.${JOB_ID}"

elif [[ "$run_where" == "interactive"   ]]; then
        folder_benchmark="/cbica/projects/fixel_db/FixelArray_benchmark"

elif [[ "$run_where" == "vmware"   ]]; then
        folder_benchmark="/home/chenying/Desktop/fixel_project/FixelArray_benchmark"

        echo "adding date to foldername"
        foldername_jobid="${foldername_jobid}.${date}"
fi

folder_jobid="${folder_benchmark}/${foldername_jobid}"
# echo "folder_jobid: ${folder_jobid}"


if [ -d ${folder_jobid} ] && [ "${overwrite}" = "TRUE" ]
then
        echo "removing existing folder:   ${folder_jobid}"
        rm -r ${folder_jobid}
fi
mkdir ${folder_jobid}
echo "output folder:   ${folder_jobid}"
# echo "output foldername for this job: foldername_jobid"

fn_output_txt="${folder_jobid}/output.txt"
# echo "fn_output_txt: ${fn_output_txt}"

# call:
bash benchmark_FixelArray.lm.sh -d $d_memrec -D $dataset_name -f $num_fixels -s $num_subj -c $num_cores -w $run_where -o ${folder_jobid} > $fn_output_txt 2>&1

