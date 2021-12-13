#!/bin/bash

# example command:
# bash qsub_wrapper_benchmark_ModelArray.lm.sh -h 20G -s 1 -D test_n50 -f 100 -S 50 -c 2 -w sge

while getopts h:s:D:f:S:c:w:O:M: flag
do
        case "${flag}" in
                h) h_vmem=${OPTARG};;   # e.g. 30G
                s) sample_sec=${OPTARG};; 
                # d) d_memrec=${OPTARG};;
                D) dataset_name=${OPTARG};;
                f) num_fixels=${OPTARG};;   # 0 as full
                S) num_subj=${OPTARG};;
                c) num_cores=${OPTARG};;
                w) run_where=${OPTARG};;    # "sge" or "interactive" or "vmware"
                O) overwrite=${OPTARG};;   # "TRUE"
                M) run_memoryProfiler=${OPTARG};;   # "TRUE" or "FALSE"
        esac
done

cmd="qsub -l h_vmem=${h_vmem} -pe threaded ${num_cores} wrapper_benchmark_ModelArray.lm.sh"
cmd+=" -s $sample_sec -D $dataset_name -f $num_fixels -S $num_subj -c $num_cores -w $run_where -O ${overwrite} -M ${run_memoryProfiler}"
echo $cmd
$cmd