#!/bin/bash

# example command:
# bash qsub_wrapper_benchmark_FixelArray.lm.sh -h 20G -d 0.01 -D test_n50 -f 100 -s 50 -c 2 -w sge   # tested, it works!

while getopts h:d:D:f:s:c:w:O: flag
do
        case "${flag}" in
                h) h_vmem=${OPTARG};;   # e.g. 30G
                d) d_memrec=${OPTARG};;
                D) dataset_name=${OPTARG};;
                f) num_fixels=${OPTARG};;   # 0 as full
                s) num_subj=${OPTARG};;
                c) num_cores=${OPTARG};;
                w) run_where=${OPTARG};;    # "sge" or "interactive" or "vmware"
                O) overwrite=${OPTARG};;   # "TRUE"
        esac
done

cmd="qsub -l h_vmem=${h_vmem} -pe threaded ${num_cores} wrapper_benchmark_FixelArray.lm.sh"
cmd+=" -d $d_memrec -D $dataset_name -f $num_fixels -s $num_subj -c $num_cores -w $run_where -O ${overwrite}"
echo $cmd
$cmd