#!/bin/bash

# examples:
# bash wrapper_benchmark_fixelcfestats.sh -S 30 -h 100 -f TRUE -F TRUE -w interactive -M TRUE -s 1

while getopts S:h:f:F:w:M:s: flag
do
    case "${flag}" in
        S) num_subj=${OPTARG};;
        h) nshuffles=${OPTARG};;
        f) ftests=${OPTARG};;   # TRUE or FALSE
        F) fonly=${OPTARG};;   # TRUE or FALSE
        w) run_where=${OPTARG};;   # TRUE or FALSE
        M) run_memoryProfiler=${OPTARG};;   # "TRUE" or "FALSE"
        s) sample_sec=${OPTARG};;
    esac
done

printf -v date '%(%Y%m%d-%H%M%S)T' -1   # $date, in YYYYmmdd-HHMMSS
echo "date variable: ${date}"

folder_main_output="/cbica/projects/fixel_db/dropbox/data_from_josiane/for_fixelcfestats/stats_FDC"
folder_output="${folder_main_output}/n-${num_subj}.nshuffles-${nshuffles}"

if [[ "$ftests" == "TRUE"   ]]; then
    folder_output+=".ftests"
fi


if [[ "$fonly" == "TRUE"   ]]; then
    folder_output+=".fonly"
fi

folder_outputs+=".${run_where}"

if [[ "$run_memoryProfiler" == "TRUE"  ]]; then
        folder_output="${folder_output}.runMemProfiler"  
else 
        folder_output="${folder_output}.noMemProfiler"
fi

folder_output+=".s-${sample_sec}sec"

folder_output+=".${date}"
echo $folder_output
mkdir -p $folder_output
fn_output_txt="${folder_output}/output.txt"



bash benchmark_fixelcfestats.sh -S $num_subj -h $nshuffles -f $ftests -F $fonly -o $folder_output -M ${run_memoryProfiler} -s ${sample_sec} > $fn_output_txt 2>&1 &