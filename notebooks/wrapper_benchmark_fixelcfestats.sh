#!/bin/bash

# examples:
# bash wrapper_benchmark_fixelcfestats.sh -S 938 -h 100 -t 2 -f TRUE -F FALSE -n FALSE -w vmware -M TRUE -s 0.1
# bash wrapper_benchmark_fixelcfestats.sh -S 938 -h 100 -t 1 -f TRUE -F FALSE -n TRUE -w interactive -M TRUE -s 1

while getopts S:h:t:f:F:n:w:M:s: flag
do
    case "${flag}" in
        S) num_subj=${OPTARG};;
        h) nshuffles=${OPTARG};;
        t) nthreads=${OPTARG};;   # number of threads in fixelcfestats
        f) ftests=${OPTARG};;   # TRUE or FALSE
        F) fonly=${OPTARG};;   # TRUE or FALSE
        n) notest=${OPTARG};;   # TRUE or FALSE
        w) run_where=${OPTARG};;   # "interactive" or "vmware" or "sge"
        M) run_memoryProfiler=${OPTARG};;   # "TRUE" or "FALSE"
        s) sample_sec=${OPTARG};;
    esac
done

echo "JOB_ID = ${JOB_ID}"

printf -v date '%(%Y%m%d-%H%M%S)T' -1   # $date, in YYYYmmdd-HHMMSS
echo "date variable: ${date}"

if [[ "$run_where" == "interactive"  ]] || [[ "$run_where" == "sge"  ]]; then
    folder_main_output="/cbica/projects/fixel_db/dropbox/data_from_josiane/for_fixelcfestats/stats_FDC"
elif [[ "$run_where" == "vmware"  ]]; then
    folder_main_output="/home/chenying/Desktop/fixel_project/data/data_from_josiane/for_fixelcfestats/stats_FDC"
fi

folder_output="${folder_main_output}/nsubj-${num_subj}.nthreads-${nthreads}"

if [[ "$ftests" == "TRUE"   ]]; then
    folder_output+=".ftests"
fi


if [[ "$fonly" == "TRUE"   ]]; then
    folder_output+=".fonly"
fi

if [[ "$notest" == "TRUE"   ]]; then
    folder_output+=".notest"
fi

folder_output+=".nshuffles-${nshuffles}"

folder_output+=".${run_where}"

if [[ "$run_memoryProfiler" == "TRUE"  ]]; then
        folder_output="${folder_output}.runMemProfiler"  
else 
        folder_output="${folder_output}.noMemProfiler"
fi

folder_output+=".s-${sample_sec}sec"

folder_output+=".${date}"

if [[ "$run_where" == "sge"  ]]; then
    folder_output+=".${JOB_ID}"
fi

echo $folder_output
mkdir -p $folder_output
fn_output_txt="${folder_output}/output.txt"



bash benchmark_fixelcfestats.sh -S $num_subj -h $nshuffles -t $nthreads -f $ftests -F $fonly -n $notest -o $folder_output -M ${run_memoryProfiler} -s ${sample_sec} -w ${run_where} > $fn_output_txt 2>&1 &