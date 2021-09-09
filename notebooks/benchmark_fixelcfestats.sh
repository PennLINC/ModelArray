#!/bin/bash

while getopts S:h:t:f:F:n:o:M:s:w: flag
do
    case "${flag}" in
        S) num_subj=${OPTARG};;
        h) nshuffles=${OPTARG};;
        t) nthreads=${OPTARG};;   # number of threads in fixelcfestats
        f) ftests=${OPTARG};;   # TRUE or FALSE
        F) fonly=${OPTARG};;   # TRUE or FALSE
        n) notest=${OPTARG};;   # TRUE or FALSE
        o) folder_output=${OPTARG};;
        M) run_memoryProfiler=${OPTARG};;   # "TRUE" or "FALSE"
        s) sample_sec=${OPTARG};;
        w) run_where=${OPTARG};;  # "interactive" or "vmware"
    esac
done


if [[ "$run_where" == "interactive"  ]]; then
    folder_main_data="/cbica/projects/fixel_db/dropbox/data_from_josiane/for_fixelcfestats"

    source /cbica/projects/fixel_db/miniconda3/etc/profile.d/conda.sh   # found by $ conda info | grep -i 'base environment'
    conda activate

elif [[ "$run_where" == "vmware"  ]]; then
    folder_main_data="/home/chenying/Desktop/fixel_project/data/data_from_josiane/for_fixelcfestats"
fi


fn_myMemProf="${folder_output}/output_myMemoryProfiler.txt"

cmd="fixelcfestats ${folder_main_data}/fdc_10_smoothed_10fwhm_new/"
cmd+=" ${folder_main_data}/list_filenames_n${num_subj}.txt"
cmd+=" ${folder_main_data}/designMatrix_lm_Age_n${num_subj}.txt"
cmd+=" ${folder_main_data}/contrastMatrix_lm_n${num_subj}.txt"
cmd+=" ${folder_main_data}/../MRtrix_conn_matrix/matrix/"
cmd+=" $folder_output"
cmd+=" -nthreads $nthreads"
cmd+=" -nshuffles $nshuffles"

if [[ "$ftests" == "TRUE"   ]]; then
    cmd+=" -ftests ${folder_main_data}/Ftest_lm_n${num_subj}.txt"
fi


if [[ "$fonly" == "TRUE"   ]]; then
    cmd+=" -fonly"
fi

if [[ "$notest" == "TRUE"   ]]; then
    cmd+=" -notest"
fi

echo $cmd
$cmd &

parent_id=$!  # get the pid of last execuated command

echo "parent id = ${parent_id}"

if [[ "$run_memoryProfiler" == "TRUE"  ]]; then
	bash myMemoryProfiler.sh -P ${parent_id} -c 1 -s ${sample_sec} -o ${folder_output} > ${fn_myMemProf} 2>&1
else
	echo "not to call myMemoryProfiler.sh to profile memory"

fi