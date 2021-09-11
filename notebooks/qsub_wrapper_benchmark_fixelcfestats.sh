#!/bin/bash

# example:
# bash qsub_wrapper_benchmark_fixelcfestats.sh -v 30G -S 30 -h 100 -t 1 -f TRUE -F FALSE -n TRUE -M TRUE -s 1

while getopts v:S:h:t:f:F:n:M:s: flag
do
    case "${flag}" in
        v) h_vmem=${OPTARG};;   # e.g. 30G
        S) num_subj=${OPTARG};;
        h) nshuffles=${OPTARG};;
        t) nthreads=${OPTARG};;   # number of threads in fixelcfestats
        f) ftests=${OPTARG};;   # TRUE or FALSE
        F) fonly=${OPTARG};;   # TRUE or FALSE
        n) notest=${OPTARG};;   # TRUE or FALSE
        M) run_memoryProfiler=${OPTARG};;   # "TRUE" or "FALSE"
        s) sample_sec=${OPTARG};;
    esac
done

cmd="qsub -l h_vmem=${h_vmem} wrapper_benchmark_fixelcfestats.sh"  # -pe threaded 1 
cmd+=" -S $num_subj -h $nshuffles -t $nthreads -f $ftests -F $fonly -n $notest"
cmd+=" -w sge -M $run_memoryProfiler -s $sample_sec"
echo $cmd
$cmd