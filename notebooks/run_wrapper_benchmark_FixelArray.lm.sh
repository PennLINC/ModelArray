#!/bin/bash

function wait_r_end {
    while :
    do
        pid_r=`ps -aux | grep "[R] --no-echo --no-restore" | awk '{print $2}'`     # [R] is to remove result of grep search; "awk" is to return pid only
        if [ -z "$pid_r" ]; then
            date
            echo "sleep for $1 sec..."
            sleep $1
            break
        fi
    done
}

# use for correction: # checkout parent.multi and child0.multi
#bash wrapper_benchmark_FixelArray.lm.sh -s 1 -D test_n50 -f 1000 -S 50 -c 4 -w vmware -M TRUE


wait_r_end 10   # in seconds
date
cmd="bash wrapper_benchmark_FixelArray.lm.sh -s 1 -D josiane -f 0 -S 938 -c 1 -w vmware -M TRUE"
echo $cmd
$cmd
date

wait_r_end 600   # in seconds
date
cmd="bash wrapper_benchmark_FixelArray.lm.sh -s 1 -D josiane -f 0 -S 938 -c 2 -w vmware -M TRUE"
echo $cmd
$cmd
date

