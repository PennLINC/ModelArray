#!/bin/bash
# This bash is to profile memory usage of a process and its child(ren)


while getopts P:c:s:o: flag
do
	case "${flag}" in
        P) parent_id=${OPTARG};;
		c) num_cores=${OPTARG};;
        s) sample_sec=${OPTARG};;
        o) output_folder=${OPTARG};;
	esac
done

#echo "parent id: ${parent_id}"
#echo "num_cores: ${num_cores}"
#echo "sample_sec: ${sample_sec}"
#echo "output folder: ${output_folder}"


source config.txt
echo $flag_where    # "cubic" or "vmware" or "dopamine"

if [[ "${flag_where}" == "vmware"   ]]; then
    cmd_wss="/home/chenying/Desktop/Apps/wss-master/wss.pl"
elif [[ "${flag_where}" == "cubic"  ]]; then
    cmd_wss="/cbica/projects/fixel_db/Apps/wss/wss.pl"
fi

# parent_id=`ps -aux | grep "[R] --no-echo --no-restore" | awk '{print $2}'`     # [R] is to remove result of grep search; "awk" is to return pid only
echo "parent id of running job = ${parent_id}"

# start to wss of parent_id, save to a file, at background
fn_parent_singlecore="${output_folder}/wss_SingleCoreStarts_parent.txt"

${cmd_wss} -C ${parent_id} ${sample_sec} > ${fn_parent_singlecore} 2>&1 &


#echo "output folder: ${output_folder}"
#echo "fn_parent_singlecore: ${fn_parent_singlecore}"

# get the id of parent's profiling
pid_wss_parent_singlecore=$!
# echo "pid of wss parent = ${pid_wss_parent_singlecore}"

# set up the fn after multiple cores start:
fn_parent_multicore="${output_folder}/wss_MultiCoreStarts_parent.txt"


fn_child_list=()
for (( i=0; i<${num_cores}; i++ ))
do
    fn_child="${output_folder}/wss_MultiCoreStarts_child${i}.txt"
    fn_child_list+=($fn_child)
    # echo "for child process # ${i}"
    
done
echo "fn for all child(ren)'s process(es): ${fn_child_list[@]}"

if [[ "${num_cores}" -gt 1  ]]; then     # more than one cores requested
    while :    # while TRUE until find all children's ids:
    do
        child_id_list_wn=`pgrep -P ${parent_id}`
        # echo "child id list with \n: ${child_id_list_wn}"
        readarray -t child_id_list <<<"$child_id_list_wn"    # remove "\n"

        
        if [  ${#child_id_list[@]} -eq ${num_cores} ]; then    # if length matches number of requested cores
            
            echo "found all child(ren)'s id(s): ${child_id_list[@]}"
            # echo "child id list - all: ${child_id_list[@]}"
            #echo "child id list[0]: ${child_id_list[0]}"
            #echo "child id list[1]: ${child_id_list[1]}"
            # echo "number of child id found = ${#child_id_list[@]}"

            break   # break while
        fi
    done

    kill -9 ${pid_wss_parent_singlecore}   # kill parent's profiling
    # we can say something at the end of ${fn_parent_singlecore}

    # then, starts profiling of both parent and children's processes, using wss
    ${cmd_wss} -C ${parent_id} ${sample_sec} > ${fn_parent_multicore} 2>&1 &
    #pid_wss_parent_multicore=$!  # actually this is not useful: when killing an R session, this wss will also terminate
    #pid_wss_child_list=()
    for (( i=0; i<${num_cores}; i++ ))
    do
        ${cmd_wss} -C ${child_id_list[$i]} ${sample_sec} > ${fn_child_list[$i]} 2>&1 &
        #pid_wss_child_list+=($!)
    done
    
    #echo "pid for wss of parent after multiple cores start: ${pid_wss_parent_multicore}"
    #echo "pid(s) for all wss of child processes: ${pid_wss_child_list[@]}"
fi