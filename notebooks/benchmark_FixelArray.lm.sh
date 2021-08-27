#!/bin/bash
while getopts d:D:f:s:c:w:o: flag
do
	case "${flag}" in
		d) d_memrec=${OPTARG};;
		D) dataset_name=${OPTARG};;
		f) num_fixels=${OPTARG};;   # 0 as full
		s) num_subj=${OPTARG};;
		c) num_cores=${OPTARG};;
		w) run_where=${OPTARG};;    # "sge" or "interactive" or "vmware"
		o) output_folder=${OPTARG};;  # generated in wrapper.sh
	esac
done

date

echo "sampling every __ sec in memrec: $d_memrec"
echo "dataset name: $dataset_name"
echo "num_fixels: $num_fixels"
echo "num_subj: $num_subj"
echo "num_cores: $num_cores"
echo "run_where: $run_where"
echo "output_folder: $output_folder"


#singularity instance start -e -B $TMPDIR:/var -B $HOME:/root    /cbica/projects/fixel_db/myr_r4.1.0forFixelArray.sif myr_r4.1.0forFixelArray
#singularity instance list
#singularity exec instance://myr_r4.1.0forFixelArray Rscript /cbica/projects/fixel_db/FixelArray/notebooks/memoryProfiling_FixelArray.lm.R

filename_memrec_output="memprofile.inMB.every${d_memrec}sec"
fn_memrec_output="${output_folder}/${filename_memrec_output}"

# memrec -d $d_memrec -M -o ${fn_memrec_output} 

echo ""
cmd="memrec -d $d_memrec -M -o ${fn_memrec_output} singularity run --cleanenv /cbica/projects/fixel_db/myr_r4.1.0forFixelArray.sif Rscript /cbica/projects/fixel_db/FixelArray/notebooks/memoryProfiling_FixelArray.lm.R $dataset_name $num_fixels $num_subj $num_cores"
echo $cmd

#start_time=$(date +%s.%3N) # however, memrec will run in the background --> the elapsed time is not full time of running memrec

$cmd

#end_time=$(date +%s.%3N)
#echo ""

#elapsed=$(echo "scale=3; $end_time - $start_time" | bc)   # millisec, with precision of three digits after floating point
#echo "elapsed time of memrec: ${elapsed}"

#echo ""

# /cbica/projects/fixel_db/FixelArray/notebooks/

#date
