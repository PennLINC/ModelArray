#!/bin/bash
temp=`ldconfig -p | grep libhdf5* | wc -l`
if [[ $temp -gt 0 ]]
then
	echo "libhdf5 exists."
else
	sudo apt-get update -y
	sudo apt-get install -y libhdf5-dev
fi
