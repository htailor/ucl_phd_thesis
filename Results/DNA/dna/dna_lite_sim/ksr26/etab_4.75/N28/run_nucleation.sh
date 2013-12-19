#!/bin/sh
INPUT_FILE=Nucleation.Input
COMPLETION_FILE=SIM_FINISHED.status

if [[ ! -f $COMPLETION_FILE ]]
then
	if [[ -z $1 ]]
	then
   		echo "Using default number of threads: 4"
		export OMP_NUM_THREADS=4
	else
   		echo "Setting number of threads: $1"
   		export OMP_NUM_THREADS=$1
	fi

	if [[ -f $INPUT_FILE ]]
	then
		echo "Found input file!!!"
   		sleep 1
   		awk < $INPUT_FILE '{print $2}' | ./Nucleation > Nucleation.Output &
	else
   	echo "No input file found."
   	exit 1
	fi
else
	echo "Results Found"
	echo "Skipping Nucleation run..."
	exit 1
fi

