#!/bin/sh
clear
N_LIST=`ls -F | grep '^N[0-9]' | sed 's/\///' | sed -e 's/.//' | sort -n`
echo "N list: `echo $N_LIST`"

for N in $N_LIST
do
	DIR="N$N"
	echo -e "\nSimulation status ($DIR)\n======================="
	cd $DIR
	SIM_START_DATE_TIME=`ls -altr | grep SIM_STARTED.status | awk '{print $6" "$7" "$8}'`
	SIM_FINISHED_DATE_TIME=`ls -altr | grep SIM_FINISHED.status | awk '{print $6" "$7" "$8}'`

	echo -e "start time:\t $SIM_START_DATE_TIME"
	echo -e "end time:\t $SIM_FINISHED_DATE_TIME"
	cd ..
done


