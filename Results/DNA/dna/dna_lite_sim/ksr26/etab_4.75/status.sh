#!/bin/sh
clear
N_LIST=`ls -F | grep '^N[0-9]' | sed 's/\///' | sed -e 's/.//' | sort -n`
#N_LIST=`ls -F | grep 'N49' | sed 's/\///' | sed -e 's/.//' | sort -n`
echo "N list: `echo $N_LIST`"

TIME=0

for N in $N_LIST
do
	DIR="N$N"
	echo -e "\nSimulation status ($DIR)\n======================="
	cd $DIR
	SIM_START_DATE_TIME=`ls -altr | grep SIM_STARTED.status | awk '{print $6" "$7" "$8}'`
	SIM_FINISHED_DATE_TIME=`ls -altr | grep SIM_FINISHED.status | awk '{print $6" "$7" "$8}'`

	echo -e "start time:\t $SIM_START_DATE_TIME"
	echo -e "end time:\t $SIM_FINISHED_DATE_TIME"

	D1=`date +%s -d "$SIM_FINISHED_DATE_TIME"`
	D2=`date +%s -d "$SIM_START_DATE_TIME"`

	if [ ! -z "$SIM_FINISHED_DATE_TIME" ];
	then
	((diff_sec=D1-D2))
	((TIME=TIME+diff_sec))
	SIM_RUNTIME=`echo | awk '{printf "%02d:%02d","'"$diff_sec"'"/(60*60),"'"$diff_sec"'"%(60*60)/60}'`
	echo -e "\nruntime:\t $SIM_RUNTIME\n"
	fi
	cd ..
done

TOTAL_SIM_RUNTIME=`echo | awk '{printf "%02d:%02d","'"$TIME"'"/(60*60),"'"$TIME"'"%(60*60)/60}'`
echo -e "\nTOTAL SIMULATION RUNTIME (HH:MM) : $TOTAL_SIM_RUNTIME\n"
