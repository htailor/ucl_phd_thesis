#!/bin/sh

MACHINE=manicorensis
MATHEMATICA_MACHINE=hendrix
NUM_THREADS=8

#KSR_RUN_LIST=(ksr24 ksr26 ksr28 ksr30 ksr32 ksr34 ksr36 ksr38 ksr40 ksr42 ksr44 ksr46 ksr48 ksr50 ksr55 ksr60 ksr65 ksr70 ksr75)
KSR_RUN_LIST=(ksr24 ksr26 ksr28 ksr30)

echo -e "\nSTARTING SIMULATIONS..."
echo "MACHINE: $MACHINE"
for KSR in ${KSR_RUN_LIST[*]}
do
	cd $KSR
	echo -e "\n($KSR)"
	ssh $MACHINE "hostname;cd $PWD;./run.sh $NUM_THREADS"
	ssh $MATHEMATICA_MACHINE "hostname;cd $PWD;./run_analysis.sh -o"
	cd ..
done

echo -e "\nSIMULATIONS COMPLETE\n"
