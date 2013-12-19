#!/bin/sh

rm -rf BATCH*
KSR_RUN_LIST=`ls -F | grep / | grep ksr | sed -e 's/\///'`

touch BATCH_JOB_STARTED
echo -e "\nSTARTING ANALYSIS..."

OPTION=$1

for KSR in $KSR_RUN_LIST
do
	cd $KSR
	echo -e "\n($KSR)"
	./run_analysis.sh $OPTION
	cd ..
done

echo -e "\nANALYSIS COMPLETE\n"
touch BATCH_JOB_FINISHED
