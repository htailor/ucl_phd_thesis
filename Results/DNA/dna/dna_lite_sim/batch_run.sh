#!/bin/sh

rm -rf BATCH*
KSR_LIST=`ls -F | grep / | grep ksr | sed -e 's/\///'`

touch BATCH_JOB_STARTED

for KSR in $KSR_LIST
do
	cd $KSR
	echo "(KSR:$KSR)"
	./run.sh 8
	cd ..
done

touch BATCH_JOB_FINISHED
