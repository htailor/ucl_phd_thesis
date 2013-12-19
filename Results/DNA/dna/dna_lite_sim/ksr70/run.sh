#!/bin/sh

#echo "Clearing status files..."
#find . -name "*.status" -exec rm -rf {} \;
#echo "COMPLETE."

NUM_THREADS=$1

RUN_FILE=`find . -name "master_run_wait.sh" | sort -n`

for FILE in $RUN_FILE
do
	DIR_PATH=${FILE%/*}
	echo $DIR_PATH
	cd $DIR_PATH
	./master_run_wait.sh $NUM_THREADS
	echo "COMPLETE"
	cd ..
done
