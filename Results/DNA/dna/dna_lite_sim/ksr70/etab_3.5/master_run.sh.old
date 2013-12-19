#!/bin/sh

NUM_THREADS=$1
DIR_LIST=`ls -F | grep '^N[0-9]' | sed 's/\///' | sort -n`

for DIR in $DIR_LIST
do
	echo "Running..$DIR";
	cd $DIR
	./run_nucleation.sh $NUM_THREADS
	cd ..
done
