#!/bin/sh

echo "Running analysis for all N enviroments..."

ETAB_DIR_LIST=`ls -F | grep etab_ | grep / | sed 's/\///'`
echo "DIR  list:  " $ETAB_DIR_LIST

OPTION=$1

for ETAB_DIR in $ETAB_DIR_LIST
do
	cd $ETAB_DIR
	./analysis.sh $OPTION
	cd ..
done

echo "COMPLETE"
