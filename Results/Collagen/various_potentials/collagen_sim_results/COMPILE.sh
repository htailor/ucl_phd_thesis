#!/bin/bash

DIR_LIST=`ls -l | grep "drw" | awk '{print $9}'`

for DIR in $DIR_LIST
do
	echo "Running..$DIR";
	cd $DIR
		cd Release
		make clean
		make -f makefile
		cd ..
        cd ..
done
