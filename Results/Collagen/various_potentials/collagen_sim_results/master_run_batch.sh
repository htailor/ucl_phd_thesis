#!/bin/bash

L=6
M=24
N=5

KAPPA=1
SIGMA=1

UMIN=0
UMAX=20

DIR_LIST=`ls -l | grep "drw" | awk '{print $9}'`

for DIR in $DIR_LIST
do
	echo "Running..$DIR";
	cd $DIR
        	./run_batch.sh $L $M $N $KAPPA $SIGMA $UMIN $UMAX &
        cd ..
done
