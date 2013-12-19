#!/bin/bash

INPUT_FILE="Nucleation.Input"
N_LIST=`ls -F | grep '^N[0-9]' | sed 's/\///' | sed -e 's/.//' | sort -n`
echo "N list: `echo $N_LIST`"

L=60.2
M=150
KAPPA=0.1
SIGMA=0.000666
UMIN=0
UMAX=70

KAPPA_SIGMA_R=`echo $KAPPA / $SIGMA | bc -l`
KSR=`echo $KAPPA_SIGMA_R | awk '{printf "%.2f",$1}'`

touch "KAPPA_SIGMA_R_$KSR"

for N in $N_LIST
do
	DIR="N$N"
	echo "Creating input file ($DIR)"
	cd $DIR
	rm -f $INPUT_FILE
	echo "L: $L" > $INPUT_FILE
	echo "m: $M" >> $INPUT_FILE
	echo "N: $N" >> $INPUT_FILE

	echo "kappa: $KAPPA" >> $INPUT_FILE
	echo "sigma: $SIGMA" >> $INPUT_FILE
	echo "umin: $UMIN" >> $INPUT_FILE
	echo "umax: $UMAX" >> $INPUT_FILE
	cd ..
done


