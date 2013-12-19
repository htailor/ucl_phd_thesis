#!/bin/bash

if [ $# -eq 7 ]; 
then
        L=$1
	M=$2
	ETA_B=$3
	KAPPA=$4
	SIGMA=$5
	UMIN=$6
	UMAX=$7
else
        echo "Incorrect number of arguments passed!"
	echo  "USAGE: $0 [L] [M] [ETA_B] [KAPPA] [SIGMA] [UMIN] [UMAX]"
        exit 1             
fi

rm -f ETA_B_*
rm -f KAPPA_SIGMA_R_*
rm -f KAPPA_*
rm -f SIGMA_*

INPUT_FILE="Nucleation.Input"
N_LIST=`ls -F | grep '^N[0-9]' | sed 's/\///' | sed -e 's/.//' | sort -n`
echo "N list: `echo $N_LIST`"



KAPPA_SIGMA_R=`echo $KAPPA / $SIGMA | bc -l`
KSR=`echo $KAPPA_SIGMA_R | awk '{printf "%.2f",$1}'`

touch "ETA_B_$ETA_B"
touch "KAPPA_$KAPPA"
touch "SIGMA_$SIGMA"
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
	echo "eta_b: $ETA_B" >> $INPUT_FILE
	echo "kappa: $KAPPA" >> $INPUT_FILE
	echo "sigma: $SIGMA" >> $INPUT_FILE
	echo "umin: $UMIN" >> $INPUT_FILE
	echo "umax: $UMAX" >> $INPUT_FILE
	cd ..
done


