#!/bin/sh

echo "Creating nucleation input files..."

rm -rf KAPPA*
rm -rf SIGMA*

L=100.25
M=200
KAPPA=0.1
KAPPA_SIGMA_R=${PWD##/*ksr}
SIGMA=`echo $KAPPA / $KAPPA_SIGMA_R | bc -l | awk '{printf "%.8f",$1}'`
UMIN=0
UMAX=60

ETAB_LIST=`ls -F | grep etab_ | grep / | sed 's/etab_//' | sed 's/\///'`
ETAB_DIR_LIST=`ls -F | grep etab_ | grep / | sed 's/\///'`
echo "ETAB list:  " $ETAB_LIST
echo "DIR  list:  " $ETAB_DIR_LIST

for ETAB_DIR in $ETAB_DIR_LIST
do
	ETAB=`echo $ETAB_DIR | sed 's/etab_//'`
	cd $ETAB_DIR
	pwd
	./create_input_files.sh $L $M $ETAB $KAPPA $SIGMA $UMIN $UMAX
	cd ..
done

touch "KAPPA_$KAPPA"
touch "SIGMA_$SIGMA"

echo "COMPLETE"
