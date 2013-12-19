#!/bin/sh
clear
rm -f ETA*
find . -name "*sim.data" -exec rm -rf {} \;

N_LIST=`ls -F | grep '^N[0-9]' | sed 's/\///' | sed -e 's/.//' | sort -n`
RESULTS_DIR="results/Intact"
MATH_SCRIPTS_DIR="/Scripts_Mathematica"
MATHEMATICA_SCRIPT_FILE="Intact_Analysis_large_etab.math"
MATHEMATICA_FILE="${PWD}$MATH_SCRIPTS_DIR/$MATHEMATICA_SCRIPT_FILE"
ETAB="$1"

echo -e "Running analysis for:\nETA_B = $ETAB"
echo -e "N list: `echo $N_LIST`"
echo -e "Mathematica script: $MATHEMATICA_SCRIPT_FILE"

sed -i -e 's%Subscript\[\\\[Eta\]\,\sB\]=.*%Subscript\[\\\[Eta\]\, B\]='${ETAB}';%' $MATHEMATICA_FILE

for N in $N_LIST
do
	DIR="N$N"
	RESULTS_LOCATION="$DIR/$RESULTS_DIR"
	echo -e "\nAnalysing ($DIR)\n================"
	echo "Results location: $RESULTS_LOCATION"
	echo "Running..."
	cd $RESULTS_LOCATION
	math < $MATHEMATICA_FILE > /dev/null
	echo "COMPLETED"
	cd ../../..
done

./grab_force_data.sh
echo -e "\nVerify data (ETA_B = $ETAB)...\n"
touch ETAB_$ETAB
cat force.data

ETAB_DIR="Analysis_etab_$ETAB"

mkdir -p $ETAB_DIR
cp KAPPA_SIGMA_R* $ETAB_DIR
cp ETAB_$ETAB $ETAB_DIR
cp force.data $ETAB_DIR
cp shear_force*.nb $ETAB_DIR
cp -r Hatch_Data ./$ETAB_DIR
