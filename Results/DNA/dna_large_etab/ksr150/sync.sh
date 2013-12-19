#!/bin/sh
clear
MASTER_DIR_SOURCE="MASTER"
INT_MATHEMATICA_FILE="results/Intact/Intact_Analysis.nb"
N_DIR_LIST=`ls -F | grep '^N[0-9]' | sed 's/\///' |  sort -n`

echo "Syncing from $MASTER_DIR_SOURCE..."
echo "N list: `echo $N_DIR_LIST`"

echo -e "\nCompiling master source..."
cd $MASTER_DIR_SOURCE
./COMPILE.sh > /dev/null
cd ..
echo -e "COMPLETE\n"

for DIR in $N_DIR_LIST
do
	cd $DIR
	echo -e "($DIR)"
	rm -r ./src
	cp -r ../$MASTER_DIR_SOURCE/src .
	cp ../$MASTER_DIR_SOURCE/Nucleation .
	cp ../$MASTER_DIR_SOURCE/$INT_MATHEMATICA_FILE ./$INT_MATHEMATICA_FILE
	cd ..
done

echo -e "\nSyncing source files COMPLETE"
echo -e "Syncing Nucleation binary COMPLETE"
echo -e "Syncing Intact_State.nb COMPLETE\n"

