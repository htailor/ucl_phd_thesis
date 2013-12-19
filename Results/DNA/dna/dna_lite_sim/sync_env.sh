#!/bin/sh

MASTER_KSR_ENV="../MASTER_KSR"
MASTER_ETA_ENV="../../MASTER_KSR/MASTER_ETA"
MASTER_DIR_SOURCE="../../MASTER_KSR/MASTER_ETA/MASTER_N"

N=(N12 N16 N20 N24 N28 N32 N50)
ETA_B=(etab_5.00)

backup_scripts(){

	count=`ls -1 *.sh 2>/dev/null | wc -l`
	if [ $count != 0 ]
	then
		rm -f *.old
		for FILE in *.sh
		do 
			mv "${FILE}" "${FILE%.sh}.sh.old"
			chmod -x "${FILE%.sh}.sh.old"
		done
	fi 
}


KSR_DIR=`ls -F | grep ksr | grep / | sed 's/\///' | sort -n`


for KSR in $KSR_DIR
do
	echo "Syncing $KSR..."
	cd $KSR

	#Create ETA_B enviroments
	echo ${ETA_B[*]}
	mkdir -p `echo ${ETA_B[*]}`		

	#Sync etab enviroments
	echo "Finding etab enviroments..."
	ETAB_DIR_LIST=`ls -F | grep etab_ | grep / | sed 's/\///'`
    echo "Found `ls -F | grep etab_ | grep / | sed 's/\///' | wc -l`"

	for ETAB in $ETAB_DIR_LIST
	do
		echo "Syncing $ETAB..."
		cd $ETAB

		#Create N enviroments
		mkdir -p `echo ${N[*]}`

		N_DIR_LIST=`ls -F | grep '^N[0-9]' | sed 's/\///' |  sort -n`

		for DIR in $N_DIR_LIST
		do
		        cd $DIR
#        		echo "($DIR)"
        		rm -rf ./src
        		cp -r ../$MASTER_DIR_SOURCE/logs .
        		cp -r ../$MASTER_DIR_SOURCE/results .
        		cp ../$MASTER_DIR_SOURCE/Nucleation .
        		cp ../$MASTER_DIR_SOURCE/*.sh .
        		rm -f COMPILE.sh
        		cd ..
		done

		backup_scripts
		cp $MASTER_ETA_ENV/*.sh .
		cd ..
	done

	#Sync remaining scripts in ksr enviroment
	echo "Backing up scripts in $KSR enviroment"	
	backup_scripts
	cp $MASTER_KSR_ENV/*.sh .
	./set_input_parameters.sh
	cd ..
done

echo "COMPLETE"


