#!/bin/sh

ANALYSIS_DIR="Analysis"
RESULTS_DIR="results"
MATH_SCRIPTS_DIR="../../Scripts_Mathematica"
MATHEMATICA_SCRIPT_FILE="DNA_Analysis.math"
MATHEMATICA_FILE="${PWD}/$MATH_SCRIPTS_DIR/$MATHEMATICA_SCRIPT_FILE"
FORCE_DATA_FILE=force.data
N_LIST=`ls -F | grep '^N[0-9]' | sed 's/\///' | sed -e 's/.//' | sort -n`


usage(){

    echo "\nUSAGE: $0 [OPTIONS]"
    echo "OPTIONS:"
	echo "          -o	Override old analysis results\n"

}

create_force_data(){

    rm -f $FORCE_DATA_FILE
    awk 'FNR==1{print ""}1' `find . -name "*sim.data" | sort -n` > test
    sed '/^$/d' test > $FORCE_DATA_FILE
    rm test

    echo -e "\nVerify data...\n"
    cat force.data

}

create_analysis_directory(){

    mkdir -p $ANALYSIS_DIR
    cp KAPPA_SIGMA_R* $ANALYSIS_DIR
    cp ETA_B_* $ANALYSIS_DIR
    cp force.data $ANALYSIS_DIR
    cp -r ../../Hatch_Data ./$ANALYSIS_DIR
    cp -r ../../MASTER_KSR/MASTER_ETA/MASTER_ANALYSIS/*.nb ./$ANALYSIS_DIR

}

delete_old_analysis_data(){

    rm -rf $ANALYSIS_DIR
    find . -name "*sim.data" -exec rm -rf {} \;

}


##########################################################################
#Main
##########################################################################



if [[ $1 == "help" ]]; then
    usage
    exit 1
fi

while getopts "o" OPTARG
do
	case $OPTARG in
		o)
			delete_old_analysis_data
			;;
	esac
done


if [ -d "Analysis" ]; then
    echo "Analysis already run!"
    exit 1
fi


echo -e "N list: `echo $N_LIST`"
echo -e "Mathematica script: $MATHEMATICA_SCRIPT_FILE"

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
	cd ../../
done

create_force_data
create_analysis_directory

