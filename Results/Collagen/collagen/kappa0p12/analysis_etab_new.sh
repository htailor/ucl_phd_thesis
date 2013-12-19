#!/bin/sh
clear
rm -f ETA*
find . -name "*sim.data" -exec rm -rf {} \;

N_LIST=`ls -F | grep '^N[0-9]' | sed 's/\///' | sed -e 's/.//' | sort -n`

MATHEMATICA_SCRIPT_FILE="./mathematica_scripts/shear_force_new.script"
MATHEMATICA_SCRIPT_FILE_RUN="./mathematica_scripts/shear_force_script.run"
ETAB="$1"

echo -e "\nRunning analysis for:\nETA_B = $ETAB"
echo -e "N list: `echo $N_LIST`"
echo -e "Mathematica script: $MATHEMATICA_SCRIPT_FILE"

sed 's%etab=.*%etab='${ETAB}';%' $MATHEMATICA_SCRIPT_FILE >$MATHEMATICA_SCRIPT_FILE_RUN

for N in $N_LIST
do
      DIR="N$N"
      cd $DIR
      echo -e "\nAnalysing ($DIR)\n================"
      echo "Running..."
      math < ../$MATHEMATICA_SCRIPT_FILE_RUN > /dev/null
      echo "COMPLETED"
      cd ..
done

./grab_force_data.sh
echo -e "\nVerify data (ETA_B = $ETAB)...\n"
touch ETAB_$ETAB
cat force.data
mv force.data new_force.data
