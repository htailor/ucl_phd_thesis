#!/bin/sh

ETAB_DIR_LIST=`ls -F | grep etab_ | grep / | sed 's/\///' | sort -n`

for ETAB_DIR in $ETAB_DIR_LIST
do
    cd $ETAB_DIR
    echo -e "\n$ETAB_DIR..."
    find . -name "SIM_FINISHED.status" | sort -n
    cd ..
done

echo "COMPLETE"
