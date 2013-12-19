#!/bin/bash

LOG_FILE="output.log"

if [ -e $LOG_FILE ]
then
	rm $LOG_FILE
fi

COMMAND="./run_simulation_batch.sh $1 $2 $3 $4 $5 $6 $7"

/bin/bash $COMMAND > output

mv output output.log

