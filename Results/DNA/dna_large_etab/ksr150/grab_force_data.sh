#!/bin/sh

FORCE_DATA_FILE=force.data

rm -f $FORCE_DATA_FILE
awk 'FNR==1{print ""}1' `find . -name "*sim.data" | sort -n` > test
sed '/^$/d' test > $FORCE_DATA_FILE
rm test
