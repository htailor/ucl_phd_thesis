#!/bin/sh

FORCE_DATA_FILE=force.data

rm -f $FORCE_DATA_FILE
awk 'FNR==1{print ""}1' `find . -name "*sim.data" | sort -n` > dump1
sed '/^$/d' dump1 > dump2 
cat dump2 | sort -n > $FORCE_DATA_FILE
rm dump1 dump2
