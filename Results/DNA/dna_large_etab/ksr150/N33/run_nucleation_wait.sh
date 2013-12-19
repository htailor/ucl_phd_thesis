#!/bin/sh
INPUT_FILE=Nucleation.Input

if [[ -f core* ]] 
then
   rm core.*
fi
if [[ -f $INPUT_FILE ]]
then
   echo "Found input file!!!"
   sleep 1
   awk < $INPUT_FILE '{print $2}' | ./Nucleation > Nucleation.Output
else
   echo "No input file found. Initialising front end GUI."
   sleep 1
   clear
   ./Nucleation
fi
