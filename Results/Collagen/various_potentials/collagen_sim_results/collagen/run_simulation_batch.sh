#!/bin/bash

BINARY_LOCATION="Release"
BINARY_FILE=`ls $BINARY_LOCATION -l | grep "\-rwx" | awk '{print $9}'`

L=$1
M=$2
N=$3

KAPPA=$4
SIGMA=$5

UMIN=$6
UMAX=$7

RUN_COMMAND="./$BINARY_LOCATION/$BINARY_FILE -L $L -m $M -N $N --kappa $KAPPA --sigma $SIGMA --umin $UMIN --umax $UMAX"

$RUN_COMMAND

exit
