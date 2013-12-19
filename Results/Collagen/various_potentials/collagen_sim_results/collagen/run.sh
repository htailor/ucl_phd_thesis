#!/bin/ksh

BINARY_LOCATION="Release"
BINARY_FILE="collagen"

L=10.25
M=20
N=1

KAPPA=0.1
SIGMA=0.001

UMIN=0
UMAX=10

RUN_COMMAND="./$BINARY_LOCATION/$BINARY_FILE -L $L -m $M -N $N --kappa $KAPPA --sigma $SIGMA --umin $UMIN --umax $UMAX"

`$RUN_COMMAND &`

exit

