#!/bin/ksh

BINARY_LOCATION="Release"
BINARY_FILE="collagen"

L=5
M=20
N=5

KAPPA=1
SIGMA=1

UMIN=0
UMAX=20

RUN_COMMAND="./$BINARY_LOCATION/$BINARY_FILE -L $L -m $M -N $N --kappa $KAPPA --sigma $SIGMA --umin $UMIN --umax $UMAX"

`$RUN_COMMAND &`

exit

