#!/bin/sh

export OMP_NUM_THREADS=8
export EIGEN_LIMIT_SMAX=5
export EIGEN_LIMIT_TMAX=5

BINARY_FILE="./Release/collagen"

rm -f KAPPA_SIGMA_R*
./clear_all.sh
find . -type d -name "N*" -exec rm -rf {} \; 2>/dev/null

L=80.1
M=400

KAPPA=1.2
SIGMA=0.00024
KAPPA_SIGMA_R=$(printf %.0f `echo "$KAPPA/$SIGMA" | bc -l`)

UMIN=0
UMAX=150

touch KAPPA_SIGMA_R_$KAPPA_SIGMA_R

$BINARY_FILE -L $L -m $M --kappa $KAPPA --sigma $SIGMA --umin $UMIN --umax $UMAX > output &
