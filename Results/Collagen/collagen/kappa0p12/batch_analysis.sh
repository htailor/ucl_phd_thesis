#!/bin/sh

ETAB=$1
ANALYSIS_DIR="Analysis_$1"

rm -rf $ANALYSIS_DIR
mkdir $ANALYSIS_DIR

LOG="./$ANALYSIS_DIR/batch_analysis.output"
echo "Running new analysis..."
. ./analysis_etab_new.sh $ETAB > $LOG
echo "Running 2nd new analysis..."
. ./analysis_etab_new2.sh $ETAB >> $LOG
echo "Running old analysis..."
. ./analysis_etab_old.sh $ETAB >> $LOG
echo "Complete"

mv *.data ./$ANALYSIS_DIR
cp shear_force_analysis.nb ./$ANALYSIS_DIR
