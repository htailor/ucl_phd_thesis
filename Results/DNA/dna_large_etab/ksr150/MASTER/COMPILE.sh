#!/bin/sh
./clear_all.sh
EXE_FILE=Nucleation
cp $EXE_FILE "$EXE_FILE.bak" 
cd ./src
clear
echo "Compiling new Nucelation file..."
echo "==============================================================================="
make -f Makefile
echo "==============================================================================="
echo "COMPLETE."
rm *.o
mv $EXE_FILE ..
cd ..
#echo "linking Nucleation file to all enviroments"
#cp $EXE_FILE ../Nucleation_1
#cp $EXE_FILE ../Nucleation_2
#cp $EXE_FILE ../Nucleation_3
#cp $EXE_FILE ../Nucleation_4
#cp $EXE_FILE ../Nucleation_5
#cp $EXE_FILE ../Nucleation_6
#echo "COMPLETE"
echo
