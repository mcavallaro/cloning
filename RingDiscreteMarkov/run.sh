#!/bin/bash 
make
if [ $# -eq 4 ]
then
mkdir -p $4
./BDMClo.exe  $1 $2 $3 $4
else
./BDMClo.exe
fi
