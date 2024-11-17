#!/bin/bash 
make
    if [ $# -eq 8 ]
then
mkdir -p $8
    ./OneSiteASEPCloning.exe  $1 $2 $3 $4 $5 $6 $7 $8
else
    ./OneSiteASEPCloning.exe
fi
