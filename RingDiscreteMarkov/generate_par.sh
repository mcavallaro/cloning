#!/bin/bash 
make
if [ $# -eq 1 ]
then
rm $1
for s in $(seq -1.0 0.05 1)
do
    echo "1000 1000 5 $s RingDiscrete directory/" >> $1
done

fi
