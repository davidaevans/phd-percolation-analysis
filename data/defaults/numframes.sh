#!/bin/bash

file=$1

numlines=$(wc -l $file | awk '{print $1}')
npart=$(head -4 $1 | tail -1)
numframes=$(awk -v numlines=$numlines -v npart=$npart 'BEGIN{print(numlines/(npart+9))}')
echo $numframes

