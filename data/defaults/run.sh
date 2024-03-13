#!/bin/bash

density=$1
pe=$2

echo $density
echo $pe

../../sc > run-rho_$1-pe_$pe.out

mkdir datafiles-$1-$2

mv *.sh *.out options  datafiles-$1-$2
mv PHS* ../files
mv datafiles* rho$1
