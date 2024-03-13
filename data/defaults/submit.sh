#!/bin/bash

# Request resources:

#SBATCH -n 1           # 1 CPU core
#SBATCH --mem=8G     # 1GB RAM
#SBATCH --time=24:0:0   #6 hours (hours:mins:secs)

#SBATCH --mail-user=david.evans@durham.ac.uk
#SBATCH --mail-type=END,FAIL


../../sc > run-analysis.out

mkdir datafiles

mv *.sh *.out options  datafiles
mv MC* PHS* ../files

#cp -r block* ../data-*
for block in $(ls -d block*);
do
	cp -r $block ../data*
done
