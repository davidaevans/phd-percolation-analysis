#!/bin/bash

#######################################################################
# This script collates and averages all the block data for LAMMPS     #
# simulations of ABP.                                                 #
# The input is the range of densities in the form of: start step end  #
# To run:                                                             #
# ./collate-and-plot.sh 0.550 0.005 0.560                             #
#######################################################################

##### STEP 1 #####
# For each block, go through and collate the data for each density

# If a filenames file exists, remove it
if [ -f filenames.dat ]; then
	rm filenames.dat
fi

# Collate data for each block
for block in $(ls -d block*);
do
	cd $block
	
	# Copy collate script into each block
	cp ../../defaults/collate.sh .
	./collate.sh $1 $2 $3
	
	# Add resulting collated data filename to file
	echo "$block/collated.dat" >> ../filenames.dat
	
	cd ../
done


##### STEP 2 #####
# Average over blocks
cp ../scripts/average.py .
python3 average.py $(<filenames.dat) > lammps-*


##### STEP 3 #####
# Gather trajectory files
#mkdir trajectories
#currentdir=$(pwd)
#cd ../
#for dir in $(ls -d PHS*);
#do
#	cp $dir/analysisdata/Video $currentdir/trajectories/$dir 
#done
#cd $currentdir
#cp ../files/* trajectories

##### STEP 3 #####
# Make plots
#datafilename=$(ls lammps-*)
#cp ../scripts/analysis-defaults/plot-percprob.gplot .
#gnuplot -e "filename='$datafilename'" plot-percprob.gplot
