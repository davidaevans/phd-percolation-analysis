#!/bin/bash
datafile="collated.dat"
dirs=$(ls -d rho*)
sequence=($(seq -w $1 $2 $3))
count=$((0))
rm $datafile
echo "#density shell xmean xrms ymean yrms xymean xyrms meanclustersize" >> $datafile
for dir in $dirs
do

	#cd $dir/datafiles*
	data=`awk '{if (NR==9) print $0}' $dir/shellseries.dat`
	echo ${sequence[$count]} $data >> $datafile
	#cd ../../
	((count+=1))
done
