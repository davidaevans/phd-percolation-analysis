#!/bin/bash

for dir in $(ls -d rho*);
do
	#cd $dir
	#mv *.out *.sh options datafiles-*
	#mv PHS* ../files
	#cd ../
	cp -r $dir/block* .
done

for block in $(ls -d block*);
do
	cp defaults/collate.sh $block
	cd $block
	./collate.sh
	cd ../
done
