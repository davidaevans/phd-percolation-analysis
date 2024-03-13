#!/bin/bash

files=$(ls files/)
count=$((0))
sequence=($(seq -w $1 $2 $3))
pe=$4
L=$5

mkdir data-L_${L}-Pe_${pe}

for file in $files;
do
	#datadir="rho0.${sequence[$count]}"
	den=${sequence[$count]}
	dendir="PHS-SQUARE-L_${L}-Dr_C-Fa_24-KT_1.5-Dt_1.5-rho_${den}-Pe_${pe}"
	
	if [ ! -d $dendir ]; then
	
		mkdir $dendir && cd $dendir
		
		#Move trajectory file to directory
		mv ../files/$file .
		echo "$file --> $dendir"

		#Copy setup and default files
		echo "Copying default files to $datadir"
		cp ../defaults/* .
		echo "filename $file" >> options 
		numframes=$(./numframes.sh $file)
		echo "sweeps $numframes" >> options
		echo "Number of frames: $numframes"
		
		#Run analysis
		#echo "./run.sh ${den} ${pe}" >> $datadir/submit.sh
		#submitfile="submit-0.${sequence[$count]}.sh"
		#mv $datadir/submit.sh $datadir/$submitfile
		
		
		#numframes=$(defaults/numframes.sh $datadir/$file)
		#echo "sweeps $numframes" >> $datadir/options
		#echo "Number of frames: $numframes"
		
		echo "Submitting..."
		sbatch submit.sh
		
		cd ../
		((count+=1))
	fi
done

cd data-L_${L}-Pe_${pe}
cp ../scripts/* .
touch lammps-L_${L}-Pe_${pe}-average.dat
