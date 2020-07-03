#!/bin/ksh
#for i in *.f90
#do
#dos2unix $i ${i%.f90}.f
#done

echo "Please give a name of the subdireactory:"
read name
mkdir $name
cp main.f90 $name
cp input.txt $name
cp TMDC.exe $name
cp tmdc.script $name
cp checkexist.sh $name
cp move.sh $name
cd $name
sbatch tmdc.script
