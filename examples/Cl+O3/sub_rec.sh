#!/bin/bash 
#$ -o outwT250nb16_rec.txt
#$ -e errwT250nb16_rec.txt
#$ -j y 
#$ -N ClO3_rec250_16
#$ -l h_rt=9996:00:00 
#$ -l mem_total=1G
#$ -l long 
#$ -pe singlenode 8 
cd 
source .bash_profile 
source /opt/intel/Compiler/11.0/074/bin/intel64/ifortvars_intel64.sh 
   
   
# Run RPMD 
cd /home/ysuleyma/RPMDrate/ 
python rpmdrate.py /home/ysuleyma/RPMDrate/Examples/Cl+O3final/input_rec16.py 250 16 -p 8
 
