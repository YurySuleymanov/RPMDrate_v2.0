#!/bin/bash

#$ -o stdout_300_1.txt
#$ -e stderr_300_1.txt
#$ -j y
#$ -N CL_H+H2
#$ -l h_rt=9996:00:00
#$ -l mem_total=1G
#$ -l long
#$ -q long5@node96
#$ -pe singlenode 48

cd
source .bash_profile

source /opt/intel/Compiler/11.0/074/bin/intel64/ifortvars_intel64.sh

# Run RPMD
cd /home/ysuleyma/RPMDrate
python rpmdrate.py /home/ysuleyma/RPMDrate/Examples/H+H2/H+H2.py 300 1 -p 48

# Turn off virtualenv when done
deactivate
