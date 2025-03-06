#!/bin/bash
#
#PBS -N 2_TB
#PBS -m a
#PBS -l walltime=01:00:00
#PBS -A 2024_079
#

STARTDIR=$PBS_O_WORKDIR

module --force purge
#module load matplotlib/3.7.2-gfbf-2023a
#module load ASE/3.23.0-gfbf-2023a
module load env/software/dodrio/cpu_rome
module load matplotlib/3.7.2-gfbf-2023a
module load ASE/3.23.0-gfbf-2023a

module list

export OMP_NUM_THREADS=12
export OMP_PROC_BIND=true

cd $STARTDIR
echo "PBS: $PBS_ID"

ls

echo "SCF Job started at : "`date`

python hopping.py

[ -d "__pycache__" ] && rm -r "__pycache__"
echo "Job ended at : "`date`
