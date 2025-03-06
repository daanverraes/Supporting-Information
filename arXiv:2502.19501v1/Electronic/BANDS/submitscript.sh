#!/bin/bash
#
#PBS -N 2_BANDS
#PBS -m a
#PBS -l walltime=04:00:00
#PBS -l nodes=1:ppn=8
#PBS -l mem=100gb
#PBS -A 2024_079
#

STARTDIR=$PBS_O_WORKDIR
export I_MPI_COMPATIBILITY=4

module purge
module load VASP/6.4.2-gomkl-2021a-VASPsol-20210413-vtst-197-Wannier90-3.1.0
module load vsc-mympirun

cd $STARTDIR
echo "PBS: $PBS_ID"
echo "loaded modules : " `module list` > out.dat

echo "Job started at : "`date` >> out.dat
mympirun vasp_std >> out.dat
echo "Job ended at : "`date` >> out.dat
