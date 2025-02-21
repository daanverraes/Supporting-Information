#!/bin/bash
#
#PBS -N make_jobs
#PBS -m a
#PBS -l walltime=01:00:00
#

template_dir="template"
pressures=($(seq 0 10 1000))

STARTDIR=$PBS_O_WORKDIR
cd $STARTDIR

# Check if the template directory exists
if [ ! -d "$template_dir" ]; then
    echo "Error: Template directory '$template_dir' does not exist."
    exit 1
fi

for pressure in "${pressures[@]}"; do
    dir_name="$pressure"
    mkdir -p "$dir_name"

    cp -r "$template_dir/"* "$dir_name/"

    # Define paths for the INCAR and submitscript.sh files in the new directory
    incar_file="$dir_name/INCAR"
    submit_script="$dir_name/submitscript.sh"

    # Update the PSTRESS value in the INCAR file
    if [ -f "$incar_file" ]; then
        sed -i "/^PSTRESS/c\PSTRESS = $pressure" "$incar_file"
    else
        echo "Warning: INCAR file not found in $dir_name. Skipping PSTRESS update."
    fi

    # Update the job name in submitscript.sh
    if [ -f "$submit_script" ]; then
        sed -i "/^#PBS -N/c\#PBS -N LNO_ACC_P_$pressure" "$submit_script"
    else
        echo "Warning: submitscript.sh file not found in $dir_name. Skipping job name update."
    fi
    
    cd $pressure
    qsub submitscript.sh
    cd ..
done

echo "Ionic optimisation jobs submitted for pressures from 40 to 240."
