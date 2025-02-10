#!/bin/bash

template_dir="template"
pressures=($(seq 0 50 1000))
temperatures=($(seq 10 10 100))

STARTDIR=$PBS_O_WORKDIR
cd $STARTDIR

# Check if the template directory exists
if [ ! -d "$template_dir" ]; then
    echo "Error: Template directory '$template_dir' does not exist."
    exit 1
fi

for pressure in "${pressures[@]}"; do
    for temperature in "${temperatures[@]}"; do
        dir_name="P${pressure}_T${temperature}"
        mkdir -p "$dir_name"

        cp -r "$template_dir/"* "$dir_name/"

        # Define paths for the INCAR and submitscript.sh files in the new directory
        incar_file="$dir_name/INCAR"
        submit_script="$dir_name/submitscript.sh"

        # Update the PSTRESS value and TEBEG/TEEND in the INCAR file
        if [ -f "$incar_file" ]; then
            sed -i "/^PSTRESS/c\PSTRESS = $pressure" "$incar_file"
            sed -i "/^TEBEG/c\TEBEG = $temperature" "$incar_file"
            sed -i "/^TEEND/c\TEEND = $temperature" "$incar_file"
        else
            echo "Warning: INCAR file not found in $dir_name. Skipping INCAR updates."
        fi

        # Update the job name in submitscript.sh
        if [ -f "$submit_script" ]; then
            sed -i "/^#PBS -N/c\#PBS -N P${pressure}_T${temperature}" "$submit_script"
        else
            echo "Warning: submitscript.sh file not found in $dir_name. Skipping job name update."
        fi

        cd "$dir_name"
        qsub submitscript.sh
        cd ..
    done
done
