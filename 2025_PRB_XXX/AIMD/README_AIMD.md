# NPT AIMD Job Submission Script

This script automates the submission of **NPT (constant pressure & temperature) AIMD simulations** for different pressure and temperature conditions using VASP.

## How the For-Loop in loop.sh submit VASP jobs

The script iterates over **a range of pressures and temperatures**, creating and submitting simulation jobs for each combination:

1. **Looping Over Parameters:**  
   - `pressures=($(seq 0 50 1000))` → Generates pressures from **0 to 100** in steps of **5 GPa**.  
   - `temperatures=($(seq 10 10 100))` → Generates temperatures from **10 K to 100 K** in steps of **10 K**.

2. **For each pressure-temperature pair:**
   - A new directory (`P<pressure>_T<temperature>`) is created.
   - The template files (including `INCAR`) are copied into this directory.
   - The script **modifies** key parameters in the `INCAR` file:
     - `PSTRESS = <pressure>` → Updates the **external pressure**.
     - `TEBEG = <temperature>` and `TEEND = <temperature>` → Sets the **initial and final temperatures**.
   - The job name in `submitscript.sh` is updated for clarity.
   - The VASP job is submitted via `qsub submitscript.sh`.

3. **VASP AIMD simulation:**

   The **INCAR** file contains VASP settings for the **NPT-AIMD simulation**. The key parameters affected by the loop:
   
   - **PSTRESS**: External pressure (in kB) is set dynamically for each job.
   - **TEBEG / TEEND**: Defines the starting and ending temperature, ensuring a constant temperature MD.
   - Other parameters (e.g., **MDALGO=3**, **SMASS=-3**, **PMASS=2000**) define the **Langevin thermostat & barostat** to control temperature and pressure during MD.
   - Other tags are explained in the INCAR
   
   This script enables efficient, automated sampling of pressure-temperature conditions for **ab initio molecular dynamics** in VASP.

## Post-processing Code

1. **Main plotting code: plotting.py:**
   - Saves geometric and thermodynamic parameters from the XDATCAR and OUTCAR files
