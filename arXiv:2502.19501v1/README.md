# Table of contents

This GitHub repository contains the input and raw output data accompanying the manuscript

**Evidence for strongly correlated superconductivity in La<sub>3</sub>Ni<sub>2</sub>O<sub>7</sub> from first principles**

by Daan Verraes, Tom Braeckevelt, Nick Bultinck, Veronique Van Speybroeck.

This work was submitted on *[arXiv](https://arxiv.org/abs/2502.19501)*.


The data presented here is licensed under the CC BY-SA 4.0 international license, a copy of which can be found [here](https://creativecommons.org/licenses/by-sa/4.0/). Under this license, you can copy and redistribute the material in any medium or format as long as you give appropriate credit, provide a link to the license, and indicate if changes were made.

Additional information concerning the data is available upon request from the authors. Please send a mail to Veronique.VanSpeybroeck@UGent.be for more information.

## Software
All simulations were performed with VASP 6.4.2 interfaced with Wannier90 3.1.0. The pre- and post-processing Python scripts are written in Python 3.12. The VASP calculations utilized POTCAR (PAW_PBE) files, which contain pseudopotentials based on the Projector Augmented-Wave method with the Perdew-Burke-Ernzerhof exchange-correlation functional. These files are proprietary and not publicly accessible.

## Structural Optimization  

These scripts automate the submission and post-processing of **ionic relaxation simulations** in VASP under varying pressure conditions. The `loop.sh` script iterates over a predefined range of pressures, creating and submitting jobs for each case.  

### Workflow:  
1. **Loop Over Pressures**: The script generates simulations for pressures from **0 to 1000 kB**.  
2. **Directory Setup**: For each pressure, a new directory is created, and necessary input files are copied from a template.  
3. **Parameter Modification**:  
   - `PSTRESS` in **INCAR** is updated to reflect the pressure condition.  
   - Job names in **submitscript.sh** are modified accordingly.  
4. **Job Submission**: Each simulation is submitted using `qsub submitscript.sh`.  

### VASP Optimization Setup:  
- **IBRION = 2** → Conjugate gradient algorithm for ionic relaxation.  
- **ISIF = 3** → Allows relaxation of both atomic positions and cell shape.  
- **PSTRESS** → External pressure is dynamically set for each run.  
- **EDIFF = 1E-08** → Tight SCF convergence for accurate energy and forces.  
- **ENCUT = 600 eV** → Plane-wave cutoff energy for reliable results.  

This setup enables efficient structural optimization across different pressures using **first-principles DFT calculations**.



## Electronic Structure and Dimer Model

## Ab initio molecular dynamics (AIMD)

These scripts automate the submission and post-processing of **NPT (constant pressure & temperature) AIMD simulations** for different pressure and temperature conditions using VASP. The loop.sh script iterates over **a range of pressures and temperatures**, creating and submitting simulation jobs for each combination:

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

4. **Main plotting code: plotting.py:**
   - Saves geometric and thermodynamic parameters from the XDATCAR and OUTCAR files
