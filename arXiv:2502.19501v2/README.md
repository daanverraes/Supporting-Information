# Table of contents

This GitHub repository contains all input data necessary to reproduce the results presented in the manuscript:

**Evidence for strongly correlated superconductivity in La<sub>3</sub>Ni<sub>2</sub>O<sub>7</sub> from first principles**

by Daan Verraes, Tom Braeckevelt, Nick Bultinck and Veronique Van Speybroeck.

This work was submitted on *[arXiv](https://arxiv.org/abs/2502.19501)*.

The data presented here is licensed under the CC BY-SA 4.0 international license, a copy of which can be found [here](https://creativecommons.org/licenses/by-sa/4.0/). Under this license, you can copy and redistribute the material in any medium or format as long as you give appropriate credit, provide a link to the license, and indicate if changes were made.

Additional information concerning the data is available upon request from the authors. Please send a mail to Veronique.VanSpeybroeck@UGent.be for more information.

## Software
All simulations were performed with VASP 6.4.2 interfaced with Wannier90 3.1.0. The pre- and post-processing Python scripts are written in Python 3.12. The VASP calculations utilized POTCAR (PAW_PBE) files: La (5s2 5p6 5d1 6s2), Ni (3d9 4s1), and O (2s2 2p4), which contain pseudopotentials based on the Projector Augmented-Wave method with the Perdew-Burke-Ernzerhof exchange-correlation functional. These files are proprietary and not publicly accessible.

## Structural Optimization  

This setup enables efficient structural optimization across different pressures using **first-principles DFT calculations**.
The available scripts automate the submission and post-processing of **ionic relaxation simulations** in VASP under varying pressure conditions with the following specifications:

- **IBRION = 2** → Conjugate gradient algorithm for ionic relaxation.  
- **ISIF = 3** → Allows relaxation of both atomic positions and cell shape.  
- **PSTRESS** → External pressure is dynamically set for each hydrostatic pressure.
- **EDIFF = 1E-08** → Tight SCF convergence for accurate energy and forces.  
- **ENCUT = 600 eV** → Plane-wave cutoff energy for reliable results.  

## Electronic Structure and Dimer Model

This setup enables DFT calculations, Wannierization and cRPA calculations (given for the unit cell optimized at 2 GPa). The calculations should be performed in this order. Some of the calculations require the output of previous ones as input files and are not always given. The correct input files for each step are given below (the asterix denotes a calculation-specific file):

**DFT groundstate calculation (SCF)**
- INCAR*
- KPOINTS
- POSCAR
- POTCAR

**DFT Long-wave limit (LWL)**
- INCAR*
- KPOINTS
- POSCAR
- POTCAR
- WAVECAR (from SCF)

**DFT Band structure (BANDS)**
- INCAR*
- KPOINTS*
- POSCAR
- POTCAR
- WAVECAR (from LWL)
- CHGCAR (from LWL)

**Wannierization MLWFS**
- INCAR*
- KPOINTS
- POSCAR
- POTCAR
- WAVECAR (from LWL)
- WAVEDER (from LWL)

**cRPA**
- INCAR*
- KPOINTS
- POSCAR
- POTCAR
- WANPROJ (from MLWFS)
- WAVECAR (from LWL)
- WAVEDER (from LWL)

**Tight-binding (TB)**
- hopping.py
- toolbox.py


The interpolation.py script in the MLWFS folder enables the plotting of the band structure superimposed by the Wannier interpolated bands. If any intermediate output files are requested, feel free to contact us.

## Ab initio molecular dynamics (AIMD)

These scripts allow **NPT (constant pressure & temperature) AIMD simulations** for different pressure and temperature conditions using VASP. Geometric and thermodynamic parameters are concatenated in the XDATCAR and OUTCAR files. The **INCAR** file contains VASP settings for the **NPT-AIMD simulation**. The key parameters:

   - **PSTRESS**: External pressure (in kB) is set dynamically for each job.
   - **TEBEG / TEEND**: Defines the starting and ending temperature, ensuring a constant temperature MD.
   - Other parameters (e.g., **MDALGO=3**, **SMASS=-3**, **PMASS=2000**) define the **Langevin thermostat & barostat** to control temperature and pressure during MD.
   - Other tags are explained as comments in the INCAR
