import numpy as np
import matplotlib.pyplot as plt
# import spglib
from toolbox import *
import os
from ase.io import read
from ase import atoms


# Move this to toolbox if working:
def point_above_or_below_plane(p1, p2, p3, p4):
    """
    Determine whether a point p4 is above or below the plane defined by points p1, p2, and p3.

    Parameters:
    p1, p2, p3, p4 : np.array
        3D points defined as numpy arrays (e.g., np.array([x, y, z])).

    Returns:
    str
        'above' if p4 is above the plane,
        'below' if p4 is below the plane,
        'on the plane' if p4 lies exactly on the plane.
    """
    # Compute vectors in the plane
    v1 = p2 - p1
    v2 = p3 - p1

    # Compute the normal vector to the plane
    normal = np.cross(v1, v2)

    # Compute vector from p1 to p4
    v_to_p4 = p4 - p1

    # Compute the dot product of the normal vector with the vector to p4
    dot_product = np.dot(normal, v_to_p4)

    # Determine the position of p4 relative to the plane
    if np.isclose(dot_product, 0):
        return "on the plane"
    elif dot_product > 0:
        return "above"
    else:
        return "below"


folder = sorted([f for f in os.listdir('.') if os.path.isdir(f) and f != '__pycache__' and f != 'figures' and f != 'template'])

printing = False

volumes = np.zeros(len(folder))

length_a = np.zeros(len(folder))
length_b = np.zeros(len(folder))
length_c = np.zeros(len(folder))

alpha = np.zeros(len(folder))
beta = np.zeros(len(folder))
gamma = np.zeros(len(folder))

angles_inter = np.zeros(len(folder))
angles_intra = np.zeros(len(folder))
angles_intra2 = np.zeros(len(folder))
angles_octa_x = np.zeros(len(folder))
angles_octa_y = np.zeros(len(folder))
angles_octa_xy = np.zeros(len(folder))
angles_octa_yx = np.zeros(len(folder))

length_octa_1 = np.zeros(len(folder))
length_octa_2 = np.zeros(len(folder))
length_octa_3 = np.zeros(len(folder))
length_octa_4 = np.zeros(len(folder))

energies = np.zeros(len(folder))
enthalpies = np.zeros(len(folder))
pressures = np.zeros(len(folder))

for p in range(len(folder)):

    pressures[p] = folder[p]
    if printing: print("\nPressure = " + str(float(folder[p]))+" GPa")

    path_to_POSCAR = folder[p]+'/CONTCAR'
    cell = read_poscar(path_to_POSCAR)

    lattice_matrix = read_poscar_lattice(path_to_POSCAR)
    a, b, c = lattice_matrix[0], lattice_matrix[1], lattice_matrix[2]

    # Lattice vector lengths and angles:
    length_a[p] = np.linalg.norm(a)
    length_b[p] = np.linalg.norm(b)
    length_c[p] = np.linalg.norm(c)
    
    alpha[p] = np.degrees(np.arccos(np.dot(b, c) / (length_b[p] * length_c[p]))) 
    beta[p] = np.degrees(np.arccos(np.dot(a, c) / (length_a[p] * length_c[p]))) 
    gamma[p] = np.degrees(np.arccos(np.dot(a, b) / (length_a[p] * length_b[p]))) 
    
    atoms = read(path_to_POSCAR)
    pos = atoms.get_positions()

    angles_inter[p] = atoms.get_angle(14, 34, 16, mic=True)
    angles_intra[p] = atoms.get_angle(14, 30, 15, mic=True)
    angles_intra2[p] = atoms.get_angle(12, 22, 13, mic=True)
    
    if point_above_or_below_plane(pos[12], pos[13], pos[12]+a, pos[22]) == "below":
        angles_intra2[p] = 360 - angles_intra2[p]
    
    
    angles_octa_x[p] = atoms.get_angle(30, 14, 33, mic=True)
    angles_octa_y[p] = atoms.get_angle(31, 14, 32, mic=True)
    angles_octa_xy[p] = atoms.get_angle(30, 14, 32, mic=True)
    angles_octa_yx[p] = atoms.get_angle(31, 14, 33, mic=True)

    length_octa_1[p] = atoms.get_distance(21, 13, mic=True)
    length_octa_2[p] = atoms.get_distance(25, 13, mic=True)
    length_octa_3[p] = atoms.get_distance(22, 13, mic=True)
    length_octa_4[p] = atoms.get_distance(27, 13, mic=True)
        
    path_to_OUTCAR = folder[p]+'/OUTCAR'
    with open(path_to_OUTCAR, 'r') as file:
        for line in file:
            if "energy(sigma->0)" in line:
                last_energy_value = float(line.split('=')[-1].strip())

            elif "enthalpy" in line:
                last_enthalpy_value = float(line.split()[-5].strip())

            elif "volume of cell" in line:
                last_volume = float(line.split()[-1].strip())

    energies[p] = last_energy_value

    if float(pressures[p]) > 1e-5:
        enthalpies[p] = last_enthalpy_value
    else:
        enthalpies[p] = last_energy_value

    volumes[p] = last_volume

    if printing: print('Total energy:', last_energy_value, 'eV')

pressures = 0.1*pressures #GPa
data = np.column_stack((pressures, enthalpies, volumes))
np.savetxt("enthalpies_2222.txt", data, header="Pressures Enthalpies Volumes", fmt="%.6f", delimiter="\t")

data = np.column_stack((pressures, length_a, length_b, length_c, angles_inter, angles_intra, angles_intra2, length_octa_1, length_octa_2, length_octa_3, length_octa_4))
np.savetxt("bilayer_geometries.txt", data, header="Pressures a b c Interlayer Intralayer_in Intralayer_out 1 2 3 4", fmt="%.6f", delimiter="\t")

fig = plt.figure(figsize=(8, 7))
#bax = brokenaxes(ylims=((5,6), (19.2,20.2)), hspace=0.1) 
plt.scatter(pressures, length_a, label='a')
plt.scatter(pressures, length_b, label='b')
plt.xlabel("External pressure [GPa]", fontsize=12)
plt.ylabel("Lattice parameter [Å]", fontsize=12)
plt.grid(True)
plt.legend()

plt.title("Lattice parameters (Z=4)")
plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
plt.savefig("figures/cell_ab.png", bbox_inches="tight")

import matplotlib.pyplot as plt

fig = plt.figure(figsize=(8, 7))

# Scatter plot
plt.scatter(pressures, np.abs(length_a-length_b), label='a')
plt.xlabel("External pressure [GPa]", fontsize=12)
plt.ylabel("Lattice parameter [Å]", fontsize=12)
plt.yscale('log')
plt.grid(True)
plt.legend()

plt.title("Lattice parameters (Z=4)")
plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
plt.savefig("figures/cell_ab_logs.png", bbox_inches="tight")
plt.show()

fig = plt.figure(figsize=(8, 7))
#bax = brokenaxes(ylims=((5,6), (19.2,20.2)), hspace=0.1) 
plt.scatter(pressures, length_c, label='c')
plt.xlabel("External pressure [GPa]", fontsize=12)
plt.ylabel("Lattice parameter [Å]", fontsize=12)
plt.grid(True)
plt.legend()
plt.ylim(19,20.5)
plt.title("Lattice parameters (Z=4)")
plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
plt.savefig("figures/cell_c.png", bbox_inches="tight")

plt.figure()
plt.scatter(pressures, alpha, label='r$\alpha$')
plt.scatter(pressures, beta, label='b')
plt.scatter(pressures, gamma, label='c')
plt.xlabel("External pressure [GPa]", fontsize=14)
# plt.ylabel("Ni-O-Ni angle [°]", fontsize=14)
plt.grid()
# plt.minorticks_on()
plt.tight_layout()
plt.savefig("figures/cell_angles", bbox_inches="tight")

plt.figure()
plt.scatter(pressures, angles_inter, linestyle="-")
plt.xlabel("External pressure [GPa]", fontsize=14)
plt.ylabel("Ni-O-Ni angle [°]", fontsize=14)
plt.grid()
plt.minorticks_on()
plt.tight_layout()
plt.title("Interlayer octahedral tilt")
plt.savefig("figures/interlayer_NiONi", bbox_inches="tight")

plt.figure()
plt.scatter(pressures, angles_intra, linestyle="-", label="in")
plt.scatter(pressures, angles_intra2, linestyle="-", label="out")
plt.xlabel("External pressure [GPa]", fontsize=14)
plt.ylabel("Ni-O-Ni angle [°]", fontsize=14)
plt.grid()
plt.legend()
plt.minorticks_on()
plt.tight_layout()
plt.title("Intralayer octahedral tilt")
plt.savefig("figures/intralayer_NiONi", bbox_inches="tight")

plt.figure()
plt.scatter(pressures, length_octa_1, linestyle="-", label='1')
plt.scatter(pressures, length_octa_2, linestyle="-", label='2')
plt.scatter(pressures, length_octa_3, linestyle="-", label='3')
plt.scatter(pressures, length_octa_4, linestyle="-", label='4')
plt.xlabel("External pressure [GPa]", fontsize=14)
plt.ylabel("Ni-O length [Å]", fontsize=14)
plt.legend()
plt.grid()
plt.minorticks_on()
plt.tight_layout()
#plt.title("Intralayer octahedral tilt")
plt.savefig("figures/intralayer_NiO_length_octa", bbox_inches="tight")

plt.figure()
plt.scatter(pressures, volumes, linestyle="-")
plt.xlabel("External pressure [GPa]", fontsize=14)
plt.ylabel(r"Unit cell volume [Å³]", fontsize=14)
plt.grid()
plt.minorticks_on()
plt.tight_layout()
plt.savefig("figures/Volumes", bbox_inches="tight")

plt.figure()
plt.scatter(pressures, energies, linestyle="-")
plt.xlabel("External pressure [GPa]", fontsize=14)
plt.ylabel(r"E$_{tot}$ [eV]", fontsize=14)
plt.grid()
plt.minorticks_on()
plt.tight_layout()
plt.savefig("figures/Energies", bbox_inches="tight")

plt.figure()
plt.scatter(pressures, energies, label='E')
plt.scatter(pressures, enthalpies, label='H')
plt.scatter(pressures, enthalpies-energies, label='PV')
plt.xlabel("External pressure [GPa]", fontsize=14)
plt.ylabel(r"H = E + PV [eV]", fontsize=14)
plt.grid()
plt.legend()
plt.minorticks_on()
plt.tight_layout()
plt.savefig("figures/Enthalpies", bbox_inches="tight")
