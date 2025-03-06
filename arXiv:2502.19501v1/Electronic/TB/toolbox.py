import numpy as np
import spglib
import random

def read_poscar(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()

    scaling_factor = float(lines[1].strip())

    lattice = []
    for i in range(2, 5):
        lattice.append([float(x) for x in lines[i].strip().split()])
    lattice = np.array(lattice) * scaling_factor

    element_counts = [int(x) for x in lines[6].strip().split()]

    numbers = []
    atom_type = 1
    for count in element_counts:
        numbers.extend([atom_type] * count)
        atom_type += 1

    fractional_positions = []
    for line in lines[8:]:
        fractional_positions.append([float(x) for x in line.strip().split()])

    return lattice, fractional_positions, numbers

def read_kpoints(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()

    mesh = [int(x) for x in lines[3].strip().split()]
    
    return mesh

def read_outcar(filename):
    system = None
    N_b = None
    
    with open(filename, 'r') as file:
        for line in file:
            if 'SYSTEM =' in line:
                system = line.split('=')[1].strip()  
            elif 'NBANDS =' in line:
                N_b = int(line.split('=')[1].strip())
    
    return system, N_b

def read_wanproj(filename):
    with open(filename, 'r') as file:
        for line in file:
            if line.strip() and not line.startswith(' #'):
                columns = line.split()
                if len(columns) >= 4:
                    N_w = columns[3]
                    break
    return int(N_w)

# WRONG
# def read_angle(filename,ion1_idx, ion2_idx, ion3_idx):
#     with open(filename, 'r') as file:
#         lines = file.readlines()
    
#     scale = 1 #float(lines[1].strip())  # Scaling factor
    
#     # Extract lattice vectors and apply scaling
#     a_vec = np.array([float(x) for x in lines[2].split()]) * scale
#     b_vec = np.array([float(x) for x in lines[3].split()]) * scale
#     c_vec = np.array([float(x) for x in lines[4].split()]) * scale
    
#     # Read atomic positions in fractional coordinates (starting after the header)
#     positions = []
#     start_idx = 8
#     for line in lines[start_idx:]:
#         if line.strip():
#             positions.append(np.array([float(x) for x in line.split()[:3]]))

#     # Convert fractional coordinates to Cartesian
#     cartesian_coords = [fractional_to_cartesian(pos, a_vec, b_vec, c_vec) for pos in positions]
    
#     # Get vectors between the ions
#     vec_ion2_ion1 = cartesian_coords[ion2_idx] - cartesian_coords[ion1_idx]
#     vec_ion3_ion1 = cartesian_coords[ion3_idx] - cartesian_coords[ion1_idx]
    
#     # Calculate the angle between vectors using the dot product
#     cos_theta = np.dot(vec_ion2_ion1, vec_ion3_ion1) / (np.linalg.norm(vec_ion2_ion1) * np.linalg.norm(vec_ion3_ion1))
#     angle_rad = np.arccos(cos_theta)
#     angle_deg = 180-np.degrees(angle_rad)
    
#     return(angle_deg)

# Function to convert fractional to Cartesian coordinates
def fractional_to_cartesian(fractional_coord, a_vec, b_vec, c_vec):
    return fractional_coord[0] * a_vec + fractional_coord[1] * b_vec + fractional_coord[2] * c_vec


def Perturb_my_POSCAR(file_path, output_file_path, max_shift):
    modified_lines = []
    with open(file_path, 'r') as file:
        lines = file.readlines()
        in_direct_coords = False  # Flag to check if we are in the Direct coordinate section

        for line in lines:
            # Detect the start of the Direct coordinates
            if line.strip().lower().startswith("direct"):
                in_direct_coords = True
                modified_lines.append(line)
                continue

            # Check if we are in the Direct coordinates section and if the line contains coordinates
            if in_direct_coords:
                # Skip empty lines or non-coordinate lines
                if line.strip() == "":
                    modified_lines.append(line)
                    continue

                # Parse the line and add a random float to each coordinate value
                coords = line.split()
                if len(coords) >= 3:  # Expecting at least three values for x, y, and z coordinates
                    modified_coords = [
                        f"{float(coord) + random.uniform(0.0, max_shift):.16f}" for coord in coords[:3]
                    ]
                    modified_line = "  ".join(modified_coords) + "\n"
                    modified_lines.append(modified_line)
                else:
                    # If line does not contain 3 coordinates, append it unchanged
                    modified_lines.append(line)
            else:
                # Before reaching Direct coordinates, add lines unmodified
                modified_lines.append(line)

    # Write the modified lines to the output file
    with open(output_file_path, 'w') as output_file:
        output_file.writelines(modified_lines)



import numpy as np

def read_poscar_lattice(file_path):
    """Reads the lattice vectors from a POSCAR file."""
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    # Get the scaling factor from line 2
    scale = float(lines[1].strip())
    
    # Read lattice vectors from lines 3 to 5 and apply scaling
    lattice_vectors = []
    for i in range(3):
        vector = np.array([float(x) for x in lines[2 + i].strip().split()])
        lattice_vectors.append(vector * scale)
    
    # Convert lattice_vectors list to a numpy array
    lattice_matrix = np.array(lattice_vectors)
    return lattice_matrix

def eigsolve_lattice_vectors(lattice_matrix):
    """Computes the eigenvalues and eigenvectors of the lattice matrix."""
    eigenvalues, eigenvectors = np.linalg.eig(lattice_matrix)
    return eigenvalues, eigenvectors

