import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib import rcParams
import re
from scipy.interpolate import griddata
from ase.io import read

rcParams.update({'font.size': 24}) 
rcParams['font.family'] = 'Times New Roman'
rcParams["mathtext.fontset"] = "cm"

def signed_angle_between_three_points(p1, p2, p3, normal=None):

    AB = p2 - p1
    BC = p3 - p2

    # normal = np.cross(AB, BC)
    # normal /= np.linalg.norm(normal)
    
    dot_product = np.dot(AB, BC)
    AB_mag = np.linalg.norm(AB)
    BC_mag = np.linalg.norm(BC)
    angle = np.arccos(dot_product / (AB_mag * BC_mag))
    
    # Cross product to determine the sign
    cross_product = np.cross(AB, BC) 
    if np.dot(normal, cross_product) < 0: 
        angle = -angle 
    
    return 180-angle*180/np.pi

def read_xdatcar(filepath):
    with open(filepath, 'r') as f:
        lines = f.readlines()

    lattice_vectors = []
    atom_positions = []
    timesteps = []

    scale_factor = float(lines[1].strip())
    species_counts = list(map(int, lines[6].split()))
    total_atoms = sum(species_counts)

    i = 7
    while i < len(lines):
        if lines[i].startswith("Direct configuration"):
            timesteps.append(lines[i].strip().split('=')[-1])

            current_lattice = np.array([
                list(map(float, lines[i - 5].split())),
                list(map(float, lines[i - 4].split())),
                list(map(float, lines[i - 3].split()))
            ]) * scale_factor
            lattice_vectors.append(current_lattice)

            positions = []
            for j in range(total_atoms):
                positions.append(list(map(float, lines[i + 1 + j].split())))
            atom_positions.append(np.array(positions))
            i += total_atoms + 1
        else:
            i += 1
    
    return lattice_vectors, atom_positions, np.array(timesteps)

def fractional_to_cartesian(lattice_vectors, fractional_positions):
    cartesian_positions = []
    for i, (lattice, positions) in enumerate(zip(lattice_vectors, fractional_positions)):
        cartesian = np.dot(positions, lattice)
        cartesian_positions.append(cartesian)
    return cartesian_positions

def angle_between_vectors(v1, v2):
    dot_product = np.dot(v1, v2)
    norms = np.linalg.norm(v1) * np.linalg.norm(v2)
    cos_theta = dot_product / norms
    theta_radians = np.arccos(np.clip(cos_theta, -1.0, 1.0))
    return np.degrees(theta_radians)

def read_outcar(outcar_path):
    ekin_values = []
    etotal_values = []
    step_counter = 0
    pressures = []
    temperatures = []
    volumes = []
    free_energies = []

    with open(outcar_path, 'r') as file:
        for line in file:

            if "kinetic energy EKIN" in line:
                ekin_values.append(float(line.split()[-1]))
            elif "total energy   ETOTAL" in line:
                etotal_values.append(float(line.split()[-2]))

            elif "FREE ENERGIE OF THE ION-ELECTRON SYSTEM" in line:
                step_counter += 1

            elif "free  energy   TOTEN  =" in line:
                match = re.search(r"free\s+energy\s+TOTEN\s+=\s+([-.\d]+)\s+eV", line)
                if match:
                    free_energies.append(float(match.group(1)))
        
            elif "temperature" in line:
                temp_match = re.search(r"\(temperature\s+([\d.]+)\s+K\)", line)
                if temp_match:
                    temperature = float(temp_match.group(1))
                    temperatures.append(temperature)

            elif "total pressure  =" in line:
                match = re.search(r"total pressure\s+=\s+([-.\d]+)\s+kB", line)
                if match:
                    pressures.append(float(match.group(1)))

            elif "volume of cell :" in line:
                match = re.search(r"volume of cell :\s+([\d.]+)", line)
                if match:
                    volumes.append(float(match.group(1)))

    return temperatures, pressures, volumes, free_energies, ekin_values, etotal_values

primary_color = (0/255, 128/255, 192/255) 
secondary_color = (139/255, 0/255, 0/255) 

root_directory = "."

avg_data = {
    "dists_x": {}, "dists_y": {}, "dists_u": {}, "dists_d": {}, "dists_O3": {}, "dists_O1": {},
    "diff_dists_xy": {}, "diff_dists_colinear": {}, "diff_dists_colinear2": {}, 
    "angles_intra": {}, "angles_intra2": {}, "angles_inter": {}, "signed_angles_inter": {}, "STD_inter": {},
    "ONiO_x": {}, "ONiO_y": {}, "ONiO_diff": {},
    "a": {}, "b": {}, "c": {}, "diff_a_b": {}, "diff_a_b_2": {},
    "alpha": {}, "beta": {}, "gamma": {},
    "temperature": {}, "pressure": {},
    "free_energy": {}, "volume": {},
    "kinetic_energy": {}, "total_energy": {},
    "timestep": {}
}

for folder in os.listdir(root_directory):
    folder_path = os.path.join(root_directory, folder)
    if os.path.isdir(folder_path):
        try:
            parts = folder.split('_')
            press = int(parts[0][1:]) // 10
            temp = int(parts[1][1:])
            
        except (IndexError, ValueError):
            print(f"Skipping folder with invalid format: {folder}")
            continue
        
        outcar_file = os.path.join(folder_path, "OUTCAR")
        xdatcar_file = os.path.join(folder_path, "XDATCAR")
        
        if not os.path.isfile(outcar_file) or not os.path.isfile(xdatcar_file) or os.path.getsize(xdatcar_file) == 0:
            for key, data_dict in avg_data.items():
                avg_data[key][(press, temp)] = np.nan
            continue

        temperatures, pressures, volumes, free_energies, ekin_values, etotal_values = read_outcar(outcar_file)
        lattice_vectors, fractional_positions, timesteps = read_xdatcar(xdatcar_file)
        cartpos = fractional_to_cartesian(lattice_vectors, fractional_positions)
        
        a = np.array([np.linalg.norm(lattice_vectors[t][0]) for t in range(len(timesteps))])
        b = np.array([np.linalg.norm(lattice_vectors[t][1]) for t in range(len(timesteps))])
        c = np.array([np.linalg.norm(lattice_vectors[t][2]) for t in range(len(timesteps))])
        
        alpha = np.array([angle_between_vectors(lattice_vectors[t][1], lattice_vectors[t][2]) for t in range(len(timesteps))])
        beta = np.array([angle_between_vectors(lattice_vectors[t][0], lattice_vectors[t][2]) for t in range(len(timesteps))])
        gamma = np.array([angle_between_vectors(lattice_vectors[t][0], lattice_vectors[t][1]) for t in range(len(timesteps))])

        dists_x = np.array([np.linalg.norm(cartpos[t][30] - cartpos[t][14]) for t in range(len(timesteps))])
        dists_y = np.array([np.linalg.norm(cartpos[t][31] - cartpos[t][14]) for t in range(len(timesteps))])
        dists_u = np.array([np.linalg.norm(cartpos[t][22] - cartpos[t][13]) for t in range(len(timesteps))])
        dists_d = np.array([np.linalg.norm(cartpos[t][25] - cartpos[t][13]) for t in range(len(timesteps))])
        
        dists_O1 = np.array([np.linalg.norm(cartpos[t][21] - cartpos[t][13]) for t in range(len(timesteps))])
        dists_O3 = np.array([np.linalg.norm(cartpos[t][27] - cartpos[t][13]) for t in range(len(timesteps))])
        
        filename = folder_path+"/XDATCAR"
        traj = read(filename, index=':') 
        
        angles_inter = np.zeros(len(traj))
        angles_intra = np.zeros(len(traj))
        angles_intra2 = np.zeros(len(traj))
        
        for p, atoms in enumerate(traj):
            angles_inter[p] = atoms.get_angle(14, 34, 16, mic=True)
            angles_intra[p] = atoms.get_angle(14, 30, 15, mic=True)
            angles_intra2[p] = atoms.get_angle(12, 22, 13, mic=True)
            
            #if point_above_or_below_plane(pos[12], pos[13], pos[12] + a, pos[22]) == "below":
            #    angles_intra2[p] = 360 - angles_intra2[p]
        
        ONiO_x = np.abs(180-np.array([signed_angle_between_three_points(cartpos[t][22], cartpos[t][13], cartpos[t][23], normal=lattice_vectors[t][0]) for t in range(len(timesteps))]))
        ONiO_y = np.abs(180-np.array([signed_angle_between_three_points(cartpos[t][23], cartpos[t][13], cartpos[t][25], normal=lattice_vectors[t][0]) for t in range(len(timesteps))]))
        ONiO_diff = np.abs(ONiO_x - ONiO_y)

        # Plot Ni-O bond lengths
        plt.figure()
        plt.plot(dists_x, linestyle="-", label='x', color=primary_color)
        plt.plot(dists_y, linestyle="-", label='y', color=secondary_color)
        plt.xlabel("time [fs]", fontsize=14)
        plt.ylabel("Ni-O length [Å]", fontsize=14)
        plt.legend()
        plt.grid()
        plt.minorticks_on()
        plt.tight_layout()
        plt.title(f"T={temp}K, P={press}GPa: Intralayer Ni-O bond lengths")
        plot_file_path = os.path.join(folder_path, "intralayer_NiO_length.png")
        plt.savefig(plot_file_path, bbox_inches="tight")
        plt.close()        
        
        # Plot Ni-O bond lengths up vs down
        plt.figure()
        plt.plot(dists_u, linestyle="-", label='u')
        plt.plot(dists_d, linestyle="-", label='d')
        plt.plot(dists_O1, linestyle="-", label='O1')
        plt.plot(dists_O3, linestyle="-", label='O3')
        plt.xlabel("time [fs]", fontsize=14)
        plt.ylabel("Ni-O length [Å]", fontsize=14)
        plt.legend()
        plt.grid()
        plt.minorticks_on()
        plt.tight_layout()
        plt.title(f"T={temp}K, P={press}GPa: Octahedral Ni-O bond lengths")
        plot_file_path = os.path.join(folder_path, "intralayer_NiO_all_lengths.png")
        plt.savefig(plot_file_path, bbox_inches="tight")
        plt.close()

        # Plot Ni-O-Ni angles (Intra and Inter)
        plt.figure()
        plt.plot(angles_intra, linestyle="-", label='Intra', color=primary_color)
        plt.plot(angles_intra2, linestyle="-", label='Intra', color=primary_color)
        plt.plot(angles_inter, linestyle="-", label='Inter', color=secondary_color)
        plt.xlabel("time [fs]", fontsize=14)
        plt.ylabel("Ni-O-Ni angle [°]", fontsize=14)
        plt.legend()
        plt.grid()
        # plt.ylim([155,180])
        plt.minorticks_on()
        plt.tight_layout()
        plt.title(f"T={temp}K, P={press}GPa: Ni-O-Ni Angles")
        plot_file_path = os.path.join(folder_path, "NiONi_angles.png")
        plt.savefig(plot_file_path, bbox_inches="tight")
        plt.close()

        # Calculate lattice parameters and angles
        a = np.array([np.linalg.norm(lattice_vectors[t][0]) for t in range(len(timesteps))])
        b = np.array([np.linalg.norm(lattice_vectors[t][1]) for t in range(len(timesteps))])
        c = np.array([np.linalg.norm(lattice_vectors[t][2]) for t in range(len(timesteps))])

        fig = plt.figure()
        # bax = brokenaxes(ylims=((5.1, 5.6), (19.5, 20.0)), hspace=0.1)
        plt.plot(a, linestyle="-", label=r'a', color=primary_color)
        plt.plot(b, linestyle="-", label=r'b', color=secondary_color)
        # plt.plot(c, linestyle="-", label=r'c', color='lightgreen')
        plt.xlabel("Time [fs]", fontsize=14)
        plt.ylabel("Lattice Parameter [Å]", fontsize=14)
        plt.grid(True)
        plt.legend(fontsize=14)
        plt.title(f"T={temp} K, P={press} GPa: Lattice Parameters", fontsize=16)
        plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
        plot_file_path = os.path.join(folder_path, "cell_lengths.png")
        plt.savefig(plot_file_path, bbox_inches="tight")
        plt.close()

        alpha = np.array([angle_between_vectors(lattice_vectors[t][1], lattice_vectors[t][2]) for t in range(len(timesteps))])
        beta = np.array([angle_between_vectors(lattice_vectors[t][0], lattice_vectors[t][2]) for t in range(len(timesteps))])
        gamma = np.array([angle_between_vectors(lattice_vectors[t][0], lattice_vectors[t][1]) for t in range(len(timesteps))])

        plt.figure()
        plt.plot(alpha, linestyle="-", label=r'$\alpha$', color=primary_color)
        plt.plot(beta, linestyle="-", label=r'$\beta$', color=secondary_color)
        plt.plot(gamma, linestyle="-", label=r"$\gamma$", color='lightgreen')
        plt.xlabel("time [fs]", fontsize=14)
        plt.ylabel("angle [°]", fontsize=14)
        plt.legend()
        plt.grid()
        plt.minorticks_on()
        plt.tight_layout()
        plt.title(f"T={temp}K, P={press}GPa: Unit cell angles")
        plot_file_path = os.path.join(folder_path, "cell_angles.png")
        plt.savefig(plot_file_path, bbox_inches="tight")
        plt.close()

        # Plot temperature, pressure, free energy, volume, kinetic energy, and total energy
        plt.figure()
        plt.plot(temperatures, label="Temperature (K)", color=primary_color)
        plt.xlabel("Timestep")
        plt.ylabel("Temperature (K)")
        plt.title(f"T={temp}K, P={press}GPa: Temperature")
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plot_file_path = os.path.join(folder_path, "temperature.png")
        plt.savefig(plot_file_path)
        plt.close()

        plt.figure()
        plt.plot(0.1*np.array(pressures), label="Pressure (GPa)", color='red')
        plt.xlabel("Timestep")
        plt.ylabel("Pressure (GPa)")
        plt.title(f"T={temp}K, P={press}GPa: Pressure")
        plt.grid(True)
        plt.tight_layout()
        plot_file_path = os.path.join(folder_path, "pressure.png")
        plt.savefig(plot_file_path)
        plt.close()

        # Plot free energy, volume, kinetic energy, total energy
        plt.figure()
        plt.plot(free_energies, color=primary_color)
        plt.xlabel("time [fs]", fontsize=14)
        plt.ylabel("Free Energy [eV]", fontsize=14)
        plt.grid()
        plt.minorticks_on()
        plt.tight_layout()
        plt.title(f"T={temp}K, P={press}GPa: Total Free energy")
        plot_file_path = os.path.join(folder_path, "free_energies.png")
        plt.savefig(plot_file_path, bbox_inches="tight")
        plt.close()

        plt.figure()
        plt.plot(volumes, linestyle="-", color=primary_color)
        plt.xlabel("time [fs]", fontsize=14)
        plt.ylabel("Volume [Å³]", fontsize=14)
        plt.grid()
        plt.minorticks_on()
        plt.tight_layout()
        plt.title(f"T={temp}K, P={press}GPa: Cell volume")
        plot_file_path = os.path.join(folder_path, "volumes.png")
        plt.savefig(plot_file_path, bbox_inches="tight")
        plt.close()

        plt.figure()
        plt.plot(ekin_values, linestyle="-", color=primary_color)
        plt.xlabel("time [fs]", fontsize=14)
        plt.ylabel("Kinetic energy [eV]", fontsize=14)
        plt.grid()
        plt.minorticks_on()
        plt.tight_layout()
        plt.title(f"T={temp}K, P={press}GPa: Kinetic Energy")
        plot_file_path = os.path.join(folder_path, "kinetic_energies.png")
        plt.savefig(plot_file_path, bbox_inches="tight")
        plt.close()

        plt.figure()
        plt.plot(etotal_values, linestyle="-", color=primary_color)
        plt.xlabel("time [fs]", fontsize=14)
        plt.ylabel("Total energy [eV]", fontsize=14)
        plt.grid()
        plt.minorticks_on()
        plt.tight_layout()
        plt.title(f"T={temp}K, P={press}GPa: Total Energy")
        plot_file_path = os.path.join(folder_path, "total_energies.png")
        plt.savefig(plot_file_path, bbox_inches="tight")
        plt.close()
 
        if len(timesteps) <= 2000:
            for key, data_dict in avg_data.items():
                avg_data[key][(press, temp)] = np.nan
            continue
            
        avg_data["dists_x"][(press, temp)] = np.mean(dists_x[2000:])
        avg_data["dists_y"][(press, temp)] = np.mean(dists_y[2000:])
        avg_data["diff_dists_xy"][(press, temp)] = np.mean(np.abs(dists_x[2000:] - dists_y[2000:]))
        
        avg_data["dists_u"][(press, temp)] = np.mean(dists_u[2000:])
        avg_data["dists_d"][(press, temp)] = np.mean(dists_d[2000:])
        avg_data["dists_O1"][(press, temp)] = np.mean(dists_O1[2000:])
        avg_data["dists_O3"][(press, temp)] = np.mean(dists_O3[2000:])
        avg_data["diff_dists_colinear"][(press, temp)] = np.mean(np.abs(dists_u[2000:] - dists_d[2000:]))
        avg_data["diff_dists_colinear2"][(press, temp)] = np.abs(np.mean(dists_u[2000:]) - np.mean(dists_d[2000:]))
        
        avg_data["angles_intra"][(press, temp)] = 180-np.mean(np.abs(180-angles_intra[2000:]))
        avg_data["angles_intra2"][(press, temp)] = 180-np.mean(np.abs(180-angles_intra[2000:]))
        avg_data["angles_inter"][(press, temp)] = 180-np.mean(np.abs(180-angles_inter[2000:]))
        
        avg_data["ONiO_x"][(press, temp)] = np.mean(ONiO_x[2000:])
        avg_data["ONiO_y"][(press, temp)] = np.mean(ONiO_y[2000:])
        avg_data["ONiO_diff"][(press, temp)] = np.mean(ONiO_diff[2000:])
        
        avg_data["signed_angles_inter"][(press, temp)] = np.mean(angles_inter[2000:])
        avg_data["STD_inter"][(press, temp)] = np.std(angles_inter[2000:])
        
        avg_data["a"][(press, temp)] = np.mean(a[2000:])
        avg_data["b"][(press, temp)] = np.mean(b[2000:])
        avg_data["diff_a_b"][(press, temp)] = np.mean(np.abs(a[2000:] - b[2000:]))
        avg_data["diff_a_b_2"][(press, temp)] = np.abs(np.mean(a[2000:]) - np.mean(b[2000:]))
        avg_data["c"][(press, temp)] = np.mean(c[2000:])
        avg_data["alpha"][(press, temp)] = np.mean(alpha[2000:])
        avg_data["beta"][(press, temp)] = np.mean(beta[2000:])
        avg_data["gamma"][(press, temp)] = np.mean(gamma[2000:])
        avg_data["temperature"][(press, temp)] = np.mean(temperatures[2000:])
        avg_data["pressure"][(press, temp)] = np.mean(pressures[2000:])
        avg_data["free_energy"][(press, temp)] = np.mean(free_energies[2000:])
        avg_data["volume"][(press, temp)] = np.mean(volumes[2000:])
        avg_data["kinetic_energy"][(press, temp)] = np.mean(ekin_values[2000:])
        avg_data["total_energy"][(press, temp)] = np.mean(etotal_values[2000:])
        
        avg_data["timestep"][(press, temp)] = len(timesteps)
        
with open("AIMD_avg_data.txt", "w") as file:
    file.write("parameter\tpressure\ttemperature\tvalue\n")
    for parameter, sub_dict in avg_data.items():
        for (pressure, temperature), value in sub_dict.items():
            file.write(f"{parameter}\t{pressure:.6f}\t{temperature:.6f}\t{value:.6f}\n")

def plot_color_map(data_dict, title, xlabel, ylabel, quantity_label):

    temperatures = sorted(set(k[1] for k in data_dict.keys()))
    pressures = sorted(set(k[0] for k in data_dict.keys()))
    
    # Add negative values for alignment if necessary
    if min(temperatures) > 0:
        temperatures.insert(0, -0.5)
    if min(pressures) > 0:
        pressures.insert(0, -0.5)

    data_grid = np.zeros((len(temperatures), len(pressures)))

    # Fill the data grid with values from the data_dict
    for i, temp in enumerate(temperatures):
        for j, press in enumerate(pressures):
            data_grid[i, j] = data_dict.get((press, temp), np.nan)

    if len(pressures) > 0:
        plt.figure(figsize=(10, 6))
        
        # Adjust the extent so the cells are centered on the correct coordinates
        plt.imshow(data_grid, aspect='auto', origin='lower',
                   extent=[min(pressures) - 0.5, max(pressures) + 0.5, 
                           min(temperatures) - 0.5, max(temperatures) + 0.5],
                   cmap='viridis')

        plt.colorbar(label=quantity_label)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        
        # Set ticks at the actual pressure and temperature points
        plt.xticks(ticks=pressures, labels=pressures)
        plt.yticks(ticks=temperatures, labels=temperatures)
        
        plt.title("Average " + title + " after 2 ps")
        plt.tight_layout()
        plt.savefig("averages/" + str(title))
        plt.close()


def plot_temps(data_dict, title, xlabel, ylabel, quantity_label):

    temperatures = sorted(set(k[1] for k in data_dict.keys()))
    pressures = sorted(set(k[0] for k in data_dict.keys()))

    plt.figure(figsize=(10, 6))
    for temp in temperatures:
        y_values = [data_dict.get((press, temp), np.nan) for press in pressures]
        filtered_pressures = [press for press, y in zip(pressures, y_values) if not np.isnan(y)]
        filtered_y_values = [y for y in y_values if not np.isnan(y)]
        if filtered_y_values:
            plt.plot(filtered_pressures, filtered_y_values, label=f"T = {temp} K")
    plt.xlabel(xlabel)
    plt.ylabel(quantity_label)
    plt.title(f"{title}")
    plt.legend(title="Temperature")
    plt.tight_layout()
    plt.savefig("constemps/"+str(title))
    plt.close()

def plot_press(data_dict, title, xlabel, ylabel, quantity_label):
    pressures = sorted(set(k[0] for k in data_dict.keys()))
    temperatures = sorted(set(k[1] for k in data_dict.keys()))

    plt.figure(figsize=(10, 6))
    for press in pressures:
        y_values = [data_dict.get((press, temp), np.nan) for temp in temperatures]
        filtered_temperatures = [temp for temp, y in zip(temperatures, y_values) if not np.isnan(y)]
        filtered_y_values = [y for y in y_values if not np.isnan(y)]
        if filtered_y_values:
            plt.plot(filtered_temperatures, filtered_y_values, label=f"P = {press} atm")

    plt.xlabel(xlabel)
    plt.ylabel(quantity_label)
    plt.title(f"{title}")
    plt.legend(title="Pressure")
    plt.tight_layout()
    plt.savefig("conspres/" + str(title))
    plt.close()

for key, data_dict in avg_data.items():

    title = f"{key.replace('_', ' ').title()}"
    xlabel = "Pressure (GPa)"  # Pressure now on the x-axis
    ylabel = "Temperature (K)"  # Temperature now on the y-axis
    quantity_label = key.replace('_', ' ').title()

    plot_color_map(data_dict, title, xlabel, ylabel, quantity_label)
    plot_temps(data_dict, title, xlabel, ylabel, quantity_label)
    plot_press(data_dict, title, xlabel, ylabel, quantity_label)


dtype = [('parameter', 'U20'), ('pressure', 'f8'), ('temperature', 'f8'), ('value', 'f8')]
loaded_data = np.loadtxt("AIMD_avg_data.txt", dtype=dtype, delimiter="\t", skiprows=1)

avg_data = {}
for entry in loaded_data:
    parameter = entry['parameter']
    pressure = entry['pressure']
    temperature = entry['temperature']
    value = entry['value']
    
    if parameter not in avg_data:
        avg_data[parameter] = {}
    avg_data[parameter][(pressure, temperature)] = value


def plot_color_map(data_dict, title, xlabel, ylabel, quantity_label):

    temperatures = sorted(set(k[1] for k in data_dict.keys()))
    pressures = sorted(set(k[0] for k in data_dict.keys()))
    
    if min(temperatures) > 0:
        temperatures.insert(0, -0.5)
    if min(pressures) > 0:
        pressures.insert(0, -0.5)

    data_grid = np.zeros((len(temperatures), len(pressures)))

    for i, temp in enumerate(temperatures):
        for j, press in enumerate(pressures):
            data_grid[i, j] = data_dict.get((press, temp), np.nan)

    if len(pressures) > 0:
        plt.figure(figsize=(10, 6))
        
        plt.imshow(data_grid, aspect='auto', origin='lower',
                   extent=[min(pressures) - 0.5, max(pressures) + 0.5, 
                           min(temperatures) - 0.5, max(temperatures) + 0.5],
                   cmap='viridis')

        plt.colorbar(label=quantity_label)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.xticks(ticks=pressures, labels=pressures)
        plt.yticks(ticks=temperatures, labels=temperatures)
        plt.tight_layout()
        plt.savefig("AIMD_" + str(title))
        plt.close()

from matplotlib.ticker import MaxNLocator

def plot_color_map_smooth(data_dict, title, xlabel, ylabel, quantity_label):

    # Extract data points and values
    points = np.array(list(data_dict.keys()))
    pressures = points[:, 0]
    temperatures = points[:, 1]
    values = np.array(list(data_dict.values()))

    # Create a finer grid for interpolation
    pressure_fine = np.linspace(min(pressures), max(pressures), 200)
    temperature_fine = np.linspace(min(temperatures), max(temperatures), 200)
    grid_x, grid_y = np.meshgrid(pressure_fine, temperature_fine)
    
    # Interpolate the values onto the grid
    grid_z = np.clip(griddata(points, values, (grid_x, grid_y), method='cubic'), 0, None)

    # Plot the contour map
    plt.figure(figsize=(10, 6))
    contour = plt.contourf(grid_x, grid_y, grid_z, levels=20, cmap='rainbow')
    
    # Add colorbar
    cbar = plt.colorbar(contour)
    cbar.ax.set_title(quantity_label, pad=15, loc='center')
    
    # Reduce the number of ticks on the colorbar
    cbar.locator = MaxNLocator(nbins=5)  # Control number of colorbar ticks
    cbar.update_ticks()  # Apply changes

    # Labels and layout
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.tight_layout(pad=0.1)

    # Sanitize the title for the filename
    sanitized_title = title.replace("*", "_").replace(":", "_").replace("/", "_").replace("\\", "_")

    # Save the plot
    plt.savefig(f"AIMD_{sanitized_title}.png", dpi=500)
    plt.close()

keys = ["angles_intra", "angles_intra2", "angles_inter", "volume", "pressure", "temperature", "dists_Ni", "timestep", "diff_a_b", "diff_a_b_2", "a"]
titles = [r"$\theta_2$ [°]", r"$\theta_3$ [°]", r"$\theta_1$ [°]", r"Cell Volume [Å³]", "Pressure [GPa]", "Temperature [K]", "Ni-Ni distance [Å]", "Time [fs]", r"$|a-b|$ [Å]", "test2", "a"]

product_key = "angles_intra*angles_inter"
product_title = r"$\theta_2$ × $\theta_1$ [°²]"
data_dict_product = {key: (avg_data["angles_intra"][key]) * avg_data["angles_inter"][key] for key in avg_data["angles_intra"].keys()}

keys.append(product_key)
titles.append(product_title)

for j in range(len(keys)):
    key, title = keys[j], titles[j]

    if key == product_key:
        data_dict = data_dict_product
    else:
        data_dict = avg_data[key]

    xlabel = "Pressure [GPa]"
    ylabel = "Temperature [K]"
    quantity_label = title

    plot_color_map_smooth(data_dict, key, xlabel, ylabel, quantity_label)