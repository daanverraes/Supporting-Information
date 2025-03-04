import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib import rcParams

jlred = tuple(np.array([0.796, 0.235, 0.2]))
jlblue = tuple(np.array([0.251, 0.388, 0.847]))
jlgreen = tuple(np.array([0.133, 0.545, 0.133]))
jlpurple = tuple(np.array([0.584, 0.345, 0.698]))
jlgrey = tuple(np.array([0.4, 0.4, 0.4]))

rcParams.update({'font.size': 26}) 
rcParams['font.family'] = 'Times New Roman'
rcParams["mathtext.fontset"] = "cm"

def read_xdatcar(filepath):
    with open(filepath, 'r') as f:
        lines = f.readlines()

    lattice_vectors = []
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
            i += total_atoms + 1
        else:
            i += 1

    return lattice_vectors, np.array(timesteps)

root_directory = "."
output_file = "lattice_params.txt"

with open(output_file, "w") as f_out:
    f_out.write("Temp(K) Pressure(GPa) Timestep a(Å) b(Å)\n") 

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

            xdatcar_file = os.path.join(folder_path, "XDATCAR")
            if not os.path.isfile(xdatcar_file) or os.path.getsize(xdatcar_file) == 0:
                continue

            lattice_vectors, timesteps = read_xdatcar(xdatcar_file)

            a = np.array([np.linalg.norm(lattice_vectors[t][0]) for t in range(len(timesteps))])
            b = np.array([np.linalg.norm(lattice_vectors[t][1]) for t in range(len(timesteps))])

            for t, a_val, b_val in zip(timesteps, a, b):
                f_out.write(f"{temp} {press} {t} {a_val:.6f} {b_val:.6f}\n")

print(f"Saved lattice parameters to {output_file}")

data = np.loadtxt(output_file, skiprows=1)

temps = data[:, 0]
pressures = data[:, 1]
timesteps = data[:, 2]
a_values = data[:, 3]
b_values = data[:, 4]

# Ensure output directory exists
output_dir = "."
os.makedirs(output_dir, exist_ok=True)

conditions = [(60, 75), (60, 0), (10, 75), (10, 0)]
colors = [jlred, jlblue, jlgreen, (255/255, 165/255, 0)]  # Red, Blue, Green, Orange

plt.figure(figsize=(12, 8))

for (temp, press), color in zip(conditions, colors):
    mask = (temps == temp) & (pressures == press)
    a_minus_b = a_values[mask] - b_values[mask]

    plt.hist(
        a_minus_b[2000:], bins=np.linspace(-0.055, 0.155, 51), 
        color=color, alpha=0.6, edgecolor='black', label=f"$T = {temp}$ K, $P = {press}$ GPa"
    )

plt.axvline(x=0, color='black', linestyle='--', linewidth=2)  # Vertical line at x=0
plt.xlabel(r"$a - b$ [Å]")
plt.legend(fontsize=20)
plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()
output_path = os.path.join(output_dir, "histograms.png")
plt.savefig(output_path, dpi=300)
print(f"Saved: {output_path}")

