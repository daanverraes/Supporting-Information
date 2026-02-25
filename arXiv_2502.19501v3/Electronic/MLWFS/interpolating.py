import re
from ase import io
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import colormaps
from matplotlib import rcParams

rcParams.update({'font.size': 20})
rcParams['font.family'] = 'Nimbus Roman'
rcParams["mathtext.fontset"] = "cm"

bands = True
project = True
windows = True
MLWF = True

def extract_fermi_energy(file_path):
    with open(file_path, "r") as f:
        for line in f:
            if "E-fermi" in line:
                return float(line.split()[2])
    raise ValueError(f"Could not find Fermi energy in {file_path}")

def extract_incar_data(file_path):
    N_bands = None
    N_kseg = None
    K_points = []
    dis_win_min = 0
    dis_win_max = 0
    dis_froz_min = 0
    dis_froz_max = 0
    kpoint_path = []
    NUM_WANN = None
    in_kpath = False

    with open(file_path, "r") as f:
        lines = f.read().splitlines()

        for line in lines:

            if "NBANDS" in line:
                N_bands = int(re.sub(r'\D', '', line.split("=")[1].strip()))

            if "dis_win_min" in line:
                dis_win_min = float(line.split("=")[1].strip())
            if "dis_win_max" in line:
                dis_win_max = float(line.split("=")[1].strip())

            if "dis_froz_min" in line:
                dis_froz_min = float(line.split("=")[1].strip())
            if "dis_froz_max" in line:
                dis_froz_max = float(line.split("=")[1].strip())
            
            if "NUM_WANN" in line:
                NUM_WANN = int(line.split("=")[1].split("#")[0].strip())


            if "begin kpoint_path" in line:
                in_kpath = True
                continue
            if "end kpoint_path" in line:
                in_kpath = False
                continue
            if in_kpath and re.match(r'^\s*[A-Za-z]+\s+[\d\.\-]+\s+[\d\.\-]+\s+[\d\.\-]+\s+[A-Za-z]+\s+[\d\.\-]+\s+[\d\.\-]+\s+[\d\.\-]+', line):
                kpoint_path.append(line.strip().split())

            if kpoint_path:
                raw_points = []
                for kp in kpoint_path:
                    raw_points.append(kp[0])
                    raw_points.append(kp[4])

                K_points = []
                for point in raw_points:
                    if not K_points or K_points[-1] != point:
                        K_points.append(point)

                K_points = [r'$\Gamma$' if p == 'G' else p for p in K_points]

            if "bands_num_points" in line:
                N_kseg = int(line.strip().split()[-1])

    return {
        'N_bands': N_bands,
        'N_kseg': N_kseg,
        'K_points': K_points,
        'dis_win_min': dis_win_min,
        'dis_win_max': dis_win_max,
        'dis_froz_min': dis_froz_min,
        'dis_froz_max': dis_froz_max,
        'NUM_WANN': NUM_WANN
    }

incar_data = extract_incar_data("INCAR")
N_bands = incar_data['N_bands']
N_kseg = incar_data['N_kseg']
K_points = incar_data['K_points']
dis_win_min = incar_data['dis_win_min']
dis_win_max = incar_data['dis_win_max']
dis_froz_min = incar_data['dis_froz_min']
dis_froz_max = incar_data['dis_froz_max']
N_w = incar_data['NUM_WANN']
N_orbs = 16

if bands: WB = open('../BANDS/OUTCAR', "r")

if bands: E_F_Bands = extract_fermi_energy("../BANDS/OUTCAR")
if MLWF: E_F_MLWF = extract_fermi_energy("../LWL/OUTCAR") 

WB = WB.read()
WB = np.array([x for x in WB.split(" ") if x != '' and x != '\n'])
k_start = np.where(WB == 'occupation')[0] + 2

atoms = io.read('POSCAR')
reciprocal_lattice = atoms.cell.reciprocal()
b1, b2, b3 = reciprocal_lattice
reciprocal_matrix = np.array([b1, b2, b3]).T

dk = []
E = np.zeros((len(k_start), N_bands))

for i in range(1, len(k_start)):
    dx = float(WB[k_start[i]-9]) - float(WB[k_start[i-1]-9])
    dy = float(WB[k_start[i]-8]) - float(WB[k_start[i-1]-8])
    dz = float(WB[k_start[i]-7]) - float(WB[k_start[i-1]-7])

    k_diff = np.array([dx, dy, dz])
    dk.append(np.linalg.norm(np.dot(reciprocal_matrix, k_diff)))

    for j in range(N_bands):
        E[i][j] = WB[k_start[i] + 3 * j]

dk = np.array(dk)/np.sum(dk)
E = np.transpose(E)

kk = np.zeros(len(k_start))
kk[1:] = np.cumsum(dk[:len(k_start)-1])
K_values = [kk[0]] + [kk[K*N_kseg - 1] for K in range(1, len(K_points))]

if project:

    with open('../BANDS/PROCAR', "r") as f:
        lines = f.readlines()

    for line in lines:
        if "k-points" in line:
            parts = line.split()
            N_k = int(parts[3])
            N_b = int(parts[7])
            N_i = int(parts[11])
            break

    Data = np.zeros((N_b, N_k, N_i, 5))

    k_idx = -1
    b_idx = -1
    line_idx = 0

    while line_idx < len(lines):
        line = lines[line_idx].strip()

        if line.startswith("k-point"):
            k_idx += 1
            b_idx = -1
            line_idx += 1
            continue

        if line.startswith("band"):
            b_idx += 1
            while not lines[line_idx].strip().startswith("ion"):
                line_idx += 1
            line_idx += 1 

            for i in range(N_i):
                parts = lines[line_idx].split()

                Data[b_idx, k_idx, i, 0] = float(parts[-1])
                Data[b_idx, k_idx, i, 1] = float(parts[5]) + float(parts[6]) + float(parts[7]) + float(parts[8]) + float(parts[9])
                Data[b_idx, k_idx, i, 2] = float(parts[2]) + float(parts[3]) + float(parts[4])
                Data[b_idx, k_idx, i, 3] = float(parts[9])
                Data[b_idx, k_idx, i, 4] = float(parts[7])

                line_idx += 1

            line_idx += 1  
            continue

        line_idx += 1

if MLWF:
    WBI = open('wannier90_band.dat', "r")
    INF = open('wannier90_band.labelinfo.dat', "r")

    WBI = WBI.read()
    WBI = np.array([x for x in WBI.split(" ")[2:] if x != '' and x != '\n'])

    INF = INF.read()
    INF = " ".join(INF.split())
    INF= np.array([x for x in INF.split(" ")])

    k_indices = np.array([int(INF[6*i+1]) for i in range(len(INF)//6)])
    k_coords = np.array([float(INF[6*i+2]) for i in range(len(INF)//6)])

    segment_lengths = np.diff(k_coords)
    total_length = np.sum(segment_lengths)
    num_segments = len(segment_lengths)
    rescaled_coords = [k_coords[0]] 

    for i in range(num_segments):
        frac_tot = kk[N_kseg*(i+1)-1]-kk[N_kseg*(i)]
        segment_points = np.linspace(0, total_length*frac_tot, k_indices[i+1] - k_indices[i])
        rescaled_coords.extend(rescaled_coords[-1] + segment_points)

    rescaled_coords = np.array(rescaled_coords)

    k_full = np.array([float(WBI[2*i]) for i in range(len(WBI)//2)])
    E_full = np.array([float(WBI[2*i+1][:-1]) for i in range(len(WBI)//2)])
    Bands = np.array_split(E_full, len(k_full))
    assert len(k_full) == len(Bands)

    Bands_full = [float(i[0]) for i in Bands]

    Bands = np.array([Bands_full[i*len(Bands_full)//N_w:i*len(Bands_full)//N_w+len(Bands_full)//N_w] for i in range(N_w)])
    k = np.linspace(0, K_values[-1], len(rescaled_coords))


# Plotting
# --------

E -= E_F_Bands

fig, axs = plt.subplots(figsize=(8, 6))
norm = plt.Normalize(0, 1)
original_cmap = colormaps.get_cmap('jet')
truncated_cmap = LinearSegmentedColormap.from_list('truncated_jet', original_cmap(np.linspace(0.1, 0.9, 256)))

if project:
    all_dot_sizes = Data[:, :, :3, 3].sum(axis=2) 
    max_size = np.max(all_dot_sizes) 

if bands:
    for i in range(len(E)):
        seg = np.array([[[kk[j], E[i][j]], [kk[min(j + 1, len(kk) - 1)], E[i][min(j + 1, len(kk) - 1)]]] for j in range(1,len(kk))])
        
        if project:
            upper_contrib = Data[i, :, 6:11, 3].sum(axis=1)
            lower_contrib = Data[i, :, 6:11, 4].sum(axis=1)
            relative_weight = upper_contrib / (upper_contrib + lower_contrib + 1e-8)            
            total_contrib = upper_contrib + lower_contrib
            whiteness = 1 - total_contrib

            dot_sizes = Data[i, :, :6, 1].sum(axis=1)
            
            if max_size > 0:
                absolute_sizes = 20 * dot_sizes 
            else:
                absolute_sizes = np.zeros_like(dot_sizes)
        
            axs.scatter(kk, E[i], s=absolute_sizes, color='darkviolet', alpha=0.7, edgecolors='none')

        lc = LineCollection(seg, cmap=truncated_cmap)  
        
        if project:
            lc.set_array(relative_weight) 
            lc.set_alpha(1 - whiteness) 
        
        lc.set_linewidth(3.0)
        axs.add_collection(lc)
        line = lc

        plain_lc = LineCollection(seg, colors='black')
        plain_lc.set_linewidth(0.5)
        axs.add_collection(plain_lc)

if MLWF:
    for i in range(len(Bands)):
        
        if i == 0:
            plt.plot(k[:], Bands[i] - E_F_MLWF, color='black', linewidth=1.0, linestyle='--', label='MLWF')
        else:
            plt.plot(k, Bands[i] - E_F_MLWF, color='black', linewidth=1.0, linestyle='--')
            
if bands and project:
    cbar = fig.colorbar(line, ax=axs, pad=0.08)
    cbar.set_ticks([line.norm.vmin, line.norm.vmax]) 
    cbar.set_ticklabels(['', ''])
    cbar.ax.annotate(r'$d_{z^2}$', xy=(0.6, -0.08), xycoords='axes fraction', ha='center', va='bottom')
    cbar.ax.annotate(r'$d_{x^2-y^2}$', xy=(0.6, 1.08), xycoords='axes fraction', ha='center', va='top')

plt.ylabel(r'$E-E_{F}$ [eV]')
plt.xlabel('K-points')
plt.ylim(-3.2, 10.0)
plt.axhline(0, linestyle=(0, (5, 5)), linewidth=1, color='g', alpha=1)

if windows:
    plt.axhline(dis_win_max - E_F_Bands, linestyle=(0, (5, 5)), linewidth=1, color='r', alpha=1)
    plt.axhline(dis_win_min - E_F_Bands, linestyle=(0, (5, 5)), linewidth=1, color='r', alpha=1)

    plt.axhline(dis_froz_max - E_F_Bands, linestyle=(0, (5, 5)), linewidth=1, color='b', alpha=1)
    plt.axhline(dis_froz_min - E_F_Bands, linestyle=(0, (5, 5)), linewidth=1, color='b', alpha=1)

plt.yticks(fontsize=18)
plt.xticks(K_values, K_points)
plt.grid(axis='x')
axs.set_xlim(kk[0], kk[-1])
plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)
plt.tight_layout(pad=0)
plt.savefig("bands.png", bbox_inches='tight', dpi=300)