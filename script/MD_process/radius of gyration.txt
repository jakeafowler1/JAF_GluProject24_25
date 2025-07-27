import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import MDAnalysis as mda

# === USER SETTINGS ===
simulation_dirs = [
    "r375k_pro9_trajout",
    "r375l_v2_pro9_trajout",
    "rw_lf_pro8_trajout",
    "w376a_v2_pro9_trajout",
    "wt_pro8_trajout"
]

# Desired label for legend (folder name → label)
label_map = {
    "r375k_pro9_trajout": "R375K",
    "r375l_v2_pro9_trajout": "R375L",
    "rw_lf_pro8_trajout": "RW_LF",
    "w376a_v2_pro9_trajout": "W376A",
    "wt_pro8_trajout": "WT"
}

# Custom color map (label → color)
color_map = {
    "R375K": "darkblue",
    "R375L": "darkgreen",
    "RW_LF": "lightsteelblue",
    "W376A": "orange",
    "WT": "black"
}

# Plot styling
y_min = 35
y_max = 45
label_fontsize = 18
tick_fontsize = 16
title_fontsize = 20
legend_fontsize = 14

# === FUNCTION: Extract total simulation time from .mdinfo ===
def read_simulation_time(mdinfo_path):
    try:
        with open(mdinfo_path, 'r') as file:
            lines = file.readlines()
        for line in lines:
            if 'TIME(PS)' in line:
                match = re.search(r"TIME\(PS\)\s*=\s*([\d\.]+)", line)
                if match:
                    return float(match.group(1))
        raise ValueError(f"TIME(PS) not found in {mdinfo_path}")
    except Exception as e:
        print(f"Error reading {mdinfo_path}: {e}")
        return None

# === FUNCTION: Compute Rg and save CSV ===
def compute_rg_csv(sim_dir, total_time_ps):
    topology = os.path.join(sim_dir, "traj.pdb")
    trajectory = os.path.join(sim_dir, "traj.dcd")
    output_csv = os.path.join(sim_dir, "radius_of_gyration_summary.csv")

    if not os.path.exists(topology) or not os.path.exists(trajectory):
        print(f"⚠️  Missing traj files in {sim_dir}")
        return None

    try:
        u = mda.Universe(topology, trajectory)
        whole_protein = u.select_atoms("protein")

        n_frames = len(u.trajectory)
        time_ns_array = np.linspace(0, total_time_ps / 1000, n_frames)

        rg_data = {
            'Time_ns': [],
            'Whole_Rg': []
        }

        for i, ts in enumerate(u.trajectory):
            rg_data['Time_ns'].append(time_ns_array[i])
            rg_data['Whole_Rg'].append(whole_protein.radius_of_gyration())

        df = pd.DataFrame(rg_data)
        df.to_csv(output_csv, index=False)
        return df

    except Exception as e:
        print(f"❌ Error processing {sim_dir}: {e}")
        return None

# === PROCESS ALL SIMULATIONS ===
max_times_ns = []
sim_data = []

for sim_dir in simulation_dirs:
    label = label_map.get(sim_dir, os.path.basename(sim_dir).upper())

    mdinfo_path = os.path.join(sim_dir, "pro9.mdinfo")
    if not os.path.exists(mdinfo_path):
        mdinfo_path = os.path.join(sim_dir, "pro8.mdinfo")
    if not os.path.exists(mdinfo_path):
        print(f"⚠️  Skipping {label}: no mdinfo found.")
        continue

    total_time_ps = read_simulation_time(mdinfo_path)
    if total_time_ps is None:
        continue
    total_time_ns = total_time_ps / 1000
    max_times_ns.append(total_time_ns)

    # Generate CSV and collect data
    df = compute_rg_csv(sim_dir, total_time_ps)
    if df is not None:
        sim_data.append((df, label, total_time_ns))

# === Determine shortest simulation time
if not max_times_ns:
    raise RuntimeError("❌ No valid simulations found.")

global_max_time_ns = min(max_times_ns)
print(f"✅ Using shortest simulation time: {global_max_time_ns:.2f} ns")

# === PLOT COMBINED ===
plt.figure(figsize=(12, 6))

for df, label, total_time_ns in sim_data:
    color = color_map.get(label, "grey")
    clipped_df = df[df['Time_ns'] <= global_max_time_ns]
    plt.plot(clipped_df['Time_ns'], clipped_df['Whole_Rg'], label=label, linewidth=2, color=color)

# === Plot styling ===
plt.xlabel("Time (ns)", fontsize=label_fontsize)
plt.ylabel("Radius of Gyration (Å)", fontsize=label_fontsize)
plt.title("Radius of Gyration (Whole Tetramer) Across Simulations", fontsize=title_fontsize, pad=15)
plt.ylim(y_min, y_max)
plt.xlim(90, global_max_time_ns)
plt.xticks(fontsize=tick_fontsize)
plt.yticks(fontsize=tick_fontsize)
plt.legend(fontsize=legend_fontsize, bbox_to_anchor=(1.01, 1), loc='upper left')
plt.tight_layout(rect=[0, 0, 0.85, 1])

# === Save plot ===
plt.savefig("combined_radius_of_gyration_plot.png", dpi=400)
plt.savefig("combined_radius_of_gyration_plot.svg", dpi=400)
plt.show()
