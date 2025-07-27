import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ---- INPUT FILES ----
topology_file = "r375k_pro9_trajout/traj.pdb"
trajectory_file = "r375k_pro9_trajout/traj.dcd"

# ---- HBOND PAIRS ----
hbonds_to_track = [
    ("K375(A)] to D71(A)",  "resid 333", "resid 29"),
    ("K375(A) to E371(A)", "resid 333", "resid 329"),
    ("K375(A) to E371(C)", "resid 333", "resid 1003"),
    ("K375(A) to E374(C)", "resid 333", "resid 1006"),
]

# ---- HBOND CUTOFFS ----
distance_cutoff = 3.5  # Angstroms
angle_cutoff = 150     # Degrees

# ---- LOAD UNIVERSE ----
u = mda.Universe(topology_file, trajectory_file)
n_frames = len(u.trajectory)

print(f"Trajectory contains {n_frames} frames.")
print(f"Hydrogens present?: {len(u.select_atoms('name H*')) > 0}")

# ---- HBOND ANALYSIS ----
occupancy_data = []

for label, donor_sel, acceptor_sel in hbonds_to_track:
    print(f"\n🔎 Checking: {label}")
    donor_group = u.select_atoms(donor_sel)
    acceptor_group = u.select_atoms(acceptor_sel)
    print(f"  Donor selection: {len(donor_group)} atoms")
    print(f"  Acceptor selection: {len(acceptor_group)} atoms")

    if len(donor_group) == 0 or len(acceptor_group) == 0:
        print("  ⚠️ Empty selection — skipping")
        occupancy_data.append((label, 0.0))
        continue

    hba = HydrogenBondAnalysis(
        universe=u,
        donors_sel=donor_sel,
        acceptors_sel=acceptor_sel
    )
    hba.distance = distance_cutoff
    hba.angle = angle_cutoff
    hba.hydrogens_sel = "name H*"
    hba.run()

    occupancy = np.zeros(n_frames, dtype=int)
    if hba.hbonds is not None:
        frames_with_hbonds = hba.hbonds[:, 0].astype(int)
        occupancy[frames_with_hbonds] = 1

    percent_occupancy = 100 * np.sum(occupancy) / n_frames
    print(f"  ✅ % Occupancy: {percent_occupancy:.2f}%")
    occupancy_data.append((label, percent_occupancy))

# ---- PREP DATA ----
labels = [x[0] for x in occupancy_data]
values = [x[1] for x in occupancy_data]
df = pd.DataFrame(values, index=labels, columns=["% Occupancy"])

# ---- TUNABLE BAR CHART ----
# --- TUNING PARAMETERS ---
fig_width = 6            # Figure width (inches)
fig_height_per_bar = 0.6 # Vertical size per bar
label_fontsize = 12
annotation_fontsize = 14
x_limit = 100            # Max X-axis (%)
text_offset = 1.5        # Text distance from bar

# --- CREATE PLOT ---
plt.figure(figsize=(fig_width, fig_height_per_bar * len(labels)))
bars = plt.barh(labels, values, color="darkblue", edgecolor="black")

plt.xlabel("Hydrogen Bond % Occupancy", fontsize=label_fontsize)
plt.title("Hydrogen Bond Occupancy Across Trajectory", fontsize=label_fontsize + 1)
plt.xlim(0, x_limit)

# Annotate bars
for bar, val in zip(bars, values):
    plt.text(val + text_offset, bar.get_y() + bar.get_height() / 2,
             f"{val:.1f}%", va='center', fontsize=annotation_fontsize)

plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

plt.tight_layout()
plt.savefig('KW_H_Bond_KW.svg', dpi=400)
plt.show()
