import pandas as pd
import os
from Bio import SeqIO, AlignIO
import matplotlib.pyplot as plt

# Define color scheme
clade_colors = {
    "A1": "#888888",   # grey
    "A2": "#000000",   # black
    "B":  "#80597B",   # muted purple
    "C":  "#7F1C1C",   # deep red
    "D":  "#A07E5F",   # brown-tan
    "E":  "#CABD9F",   # beige
    "F":  "#C4D94B",   # lime green
    "G":  "#7DC425",   # vivid green
    "H":  "#5C8124",   # olive green
    "I":  "#2C6A6A",   # teal
    "J":  "#255986",   # dark blue
    "K":  "#4746C7",   # medium blue
    "L":  "#6352C5",   # purple-blue
}

# --- Set your directory path ---
input_dir = "GluR0_tree_clades/clades"

# Load reference sequence from ssGluR0_aln.txt
ref_file = os.path.join(input_dir, "ssGluR0_aln.txt")
ref_record = list(SeqIO.parse(ref_file, "fasta"))[0]
ref_seq = str(ref_record.seq)

# Get positions in the alignment where reference has an amino acid
occupied_positions = [i for i, aa in enumerate(ref_seq) if aa != '-']
ref_aa_labels = [ref_seq[i] for i in occupied_positions]
ref_position_map = list(range(1, len(occupied_positions) + 1))  # 1-indexed

# Build table mapping reference position → alignment column number
column_mapping_df = pd.DataFrame({
    "Ref_AA_Position": ref_position_map,
    "Alignment_Column": occupied_positions,
    "Reference_AA": ref_aa_labels
})
column_mapping_df.to_csv("reference_position_mapping.csv", index=False)

# Get all clade alignment files (excluding the reference file)
clade_files = sorted(f for f in os.listdir(input_dir)
                     if f.endswith(".txt") and not f.startswith("ssGluR0"))
clade_labels = [f.replace("Clade ", "").replace(".txt", "") for f in clade_files]

# Store conservation data per clade
clade_conservation = {}
for file, label in zip(clade_files, clade_labels):
    aln = AlignIO.read(os.path.join(input_dir, file), "fasta")
    sequences = [str(rec.seq) for rec in aln]
    
    conservation = []
    for i in occupied_positions:
        ref_aa = ref_seq[i]
        matches = sum(1 for seq in sequences if seq[i] == ref_aa)
        percent = (matches / len(sequences)) * 100
        conservation.append(percent)
    
    clade_conservation[label] = conservation

# --- Create CSV output with all conservation data ---
# Build DataFrame with conservation data
conservation_df = pd.DataFrame({
    "Ref_AA_Position": ref_position_map,
    "Reference_AA": ref_aa_labels
})

# Add conservation percentages for each clade
for clade, conservation in clade_conservation.items():
    conservation_df[f"Clade_{clade}_Conservation_%"] = conservation

# Save the conservation data to CSV
conservation_df.to_csv("clade_conservation_data.csv", index=False)
print(f"Conservation data saved to 'clade_conservation_data.csv'")
print(f"Shape: {conservation_df.shape[0]} positions × {conservation_df.shape[1]} columns")

# --- Plotting ---
# Font size settings - adjust these values as needed
main_title_fontsize = 20       # Main plot title
subplot_title_fontsize = 20     # Individual subplot titles
axis_label_fontsize = 16        # X and Y axis labels
tick_label_fontsize = 16        # Tick labels (numbers/letters on axes)
legend_fontsize = 20            # Legend text
legend_title_fontsize = 20      # Legend title

# Line appearance settings
line_width = 1.8               # Line thickness (default is usually 1.0)

# Legend appearance settings
legend_marker_size = 6          # Size of legend markers (default is usually 6)
legend_marker_width = 8         # Width of legend line markers (default is usually 2)

num_positions = len(ref_position_map)
chunk_size = num_positions // 3 + (num_positions % 3 > 0)
fig, axs = plt.subplots(3, 1, figsize=(24, 12), sharey=True)
fig.subplots_adjust(hspace=0.4)

for i in range(3):
    start = i * chunk_size
    end = min(start + chunk_size, num_positions)
    x_vals = ref_position_map[start:end]
    x_labels = ref_aa_labels[start:end]
    
    for clade, conservation in clade_conservation.items():
        # Use the custom color if available, otherwise use default matplotlib color
        color = clade_colors.get(clade, None)
        axs[i].plot(x_vals, conservation[start:end], label=clade, color=color, linewidth=line_width)
    
    axs[i].set_title(f"Positions {x_vals[0]}–{x_vals[-1]}", fontsize=subplot_title_fontsize)
    axs[i].set_ylabel("Conservation % (match to reference)", fontsize=axis_label_fontsize)
    axs[i].set_xticks(x_vals)
    axs[i].set_xticklabels(x_labels, rotation=0, fontsize=tick_label_fontsize)
    
    # Set Y-axis tick label font size
    axs[i].tick_params(axis='y', labelsize=tick_label_fontsize)
    
    # Remove padding at the ends so lines intersect y-axis
    axs[i].set_xlim(x_vals[0], x_vals[-1])

axs[-1].set_xlabel("Reference sequence (1–222, amino acid at each position)", fontsize=axis_label_fontsize)
legend = axs[0].legend(title="Clade", bbox_to_anchor=(1.05, 1), loc='upper left', 
                      fontsize=legend_fontsize, title_fontsize=legend_title_fontsize,
                      markerscale=legend_marker_size, markerfirst=True)
# Make legend markers thicker/more visible
for line in legend.get_lines():
    line.set_linewidth(legend_marker_width)
fig.suptitle("Per-position Conservation by Clade (Reference-matched Positions Only)", fontsize=main_title_fontsize)
#plt.tight_layout(rect=[0, 0, 1, 0.97])
plt.savefig("clade_conservation_split.png", dpi=300)
plt.show()