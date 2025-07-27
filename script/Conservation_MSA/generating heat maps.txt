from Bio import AlignIO
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import matplotlib.patches as patches

# Define the amino acid groups based on the image reference
amino_acid_groups = {
    'Aliphatic': ['L', 'V', 'I'],
    'Glycine': ['G'],
    'Proline_Alanine': ['P', 'A'],
    'Aromatic': ['W', 'Y'],
    'Amidic': ['N', 'Q'],
    'Hydroxyl': ['S', 'T'],
    'Charged Positive': ['R', 'K'],
    'Charged Negative': ['D', 'E'],
    'Sulfur Containing': ['C']
}

# Input list of positions and corresponding reference amino acids with the label (e.g., 'N51')
positions_reference = [
    #('F54', 93),   # Y426 -> Alignment Position 34
    #('Y445', 377),  #Y445
    #('Y129', 772),   # F516 -> Alignment Position 16
    #'('F130', 813), # M517 -> Alignment Position 661
    #('Y342', 1606), # Y753 -> Alignment Position 664
    #('W788', 1862), # W788 -> Alignment Position 696
    
    ('G67', 376),   # Y426 -> Alignment Position 34
    ('D71', 380),   # Y426 -> Alignment Position 34
    ('R375', 1859),   # Y426 -> Alignment Position 34
    ('W376', 1860),   # Y426 -> Alignment Position 34


]

# Function to extract specified columns from the MSA and calculate conservation percentage
def analyze_alignment(file_name, positions_reference):
    alignment = AlignIO.read(file_name, "fasta")
    results = []
    
    for label, pos in positions_reference:
        column_data = [record.seq[pos-1] for record in alignment]  # Extract column data for position (1-indexed)
        total_sequences = len(column_data)
        reference_aa = label[0]  # Extract the actual amino acid (first character of the label)
        
        # Identify the group that the reference amino acid belongs to, or just use the amino acid itself if no group
        aa_group = [reference_aa]
        
        # Expand to the whole group if the reference belongs to one
        for group, members in amino_acid_groups.items():
            if reference_aa in members:
                aa_group = members
                break
        
        # Count how many amino acids in the column match the group (either exact amino acid or any from the group)
        matches_in_group = sum(1 for aa in column_data if aa in aa_group)
        
        # Calculate the percentage of sequences in this column that belong to the same group
        percentage_conservation = (matches_in_group / total_sequences) * 100
        results.append(percentage_conservation)
    
    return results

# List of alignment files (you can provide more files here)
alignment_files = [
    r"C:\Users\jakea\Desktop\Conservation_plotting\GluR0_tree_clades\clades\Clade A1.txt",
    r"C:\Users\jakea\Desktop\Conservation_plotting\GluR0_tree_clades\clades\Clade A2.txt",
    r"C:\Users\jakea\Desktop\Conservation_plotting\GluR0_tree_clades\clades\Clade B.txt",
    r"C:\Users\jakea\Desktop\Conservation_plotting\GluR0_tree_clades\clades\Clade C.txt",
    r"C:\Users\jakea\Desktop\Conservation_plotting\GluR0_tree_clades\clades\Clade D.txt",
    r"C:\Users\jakea\Desktop\Conservation_plotting\GluR0_tree_clades\clades\Clade E.txt",
    r"C:\Users\jakea\Desktop\Conservation_plotting\GluR0_tree_clades\clades\Clade F.txt",
    r"C:\Users\jakea\Desktop\Conservation_plotting\GluR0_tree_clades\clades\Clade G.txt",
    r"C:\Users\jakea\Desktop\Conservation_plotting\GluR0_tree_clades\clades\Clade H.txt",
    r"C:\Users\jakea\Desktop\Conservation_plotting\GluR0_tree_clades\clades\Clade I.txt",
    r"C:\Users\jakea\Desktop\Conservation_plotting\GluR0_tree_clades\clades\Clade J.txt",
    r"C:\Users\jakea\Desktop\Conservation_plotting\GluR0_tree_clades\clades\Clade K.txt",
    r"C:\Users\jakea\Desktop\Conservation_plotting\GluR0_tree_clades\clades\Clade L.txt"
    
     #r"C:\Users\jakea\Desktop\Conservation_plotting\GluR0_tree_clades\Clade A breakdown\Clade A1_1.txt",
    # r"C:\Users\jakea\Desktop\Conservation_plotting\GluR0_tree_clades\Clade A breakdown\Clade A1_2.txt",
     #r"C:\Users\jakea\Desktop\Conservation_plotting\GluR0_tree_clades\Clade A breakdown\Clade A1_3.txt",
     #r"C:\Users\jakea\Desktop\Conservation_plotting\GluR0_tree_clades\Clade A breakdown\Clade A1_4.txt",
     #r"C:\Users\jakea\Desktop\Conservation_plotting\GluR0_tree_clades\Clade A breakdown\Clade A1_5.txt",
     #r"C:\Users\jakea\Desktop\Conservation_plotting\GluR0_tree_clades\Clade A breakdown\Clade A1_6.txt",
     #r"C:\Users\jakea\Desktop\Conservation_plotting\GluR0_tree_clades\Clade A breakdown\Clade A1_7.txt",
     #r"C:\Users\jakea\Desktop\Conservation_plotting\GluR0_tree_clades\Clade A breakdown\Clade A1_8.txt",
     #r"C:\Users\jakea\Desktop\Conservation_plotting\GluR0_tree_clades\Clade A breakdown\Clade A1_9.txt",
     #r"C:\Users\jakea\Desktop\Conservation_plotting\GluR0_tree_clades\Clade A breakdown\Clade A1_10.txt",
     #r"C:\Users\jakea\Desktop\Conservation_plotting\GluR0_tree_clades\Clade A breakdown\Clade A1_11.txt",
     #r"C:\Users\jakea\Desktop\Conservation_plotting\GluR0_tree_clades\Clade A breakdown\Clade A1_12.txt",
     #r"C:\Users\jakea\Desktop\Conservation_plotting\GluR0_tree_clades\Clade A breakdown\Clade A1_13.txt"
    
    
]

# Analyze each alignment file and create a list of results
heatmap_data = []
file_labels = []
for alignment_file in alignment_files:
    heatmap_data.append(analyze_alignment(alignment_file, positions_reference))
    
    # Extract the relevant part of the filename (e.g., A1 from Clade A1.txt)
    file_label = os.path.splitext(os.path.basename(alignment_file))[0].split()[-1]
    file_labels.append(file_label)

# Convert the data into a DataFrame for visualization, using the custom labels for positions (e.g., 'N51')
df = pd.DataFrame(heatmap_data, columns=[label for label, _ in positions_reference],
                  index=file_labels)

# Plotting the heatmap without annotations
plt.figure(figsize=(10, 6))  # Adjusting figure size to control the square size
ax = sns.heatmap(df, annot=False, cmap="coolwarm", cbar_kws={'label': 'Conservation (%)'}, 
            vmin=0, vmax=100, linewidths=0.5, square=True)  # 'square=True' makes each cell square

# Customizing font sizes
plt.xticks(rotation=35, fontsize=12)  # X-axis tick labels font size
plt.yticks(rotation=360, fontsize=12)               # Y-axis tick labels font size
plt.xlabel('ssGluR0 Hinge Residues', fontsize=14, labelpad=20)  # X-axis title font size
plt.ylabel('Clade', fontsize=14, labelpad=20)      # Y-axis title font size

# Customizing the legend
cbar = ax.collections[0].colorbar  # Access the color bar (legend)
cbar.ax.tick_params(labelsize=12)  # Legend tick labels font size
cbar.set_label('Conservation (%)', fontsize=14)  # Legend title font size

# Add a black border around the heatmap
for _, spine in ax.spines.items():
    spine.set_visible(True)
    spine.set_color('black')
    spine.set_linewidth(1.5)

# Add a black border around the colorbar (legend)
cbar_outline = patches.Rectangle(
    (0, 0), 1, 1, transform=cbar.ax.transAxes,
    linewidth=1.5, edgecolor='black', facecolor='none'
)
cbar.ax.add_patch(cbar_outline)

# Show the final plot
plt.title('Anchor Conservation Heatmap', fontsize=16, pad=20)  # Main title font size
plt.tight_layout()  # Adjust layout to prevent overlap
plt.savefig("Anchor_heatmap_blackborder_4_2_CladesA_L.svg", dpi=400)
plt.show()