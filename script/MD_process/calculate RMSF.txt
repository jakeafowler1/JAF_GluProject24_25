import pandas as pd
import matplotlib.pyplot as plt
import glob
import os

# Function to read data from a single .agr file and plot the B factor fluctuations
def plot_RMSF(file_path, ax, color, linewidth=2):
    try:
        # Read the data, ignoring lines starting with '@'
        data = pd.read_csv(file_path, delim_whitespace=True, comment='@', header=None, names=['Res', 'RMSF'])
        
        # Plot the B factor fluctuations
        ax.plot(data['Res'], data['RMSF'], label=os.path.relpath(file_path), color=color, linewidth=linewidth)
        
        return max(data['Res'])
    except Exception as e:
        print(f"Error plotting file: {file_path}. Error message: {str(e)}")
        return 0

# Get a list of all .agr files in the current directory and subdirectories
file_paths = glob.glob('**/*fluctpro.agr', recursive=True)

# Initialize lists to store all RMSF data
all_RMSF = []

# Loop over each file and collect all RMSF data
for file_path in file_paths:
    try:
        data = pd.read_csv(file_path, delim_whitespace=True, comment='@', header=None, names=['Res', 'RMSF'])
        all_b_factors.extend(data['RMSF'].tolist())
    except Exception as e:
        print(f"Error reading file: {file_path}. Error message: {str(e)}")

# Calculate y-axis limits based on the range of all B factor data
y_min = 0
y_max = 20
x_min = 0
x_max = 1348

# Create a figure and axis for plotting with a larger size
fig, ax = plt.subplots(figsize=(25, 6))

# Loop over each file and plot its B factor data with a specified color
colors = ['darkblue', 'darkgreen', 'lightsteelblue', 'orange', 'black']
max_res = 0
for i, file_path in enumerate(file_paths):
    res = plot_RMSF(file_path, ax, color=colors[i % len(colors)], linewidth=2)
    max_res = max(max_res, res)

# Add vertical dotted lines indicating chain breaks at intervals of 337
for x in range(337, int(max_res) + 337, 337):
    ax.axvline(x=x, color='k', linestyle='--', linewidth=1.5)

# Set y-axis limits to span the range of B factor values
ax.set_ylim(y_min, 18)
ax.set_xlim(x_min, x_max)

# Add labels and title with bold and larger font sizes
ax.set_xlabel('Residue Number', fontsize=22, fontweight='bold')
ax.set_ylabel('RMSF, Å', fontsize=22, fontweight='bold')
ax.set_title('RMSF', fontsize=25, fontweight='bold')

# Adjust font size of the ticks
ax.tick_params(axis='both', which='major', labelsize=22)

# Position the legend outside the plot
ax.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=16)

# Adjust layout to make room for the legend
plt.tight_layout(rect=[0, 0, 0.75, 1])

# Save the plot as a PNG file with 500 DPI resolution
plt.savefig(' RMSF_line_plot_2.svg', dpi=500)

# Show the plot
plt.show()