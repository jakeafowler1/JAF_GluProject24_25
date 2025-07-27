import pandas as pd
import matplotlib.pyplot as plt
import glob
import os
import re
import numpy as np

# Function to read the total simulation time from the second line of the pro3.mdinfo file
def read_simulation_time(mdinfo_path):
    try:
        with open(mdinfo_path, 'r') as file:
            file.readline()  # Skip the first line
            line = file.readline()  # Read the second line

            # Use a regular expression to extract the numeric value for TIME(PS)
            time_match = re.search(r"TIME\(PS\)\s*=\s*([\d\.]+)", line)
            if time_match:
                total_time = float(time_match.group(1))  # Convert the matched string to a float
                return total_time
            else:
                raise ValueError(f"TIME(PS) not found in {mdinfo_path}")
    except Exception as e:
        print(f"Error reading {mdinfo_path}: {e}")
        return None

# Function to read and plot temperature data on separate graphs
def plot_temperature_data(temp_path, total_time, output_dir, max_time, label, color='black', linewidth=2, ax=None, save_plot=False):
    try:
        # Read the temperature data
        data = pd.read_csv(temp_path, delim_whitespace=True, comment='#', names=['Time', 'MDOUT_TEMP'])

        if data.empty:
            print(f"No data found in {temp_path}")
            return

        # Calculate the time for each reading by dividing the total time by the number of readings
        num_readings = len(data)
        time_per_reading = total_time / num_readings
        data['Time'] = [(i + 1) * time_per_reading for i in range(num_readings)]

        # Convert time to nanoseconds by dividing by 1000
        data['Time'] = data['Time'] / 1000

        # If an axis is provided (for combined plot), use it, otherwise create a new plot
        if ax is None:
            # Create a figure and axis for plotting
            fig, ax = plt.subplots(figsize=(10, 5))

        # Plot the raw temperature data
        ax.plot(data['Time'], data['MDOUT_TEMP'], label=label, color=color, linewidth=linewidth)

        # Set x-axis limits to start from 0 and end at the global minimum max_time
        x_min = 0
        x_max = 944
        ax.set_xlim(x_min, max_time)

        # Set y-axis limits based on the data
        y_min = -655000 
        y_max = -635000
        ax.set_ylim(y_min, y_max)

        # Set x and y ticks to be more spaced out
        ax.set_xticks(np.arange(0, int(max_time) + 1, 100))  # Set ticks every 100 ns on the x-axis
        ax.set_yticks(np.linspace(y_min, y_max, 5))  # Set 5 evenly spaced ticks on the y-axis

        # Add labels and title with bold and larger font sizes and custom padding
        ax.set_xlabel('Time (ns)', fontsize=18, labelpad=20)  # Adjusted label padding
        ax.set_ylabel('Energy (kJ mol ^-1)', fontsize=18, labelpad=20)
        ax.set_title(f'Total Energy vs Time\n{label}', fontsize=20, fontweight='bold', pad=20)

        # Adjust font size of the ticks
        ax.tick_params(axis='both', which='major', labelsize=14)

        # Save the individual plot as a PNG file if required
        if save_plot:
            plt.tight_layout()
            output_file = os.path.join(output_dir, f'{label}_energy_vs_time.png')
            plt.savefig(output_file, dpi=300)
            plt.close()

    except Exception as e:
        print(f"Error plotting file: {temp_path}. Error message: {str(e)}")

# Get a list of all .dat files in the current directory and subdirectories
temp_file_paths = glob.glob('**/*etot.dat', recursive=True)

# Initialize a list to track the maximum time values from all files
max_times = []

# First loop: Find the smallest max time across all files
for temp_path in sorted(temp_file_paths):
    directory = os.path.dirname(temp_path)
    mdinfo_path = os.path.join(directory, 'pro9.mdinfo')
    total_time = read_simulation_time(mdinfo_path)

    if total_time is not None:
        data = pd.read_csv(temp_path, delim_whitespace=True, comment='#', names=['Time', 'MDOUT_TEMP'])
        if not data.empty:
            num_readings = len(data)
            time_per_reading = total_time / num_readings
            max_time = (num_readings * time_per_reading) / 1000  # Convert to nanoseconds
            max_times.append(max_time)

# Find the global minimum max_time
if max_times:
    global_max_time = min(max_times)
else:
    global_max_time = None
    print("No valid temperature data found.")

# Ensure the output directory exists
output_dir = 'temperature_plots'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Initialize a figure for the combined plot
fig_combined, ax_combined = plt.subplots(figsize=(14, 6))

# Function to extract label from directory name and capitalize
def get_label_from_directory(directory):
    base_dir = os.path.basename(directory)
    return base_dir.upper()  # Convert directory name to uppercase for legend

# Second loop: Plot each file with the global minimum max_time, individually and on a combined plot
colors = ['darkblue', 'darkgreen', 'lightsteelblue', 'orange', 'black', 'darkred', 'k']  # Define a set of colors for the lines
for i, temp_path in enumerate(sorted(temp_file_paths)):
    directory = os.path.dirname(temp_path)
    output_dir_name = os.path.basename(directory)
    mdinfo_path = os.path.join(directory, 'pro9.mdinfo')
    total_time = read_simulation_time(mdinfo_path)

    # Extract the label for the legend
    label = get_label_from_directory(directory)

    if total_time is not None and global_max_time is not None:
        # Plot individual plots and save them
        plot_temperature_data(temp_path, total_time, output_dir, global_max_time, label=label, color=colors[i % len(colors)], linewidth=1.4, save_plot=True)

        # Plot on the combined graph
        plot_temperature_data(temp_path, total_time, output_dir, global_max_time, label=label, color=colors[i % len(colors)], linewidth=1.4, ax=ax_combined)

# Customize the combined plot
ax_combined.set_xlabel('Time (ns)', fontsize=20, labelpad=20)
ax_combined.set_ylabel('Energy (kJ mol ^-1)', fontsize=20, labelpad=20)
ax_combined.set_title('Total Energy vs Time for All Simulations', fontsize=22, fontweight='bold', pad=20)

# Adjust x and y ticks for the combined plot
ax_combined.set_xticks(np.arange(0, int(global_max_time) + 1, 100))  # Set x ticks every 100 ns
y_min_combined, y_max_combined = ax_combined.get_ylim()
ax_combined.set_yticks(np.linspace(y_min_combined, y_max_combined, 5))  # Set 5 evenly spaced ticks on the y-axis

# Adjust font size of the ticks
ax_combined.tick_params(axis='both', which='major', labelsize=18)

# Add a legend to the combined plot with custom font size
ax_combined.legend(loc='best', bbox_to_anchor=(1, 1), fontsize=12)

# Adjust layout and save the combined plot
plt.tight_layout(rect=[0, 0, 0.75, 1])
combined_output_file = os.path.join(output_dir, 'combined_energy_vs_time_2.svg')
plt.savefig(combined_output_file, dpi=400)
plt.show()
