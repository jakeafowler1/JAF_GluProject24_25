import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# Define a function to convert RGB values to hex
def rgb_to_hex(rgb):
    return "#{:02x}{:02x}{:02x}".format(int(rgb[0]*255), int(rgb[1]*255), int(rgb[2]*255))

# Read data from file skipping lines starting with '@'
data = []
with open('r375k_pro9_trajout/fluctpro.agr', 'r') as file:
    for line in file:
        if not line.startswith('@'):
            data.append([float(x) for x in line.split()])

# Convert data to numpy array
data_array = np.array(data)

# Extract residue numbers and RMSF values
residue_numbers = data_array[:, 0].astype(int)
rmsf_values = data_array[:, 1]

# Define the protein sequence (replace this with your actual sequence)
protein_sequence = "TLKVGVVGNPPFVFYGEGKNAAFTGISLDVWRAVAESQKWNSEYVRQNSISAGITAVAEGELDILIGPISVTPERAAIEGITFTQPYFSSGIGLLIPGKPVSLWERFSPFFGIAALSSAGVLTLLLFLVGNLIWLAEHRKNPEQFSPHYPEGVQNGMWFALVTLTTVGYGDRSPRTKLGQLVAGVWMLVALLSFSSITAGLASAFSTALSEASATPLFRSVGDLKNKEVAVVRDTTAVDWANFYQADVRETNNLTAAITLLQKKQVEAVMFDRPALIYYTRQNPNLNLEVTEIRVSLEPYGFVLKENSPLQKTINVEMLNLLYSRVIAEFTERWLGPTLKVGVVGNPPFVFYGEGKNAAFTGISLDVWRAVAESQKWNSEYVRQNSISAGITAVAEGELDILIGPISVTPERAAIEGITFTQPYFSSGIGLLIPGKPVSLWERFSPFFGIAALSSAGVLTLLLFLVGNLIWLAEHRKNPEQFSPHYPEGVQNGMWFALVTLTTVGYGDRSPRTKLGQLVAGVWMLVALLSFSSITAGLASAFSTALSEASATPLFRSVGDLKNKEVAVVRDTTAVDWANFYQADVRETNNLTAAITLLQKKQVEAVMFDRPALIYYTRQNPNLNLEVTEIRVSLEPYGFVLKENSPLQKTINVEMLNLLYSRVIAEFTERWLGPTLKVGVVGNPPFVFYGEGKNAAFTGISLDVWRAVAESQKWNSEYVRQNSISAGITAVAEGELDILIGPISVTPERAAIEGITFTQPYFSSGIGLLIPGKPVSLWERFSPFFGIAALSSAGVLTLLLFLVGNLIWLAEHRKNPEQFSPHYPEGVQNGMWFALVTLTTVGYGDRSPRTKLGQLVAGVWMLVALLSFSSITAGLASAFSTALSEASATPLFRSVGDLKNKEVAVVRDTTAVDWANFYQADVRETNNLTAAITLLQKKQVEAVMFDRPALIYYTRQNPNLNLEVTEIRVSLEPYGFVLKENSPLQKTINVEMLNLLYSRVIAEFTERWLGPTLKVGVVGNPPFVFYGEGKNAAFTGISLDVWRAVAESQKWNSEYVRQNSISAGITAVAEGELDILIGPISVTPERAAIEGITFTQPYFSSGIGLLIPGKPVSLWERFSPFFGIAALSSAGVLTLLLFLVGNLIWLAEHRKNPEQFSPHYPEGVQNGMWFALVTLTTVGYGDRSPRTKLGQLVAGVWMLVALLSFSSITAGLASAFSTALSEASATPLFRSVGDLKNKEVAVVRDTTAVDWANFYQADVRETNNLTAAITLLQKKQVEAVMFDRPALIYYTRQNPNLNLEVTEIRVSLEPYGFVLKENSPLQKTINVEMLNLLYSRVIAEFTERWLGP"  # Example sequence

# Check if the sequence length matches the number of residues
if len(protein_sequence) != len(rmsf_values):
    raise ValueError("The length of the protein sequence and the number of RMSF values must match.")

# Input: Specify the range of residues to color (1-based index)
start_residue = 1  # Example start residue (inclusive)
end_residue = 337    # Example end residue (inclusive)

# Define the target RMSF range
min_rmsf = 2
max_rmsf = 15

# Normalize RMSF values to [1, 15] range
normalized_rmsf = np.clip((rmsf_values - rmsf_values.min()) / (rmsf_values.max() - rmsf_values.min()) * (max_rmsf - min_rmsf) + min_rmsf, min_rmsf, max_rmsf)

# Get colors from the 'BuPu' colormap
cmap = plt.get_cmap('BuPu')
colors = cmap((normalized_rmsf - min_rmsf) / (max_rmsf - min_rmsf))

# Create an HTML file to display the sequence with colored highlights
with open('R375K_rmsf_colored_sequence.html', 'w') as html_file:
    html_file.write("<html><body>\n")
    html_file.write("<h1>Protein Sequence Colored by RMSF (Residues {} to {})</h1>\n".format(start_residue, end_residue))
    html_file.write("<pre style='font-size: 18px;'>\n")

    for i, residue in enumerate(protein_sequence):
        residue_index = i + 1  # Convert to 1-based index
        
        # If the residue is in the specified range, color it based on RMSF
        if start_residue <= residue_index <= end_residue:
            color = colors[i]
            color_hex = rgb_to_hex(color[:3])  # Convert RGB to hex for HTML
            html_file.write(f"<span style='background-color:{color_hex};'>{residue}</span>")
        else:
            # Default color for residues outside the specified range
            html_file.write(f"<span style='background-color:#FFFFFF;'>{residue}</span>")  # White background
        
        # Optionally, wrap lines every 50 residues for better readability
        if (i + 1) % 50 == 0:
            html_file.write("\n")
    
    html_file.write("</pre>\n")
    html_file.write("</body></html>\n")

print("HTML file 'RW_LF_rmsf_colored_sequence.html' has been created. Open it in a browser to view the highlighted sequence.")
