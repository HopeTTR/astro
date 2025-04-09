import os
import shutil

# Paths
runs_path = "/home/IITB/CompnalAstrphyNRelvity/23n0315/MESA/mesa-24.08.1/binary_pop/runs"
output_folder = "/home/IITB/CompnalAstrphyNRelvity/23n0315/MESA/all_binary_histories"

# Create output folder if it doesn't exist
os.makedirs(output_folder, exist_ok=True)

# Loop through all binXX directories
for i in range(1, 1001):  # Adjust the range based on how many bins you have
    bin_name = f"bin{i:02d}"
    bin_path = os.path.join(runs_path, bin_name)
    source_file = os.path.join(bin_path, "binary_history.data")

    if os.path.isfile(source_file):
        dest_file = os.path.join(output_folder, f"binary_history_{bin_name}.data")
        shutil.copy(source_file, dest_file)
        print(f"Copied from {bin_name} to {dest_file}")
    else:
        print(f"No binary_history.data found in {bin_name}, skipping.")
