import os  # Module to work with file system paths and directories
import shutil  # Module to copy or delete entire folders and files
import numpy as np  # NumPy library used for generating random numbers
import subprocess  # Module to run shell commands from Python

# Path to the MESA project template directory
#template_path = "/Users/hope/Downloads/MSc_project/mesa_projects/binary_pop/template"

# Path where new binary run directories will be created
#runs_path = "/Users/hope/Downloads/MSc_project/mesa_projects/binary_pop/runs"

template_path = "/home/IITB/CompnalAstrphyNRelvity/23n0315/MESA/mesa-24.08.1/binary_pop/template"
runs_path = "/home/IITB/CompnalAstrphyNRelvity/23n0315/MESA/mesa-24.08.1/binary_pop/runs"
#slurm_script_path = os.path.join(runs_path, "run_all_bins.slurm")


# Number of binary systems you want to simulate
n_binary = 10

# Loop from 1 to n_binary (inclusive) to create bin01, bin02, ..., bin10
for i in range(1, n_binary + 1):
    bin_name = f"bin{i:02d}"  # Format binary folder name with 2 digits, e.g., 'bin01'
    bin_path = os.path.join(runs_path, bin_name)  # Construct the full path for this binary run
    sh_script_path = os.path.join(bin_path, "binary_job.sh")
    
    if os.path.exists(bin_path):  # Check if the directory already exists
        shutil.rmtree(bin_path)  # If it exists, delete it to avoid conflicts

    shutil.copytree(template_path, bin_path)  # Copy the entire template folder into binXX

    m1 = np.random.uniform(0.8, 8)  # Randomly choose m1 between 0.8 and 8 (primary star mass)
    m2 = np.random.uniform(0.8, m1)  # Randomly choose m2 between 0.8 and m1 (secondary star mass)
    period = np.random.uniform(1, 100)  # Randomly choose orbital period between 1 and 100 days

    print(m1, m2, period)  # Print the generated values for checking

    # Path to the inlist_project file inside this binary run folder
    inlist_path = os.path.join(bin_path, "inlist_project")

    # Open the inlist_project file in read mode and store its lines in a list
    with open(inlist_path, "r") as file:
        content = file.readlines()

    # Open the same file in write mode to overwrite it with updated content
    with open(inlist_path, "w") as file:
        for line in content:  # Loop through each line of the original file
            if "m1 =" in line:  # If the line contains m1 assignment
                file.write(f"   m1 = {m1:.6f}d0  ! donor mass in Msun\n")  # Replace it with new m1
            elif "m2 =" in line:  # If the line contains m2 assignment
                file.write(f"   m2 = {m2:.6f}d0  ! companion mass in Msun\n")  # Replace with new m2
            elif "initial_period_in_days" in line:  # If the line has period assignment
                file.write(f"   initial_period_in_days = {period:.6f}d0\n")  # Replace with new period
            else:
                file.write(line)  # Otherwise, write the line unchanged

    # First compile the MESA binary project using ./mk
    try:
        subprocess.run(["./mk"], cwd=bin_path, check=True)  # Compile in the binXX directory
        print(f"Compilation successful in {bin_name}")
    except subprocess.CalledProcessError as e:
        print(f"Compilation failed in {bin_name}: {e}")
        continue  # Skip running ./rn if compilation fails
'''

# -----------------------------------------#

    with open(sh_script_path, "w") as f:
        f.write(f"""#!/bin/bash
#SBATCH --job-name=MESA_Binary
#SBATCH --output=mesa_binary_out.log
#SBATCH --error=mesa_binary_err.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --partition=small

# Load your Conda environment (uncomment if needed)
# module load conda
# conda activate mesa_env  # Replace with your actual MESA conda environment

echo "Submitting MESA binary simulation..."

cd {bin_path}  # Navigate to the specific binary run directory

./rn  # Run the binary simulation

echo "Simulation completed."
""")

        print(f"Created binary_job.sh in {bin_name}")

#-------------------------------------#
    job_script = os.path.join(bin_path, "binary_job.sh")  # Full path to the SLURM job script inside that bin folder

    # Check if binary_job.sh exists in the current bin folder
    if not os.path.isfile(job_script):
        print(f" Job script missing in {bin_name}, skipping...")  # Warn if the script doesn't exist
        continue  # Skip this bin and move to the next
        '''
