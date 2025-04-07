import os
import shutil
import numpy as np
import subprocess  # Module to run shell commands from Python

template_path = "/home/IITB/CompnalAstrphyNRelvity/23n0315/MESA/mesa-24.08.1/binary_pop/template"
runs_path = "/home/IITB/CompnalAstrphyNRelvity/23n0315/MESA/mesa-24.08.1/binary_pop/runs"
slurm_script_path = os.path.join(runs_path, "run_all_bins.slurm")
n_binary = 10

for i in range(1, n_binary + 1):
    bin_name = f"bin{i:02d}"
    bin_path = os.path.join(runs_path, bin_name)

    if os.path.exists(bin_path):
        shutil.rmtree(bin_path)  # remove existing directory

    shutil.copytree(template_path, bin_path)  # now safe to copy

    m1=np.random.uniform(0.8,8)
    m2=np.random.uniform(0.8,m1)
    period=np.random.uniform(1,100)
#    print(type(m1),'/n',m2,'/n',period)

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

    # Then run the compiled binary using ./rn
    try:
       subprocess.run(["./rn"], cwd=bin_path, check=True)# Run simulation
       print(f"Simulation started in {bin_name}")
    except subprocess.CalledProcessError as e:
       print(f"Error running simulation in {bin_name}: {e}")
