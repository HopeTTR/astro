import os
import shutil
import subprocess

import numpy as np
import pandas as pd

# -------------------------------
# Mass‐distribution routines
# -------------------------------

# Gravitational constant in R_sun^3 / M_sun / day^2
G = 2945.49  

def stellar_radius(M):
    """Main‐sequence radius law: R ∝ M^0.8 (in R_sun)."""
    return M**0.8

def sample_power_law(alpha, x_min, x_max, size=1):
    """
    Sample from a power‐law PDF: dN/dM ∝ M^(-alpha),
    using inverse‐transform sampling.
    """
    r = np.random.uniform(0, 1, size)
    exponent = 1.0 - alpha
    return ((x_max**exponent - x_min**exponent) * r + x_min**exponent) ** (1.0 / exponent)

def generate_CE_binary():
    """
    Generate one binary likely to undergo CE:
    - M1 from Salpeter‐like (α=2.3) between 0.8–8 Msun
    - M2 = q * M1, q∈[0.1,1]
    - Compute radii, pick separation factor f, get a & P via Kepler’s law
    Returns: M1, M2, R1, R2, f, a, P
    """
    M1 = sample_power_law(alpha=2.3, x_min=0.8, x_max=8.0, size=1)[0]
    q  = np.random.uniform(0.1, 1.0)
    M2 = q * M1

    R1 = stellar_radius(M1)
    R2 = stellar_radius(M2)

    f  = np.random.uniform(5.0, 40.0)
    a  = f * (R1 + R2)  # separation in R_sun

    P  = 2 * np.pi * np.sqrt(a**3 / (G * (M1 + M2)))  # days

    return M1, M2, R1, R2, f, a, P

# -------------------------------
# MESA‐run boilerplate
# -------------------------------

# Paths (adjust to your setup)
template_path    = "/home/IITB/CompnalAstrphyNRelvity/23n0315/MESA/mesa-24.08.1/binary_pop/template"
runs_path        = "/home/IITB/CompnalAstrphyNRelvity/23n0315/MESA/mesa-24.08.1/binary_pop/runs2"

# How many binaries to set up
n_binary = 10

for i in range(1, n_binary + 1):
    bin_name    = f"bin{i:02d}"
    bin_path    = os.path.join(runs_path, bin_name)
    sh_script   = os.path.join(bin_path, "binary_job.sh")

    # (Re)create the run directory
    if os.path.exists(bin_path):
        shutil.rmtree(bin_path)
    shutil.copytree(template_path, bin_path)

    # Generate masses & period via your distribution
    M1, M2, R1, R2, f, a, period = generate_CE_binary()

    # Edit inlist_project
    inlist_path = os.path.join(bin_path, "inlist_project")
    with open(inlist_path, "r") as file:
        content = file.readlines()
    with open(inlist_path, "w") as file:
        for line in content:
            if "m1 =" in line:
                file.write(f"   m1 = {M1:.6f}d0  ! donor mass in Msun\n")
            elif "m2 =" in line:
                file.write(f"   m2 = {M2:.6f}d0  ! companion mass in Msun\n")
            elif "initial_period_in_days" in line:
                file.write(f"   initial_period_in_days = {period:.6f}d0\n")
            else:
                file.write(line)

    # Compile
    try:
        subprocess.run(["./mk"], cwd=bin_path, check=True)
        print(f"[{bin_name}] Compilation successful")
    except subprocess.CalledProcessError as e:
        print(f"[{bin_name}] Compilation failed: {e}")
        continue

    # Write SLURM script
    with open(sh_script, "w") as f:
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

echo "Starting MESA binary simulation in {bin_name}"
cd {bin_path}
./rn
echo "Finished {bin_name}"
""")
    print(f"[{bin_name}] SLURM script written")

    # Submit
    try:
        subprocess.run(["sbatch", "binary_job.sh"], cwd=bin_path, check=True)
        print(f"[{bin_name}] Job submitted")
    except subprocess.CalledProcessError as e:
        print(f"[{bin_name}] Submission failed: {e}")
