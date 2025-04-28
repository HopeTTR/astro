import os
import shutil
import subprocess

import numpy as np
import pandas as pd

G = 2945.49  # Gravitational constant in units consistent with R_sun, M_sun, days

def stellar_radius(M):
    """Main-sequence radius law: R \u221d M^0.8 (in R_sun)."""
    return M**0.8

def roche_lobe_radius_fraction(q):
    """
    Eggleton (1983) approximation for the primary's Roche-lobe fraction:
      RL1/a = 0.49 q^(2/3) / [0.6 q^(2/3) + ln(1 + q^(1/3))]
    where q = M1/M2.
    """
    return 0.49 * q**(2/3) / (0.6 * q**(2/3) + np.log(1 + q**(1/3)))

def compute_CE_binary_params(M1, M2):
    """
    Given M1, M2 (Msun), compute:
      R1, R2      -- stellar radii (R_sun)
      a           -- separation (R_sun) so primary exactly fills Roche lobe
      P           -- orbital period (days) via Kepler's third law
    Returns: R1, R2, a, P
    """
    R1 = stellar_radius(M1)
    R2 = stellar_radius(M2)

    q = M1 / M2
    RL1_frac = roche_lobe_radius_fraction(q)

    a = R1 / RL1_frac
    P = 2 * np.pi * np.sqrt(a**3 / (G * (M1 + M2)))

    return R1, R2, a, P

def generate_CE_binary():
    """
    Generate one CE-onset binary:
    - M1 from 2-4 Msun (uniform)
    - q from 0.1 to 1.0
    - M2 = q * M1
    - Compute R1, R2, a, P
    Returns: M1, M2, R1, R2, a, P
    """
    M1 = np.random.uniform(2, 4)
    q = np.random.uniform(0.1, 1.0)
    M2 = q * M1

    R1, R2, a, P = compute_CE_binary_params(M1, M2)

    return M1, M2, R1, R2, a, P
# # Paths (adjust to your setup)
template_path    = "/home/IITB/CompnalAstrphyNRelvity/23n0315/MESA/mesa-24.08.1/binary_pop/template"
runs_path        = "/home/IITB/CompnalAstrphyNRelvity/23n0315/MESA/mesa-24.08.1/binary_pop/BPS/runs/2to4"

# How many binaries to set up
n_binary = 10

# List to hold all generated binaries
tuples_all = []

for i in range(1, n_binary + 1):
    M1, M2, R1, R2, a, P = generate_CE_binary()
    tuples_all.append((M1, M2, P, R1, R2, a))

# Save to JSON file
output_filename = "parameters_2to4.json"

with open(output_filename, "w") as f:
    json.dump(
        [
            {"M1": float(m1), "M2": float(m2), "P": float(p), "R1": float(r1), "R2": float(r2), "a": float(a_)}
            for m1, m2, p, r1, r2, a_ in tuples_all
        ],
        f,
        indent=2,
    )

print(f"Saved {n_binary} binaries to {output_filename}")
# -------------------------------
# MESA‐run boilerplate
# -------------------------------

# How many binaries to set up
# n_binary = 10

# for i in range(1, n_binary + 1):
    bin_name    = f"bin{i:02d}"
    bin_path    = os.path.join(runs_path, bin_name)
    sh_script   = os.path.join(bin_path, "binary_job.sh")

    # (Re)create the run directory
    if os.path.exists(bin_path):
        shutil.rmtree(bin_path)
    shutil.copytree(template_path, bin_path)

    # Generate masses & period via your distribution
    M1, M2, R1, R2, a, P = generate_CE_binary()

  
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
                file.write(f"   initial_period_in_days = {P:.6f}d0\n")
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
        print(f"[{bin_name}] Submission failed: {e}")'''
