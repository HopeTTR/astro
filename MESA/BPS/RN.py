# ----------------------------------------------------------------------
#  MESA population builder
#  – Draw 200 Monte-Carlo binaries (M1, M2, P_Roche)
#  – Keep the first 20 whose periods lie in 1–5 days
#  – Prepare bin01…bin20, compile, and write binary_job.sh
# ----------------------------------------------------------------------

# --------------------------------- imports -------------------------------
import os, shutil, subprocess, json
import numpy as np

# --------------------------------- paths --------------------------------
template_path = (
    "/home/IITB/CompnalAstrphyNRelvity/23n0315/MESA/mesa-24.08.1/binary_pop/template"
)
runs_path = (
    "/home/IITB/CompnalAstrphyNRelvity/23n0315/MESA/mesa-24.08.1/binary_pop/BPS/runs/RN"
)

# --------------------------------- constants ----------------------------
# 2π √(R☉³ / G M☉) expressed in **days**
RL_CONST = 0.11583331125440187

# ---------------------- helper: Roche-lobe period -----------------------
def roche_period_days(m1, m2, radius_rsun=None):
    """
    Orbital period (days) at which star-1 (mass m1) fills its Roche lobe.

    Parameters
    ----------
    m1, m2        :  primary & secondary masses  (solar units)
    radius_rsun   :  stellar radius of star-1 in R☉ (optional; if None,
                     uses ZAMS power-law R ≈ M^0.8 R☉)

    Returns
    -------
    P_days        :  orbital period in **days**
    """
    if radius_rsun is None:
        radius_rsun = m1 ** 0.80  # crude ZAMS mass-radius law

    q = m1 / m2  # Eggleton uses q = M1 / M2
    f_q = 0.49 * q ** (2 / 3) / (0.6 * q ** (2 / 3) + np.log(1 + q ** (1 / 3)))
    a_rsun = radius_rsun / f_q  # separation in R☉

    return RL_CONST * np.sqrt(a_rsun ** 3 / (m1 + m2))


# ---------------------- Monte-Carlo population --------------------------
def draw_binary_population(N, rng=None):
    if rng is None:
        rng = np.random.default_rng(seed=42)  # reproducible randomness

    # -- primary mass  (Salpeter IMF 0.8–8 M☉)
    alpha, m1_min, m1_max = 2.3, 4, 8
    u = rng.random(N)
    m1 = ((m1_max ** (1 - alpha) - m1_min ** (1 - alpha)) * u + m1_min ** (1 - alpha)) ** (
        1 / (1 - alpha)
    )

    # -- mass ratio  (flat 0.1–1.0)
    q = rng.uniform(0.25, 1.0, N)
    m2 = q * m1

    # -- Roche-lobe period
    P = roche_period_days(m1, m2)

    return list(zip(m1, m2, P))


# ---------------------- 1) generate & store 200 systems -----------------
N_total = 200
tuples_all = draw_binary_population(N_total)

with open("MC_RN.json", "w") as f:
    json.dump(
        [{"m1": float(a), "m2": float(b), "P_days": float(c)} for a, b, c in tuples_all],
        f,
        indent=2,
    )



# # ---------------------- 2) use first 20 for this batch ------------------
tuples = tuples_all[:20]  # bin01 … bin20

'''# ---------------- 2) select first 20 with 1–5-day periods -------------
P_MIN, P_MAX = 1.0, 5.0     # days
filtered = [t for t in tuples_all if P_MIN <= t[2] <= P_MAX]

if len(filtered) < 20:
    print(f"⚠  Only {len(filtered)} systems have P in [{P_MIN},{P_MAX}] d; "
          "bins will be fewer than 20.")

tuples = filtered[:20]      # these feed bin01…bin20
'''

# ---------------------- 3) loop over those 20 systems -------------------
for i, (m1, m2, period) in enumerate(tuples, start=1):
    bin_name = f"bin{i:02d}"
    bin_path = os.path.join(runs_path, bin_name)
    sh_script_path = os.path.join(bin_path, "binary_job.sh")

    # fresh copy of template
    if os.path.exists(bin_path):
        shutil.rmtree(bin_path)
    shutil.copytree(template_path, bin_path)

    # --- edit inlist_project -------------------------------------------
    inlist_path = os.path.join(bin_path, "inlist_project")
    with open(inlist_path, "r") as f:
        lines = f.readlines()
    with open(inlist_path, "w") as f:
        for line in lines:
            if "m1 =" in line:
                f.write(f"   m1 = {m1:.6f}d0  ! donor mass in Msun\n")
            elif "m2 =" in line:
                f.write(f"   m2 = {m2:.6f}d0  ! companion mass in Msun\n")
            elif "initial_period_in_days" in line:
                f.write(f"   initial_period_in_days = {period:.6f}d0\n")
            else:
                f.write(line)

    # --- compile -------------------------------------------------------
    try:
        subprocess.run(["./mk"], cwd=bin_path, check=True)
        print(f"[{bin_name}] compilation OK")
    except subprocess.CalledProcessError as e:
        print(f"[{bin_name}] compilation failed: {e}")
        continue

    # --- write SLURM script -------------------------------------------
    with open(sh_script_path, "w") as f:
        f.write(
            f"""#!/bin/bash
#SBATCH --job-name=MESA_Binary
#SBATCH --output=mesa_binary_out.log
#SBATCH --error=mesa_binary_err.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --partition=small

cd {bin_path}
./rn
"""
        )

    # optional auto-submit
    # subprocess.run(["sbatch", sh_script_path], check=True)

print("Finished preparing bin01–bin20; full population saved to montecarlo_population.json.")
