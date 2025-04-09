#!/bin/bash
#SBATCH --job-name=binary_array
#SBATCH --output=slurm_out_%A_%a.out
#SBATCH --error=slurm_err_%A_%a.err
#SBATCH --array=1-8        # Adjust range to match number of binaries
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --partition=small

# Load MESA environment
#source /path/to/mesasdk_init.sh  # Replace with your real environment setup

# Compute the bin folder name (e.g., bin01, bin02, ...)
BIN_ID=$(printf "bin%02d" ${SLURM_ARRAY_TASK_ID})
BIN_DIR="/home/IITB/CompnalAstrphyNRelvity/23n0315/MESA/mesa-24.08.1/binary_pop/runs/${BIN_ID}"

echo "Running ${BIN_ID} at $(date)"
cd "$BIN_DIR"

./rn  # Run the compiled binary simulation
echo "Finished ${BIN_ID} at $(date)"
