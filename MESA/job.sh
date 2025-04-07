#!/bin/bash
#SBATCH --job-name=setup_mesa_bins
#SBATCH --output=setup_mesa_binsout.log
#SBATCH --error=setup_mesa_bins_err.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --partition=small

# Load environment (optional)
# module load python
# conda activate your_env  # if using conda

echo "Starting Python MESA automation script..."
python /home/IITB/CompnalAstrphyNRelvity/23n0315/MESA/mesa-24.08.1/binary_pop/Binary_pop.py
echo "Finished submitting all binary jobs."
