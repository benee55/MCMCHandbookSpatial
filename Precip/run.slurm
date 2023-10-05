#!/bin/sh

#SBATCH --job-name=DR_BS
#SBATCH --partition=normal
#SBATCH --output=/scratch/slee287/logFiles/%j.out  # Output file
#SBATCH --error=/scratch/slee287/logFiles/%j.err   # Error file

#SBATCH --mem-per-cpu=10G        # Total memory needed per task (units: K,M,G,T)
#SBATCH --time=02-00:00  # Total time needed for job: Days-Hours:Minutes
#SBATCH --ntasks 1          # Number of processes you plan to launch

module load r
cd /home/slee287/MCMCHandbook/Precip/run_bisquare
Rscript B_MCMC.R


