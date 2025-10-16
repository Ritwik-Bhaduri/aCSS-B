#!/bin/bash
#SBATCH -J aCSS  # A single job name for the array
#SBATCH -p shared,sapphire # Partitions
#SBATCH -c 1 # number of cores
#SBATCH -t 1-00:00  # Running time in the format - D-HH:MM
#SBATCH --mem 2000 # Memory request - 1000 corresponds to 1GB
#SBATCH -o ./log_files/output_file_%A_%a.out # Standard output
#SBATCH -e ./log_files/error_file_%A_%a.err # Standard error

module load R
export R_LIBS_USER="$HOME/apps/R/R_4.2.2:${R_LIBS_USER:-}"
Rscript Run.R #SLURM_ARRAY_TASK_ID
