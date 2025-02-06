#!/bin/bash
#SBATCH --job-name=polymer_sim
#SBATCH --output=logs/output_%j.log
#SBATCH --error=logs/error_%j.log
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G


cd $SLURM_SUBMIT_DIR


module load miniconda/3
module load gromacs/2019.3
conda activate md_env


python main.py
