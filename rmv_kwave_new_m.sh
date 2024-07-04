#!/bin/bash
#SBATCH --nodelist=worker8
#SBATCH --output="slurm-%j.out"
#SBATCH --time=10000000
#SBATCH --partition=thinkstation-p360
#SBATCH --gres=gpu:1

srun /usr/local/MATLAB/R2023b/bin/matlab -nosplash -nodesktop -nodisplay -r "kwave_dataset_AC; exit"