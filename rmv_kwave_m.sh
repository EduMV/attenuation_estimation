#!/usr/bin/bash
#SBATCH --gpus-per-node=1
#SBATCH --nodes=1
#SBATCH --partition=thinkstation-p360
#SBATCH --nodelist=worker10
#SBATCH --output="log.out"

srun matlab -nosplash -nodesktop -nodisplay -r "kwave_dataset_AC; exit"