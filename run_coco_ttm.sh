#!/bin/bash

# Slurm sbatch options
#SBATCH -o output.sh.log-%j
#SBATCH --exclusive
#SBATCH -n 48
#SBATCH -N 1

# Load any required modules
source /etc/profile
#module load anaconda/2021b

# Run Script
#matlab -nodisplay -r "cd('$HOME/TTM_dev'); generateMultipleMongrelsFromListParallel('coco_test/coco_test.txt','default_fulliter.job'); exit;"
matlab -nodisplay -r "cd('$HOME/TTM_dev'); generateMultipleMongrelsFromListParallel('coco_test/coco_foveated_edge.txt','default_fulliter.job'); exit;"


##run this with:
##LLsub run_coco_ttm.sh -N
##OR
##sbatch run_coco_ttm.sh
