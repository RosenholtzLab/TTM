#!/bin/bash

# Slurm sbatch options
#SBATCH -o output.sh.log-%j
#SBATCH -n 8

# Load any required modules
source /etc/profile
#module load anaconda/2021b

# Run Script
#matlab -nodisplay -r "cd('$HOME/TTM_dev'); generateMultipleMongrelsFromListParallel('coco_test/coco_test.txt','default_fulliter.job'); exit;"
/state/partition1/llgrid/pkg/matlabr2021b/bin/matlab -nodisplay -r "cd('/home/gridsan/groups/RosenholtzLab/TTM'); generateMultipleMongrelsFromListParallel('coco_test/debug.txt','default_rediter_60olap.job'); exit;"


##run this with:
##LLsub run_coco_ttm.sh -N
##OR
##sbatch run_coco_ttm.sh
