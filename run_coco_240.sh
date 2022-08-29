#!/bin/bash
#SBATCH -o coco_240_output_files/coco_output.sh.log-%j-%a
#SBATCH -n 48
#SBATCH -N 1
#SBATCH --array=1-10%2

# Load any required modules
source /etc/profile
#module load anaconda/2021b

# Run Script
#matlab -nodisplay -r "cd('$HOME/TTM_dev'); generateMultipleMongrelsFromListParallel('coco_test/coco_test.txt','default_fulliter.job'); exit;"
#matlab -nodisplay -r "cd('$HOME/TTM_dev'); generateMultipleMongrelsFromListParallel('coco_test/coco_foveated_edge.txt','default_fulliter.job'); exit;"

echo "SLURMARRAYTASK: "$SLURM_ARRAY_TASK_ID
job_txt_file=$(ls ./coco_testset_lists_240/*_${SLURM_ARRAY_TASK_ID}.txt)

matlab -nodisplay -r "cd('/home/gridsan/groups/RosenholtzLab/TTM'); 
generateMultipleMongrelsFromListParallel('$job_txt_file','default_fulliter_60olap.job'); exit;"


# #run this with:
# #LLsub run_coco_ttm.sh -N
# #OR
# #sbatch run_coco_ttm.sh
