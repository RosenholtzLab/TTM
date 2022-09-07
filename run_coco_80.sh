#!/bin/bash
#SBATCH -o coco_80_output_files/coco_output.sh.log-%j
#SBATCH -n 24
#SBATCH --array=1-89
#SBATCH --mem=75GB
# Load any required modules
source /etc/profile
#module load anaconda/2021b

# Run Script
#matlab -nodisplay -r "cd('$HOME/TTM_dev'); generateMultipleMongrelsFromListParallel('coco_test/coco_test.txt','default_fulliter.job'); exit;"
#matlab -nodisplay -r "cd('$HOME/TTM_dev'); generateMultipleMongrelsFromListParallel('coco_test/coco_foveated_edge.txt','default_fulliter.job'); exit;"

echo "SLURMARRAYTASK: "$SLURM_ARRAY_TASK_ID
job_txt_file=$(ls ../coco_testset_lists_80/*_${SLURM_ARRAY_TASK_ID}.txt)
echo $job_txt_file
matlab -nodisplay -r "cd('/home/gridsan/groups/RosenholtzLab/TTM'); generateMultipleMongrelsFromListParallel('$job_txt_file','default_fulliter_60olap.job'); exit;"


# #run this with:
# #LLsub run_coco_ttm.sh -N
# #OR
# #sbatch run_coco_ttm.sh
