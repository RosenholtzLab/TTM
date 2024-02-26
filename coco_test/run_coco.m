addpath(genpath('/data/scratch/vasha/TTM'))

tic

generateMultipleMongrelsFromListParallel('/data/scratch/vasha/TTM/coco_test/coco_list.txt', 'default_rediter.job');

toc

t = toc-tic