import h5py
import os
import torch

filename = "test2.h5"

try:
    hf = h5py.File(filename, 'r')
except OSError:
    os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
    hf = h5py.File(filename, 'r')
    
stats_key = list(hf.keys())[0]
stats = torch.tensor(hf.get(stats_key))

#you can read the corresponding text file to index specific stats

    
