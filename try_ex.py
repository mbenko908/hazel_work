import numpy as np
import matplotlib.pyplot as pl
import hazel

iterator = hazel.Iterator(use_mpi=True)
mod = hazel.Model('/mnt/RAIDofsF/mbenko/results/160823/inversion_he01/configuration/conf_single_he_mpi.ini', rank=iterator.get_rank(), working_mode='inversion', verbose=3)#, rank=rank, verbose=3)
iterator.use_model(model=mod)
iterator.run_all_pixels()
