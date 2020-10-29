import numpy as np

N_res = 100000
mtot = 10.
nblob = 100
blob_n = N_res//nblob
blob_m = mtot/nblob
blobrad = 0.005
blobtemp = 1.e3

data=np.loadtxt("data/frame000098.dat")

