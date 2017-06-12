print("Importing")

# Import libraries to do our magic
import numpy as np
import sys
import os
#from joblib import Parallel, delayed

#from sys import path
#path.append("../src/")
import sph_frame

print("Running")

run_id = "1003"
output_dirs = ["L_abstest"+str(x) for x in range(7)]

snapfs = np.zeros(len(output_dirs))

for idir in range(len(output_dirs)):
    output_dir = output_dirs[idir]
    snapi = 0
    fullDir = "/export/1/djw/gizmos/"+run_id+"/"+output_dir

    fnames = os.listdir(fullDir)
    fnames = np.array(fnames)
    fnames.sort()
    snapshotfilebools = np.array([x.startswith("snapshot") for x in fnames])
    snapshotfiles = fnames[snapshotfilebools]
    snapf = 0
    ctime = os.path.getmtime(fullDir+"/snapshot_000.hdf5")
    for fname in snapshotfiles[1:]:
        new_snapf = int(fname[9:12])
        new_ctime = os.path.getmtime(fullDir+"/"+fname)
        if ( new_ctime>ctime ) :
            ctime = new_ctime
            snapf = new_snapf
    snapfs[idir] = snapf

snapf_min = np.min(snapfs)


for idir in range(len(output_dirs)):
    output_dir = output_dirs[idir]
    fullDir = "/export/1/djw/gizmos/"+run_id+"/"+output_dir
    infile = fullDir+"/snapshot_"+("%03d" % snapf_min)+".hdf5"
    outfile = "../pics/final_sphplot"+run_id+output_dir+"_%03d.png"%snapf_min
    sph_frame.makesph_trhoz_frame(infile,outfile,cmap="plasma")
