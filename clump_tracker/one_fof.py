import numpy as np
import h5py
from fof import fof

infile = "/export/1/djw/gizmos/1032/goodsurf_fat_fixedtable/snapshot_270.hdf5"

f = h5py.File(infile,"r")

xyz = np.array(f["/PartType0/Coordinates"])
xyz*= 1.e3 # to pc

xyz = xyz[:5000,:]

#xyz[:,2] = 0.

grp = fof.calc_fof(xyz,.1)

print(grp)

print(np.unique(grp))

uniq_grps = np.unique(grp)

np.savetxt("all",xyz)

for igrp in uniq_grps:
    xyz_out = xyz[(grp==igrp),:]
#     if ( xyz_out.size//3>2 ):
    print("output{} size={}".format(igrp,xyz_out.size//3))
    np.savetxt("group{}".format(igrp),xyz_out)
