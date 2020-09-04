import h5py
import numpy as np
from scipy.io import FortranFile


# f = h5py.File("/srv/djw1g16/gizmos/ICs/open_binary_test_snapshot_238.hdf5","r")
f = h5py.File("/srv/djw1g16/gizmos/ICs/open_binary_inslice_cont_snapshot_147.hdf5","r")
r = np.array(f["/PartType0/Coordinates"])
v = np.array(f["/PartType0/Velocities"])
u = np.array(f["/PartType0/InternalEnergy"])
print(r.shape,v.shape,u.shape)

print(r.T[0,:])
print(v.T[0,:])
print(u[0:10])

# f = FortranFile("/srv/djw1g16/gizmos/ICs/open_binary_test_snapshot_238.dat",'w')
f = FortranFile("/srv/djw1g16/gizmos/ICs/open_binary_inslice_cont_snapshot_147.dat",'w')
f.write_record(r[:,0])
f.write_record(r[:,1])
f.write_record(r[:,2])
f.write_record(v[:,0])
f.write_record(v[:,1])
f.write_record(v[:,2])
f.write_record(u)
f.close()