print("Importing")

# Import libraries to do our magic
import numpy as np
import h5py
import sys

print("Running")

def vec_kernel(r,h):
    x = r/h
    k = np.zeros(x.size)
    # k[x<.5] = 1.-2.*x[x<.5]**2
#     k[(x<1.)&(x>=.5)] = 2.*(x[(x<1.)&(x>=.5)]-1.)**2
    #k = 1.-x
    #k = 1./(x**2.)
    k[x<.5] = 1.+6.*(x[x<.5]-1.)*x[x<.5]**2
    k[(x<1.)&(x>=.5)] = 2.*(1.-x[(x<1.)&(x>=.5)])**3
    # normalise
    k*=8./np.pi/h**3
    return k


f = h5py.File("/export/1/djw/gizmo_public/disc_nocool_midres_out/snapshot_050.hdf5","r")

r_p = np.array(f["/PartType0/Coordinates"])

h_p = np.array(f["/PartType0/SmoothingLength"])

#tograd_p = np.array(f["/PartType0/Density"])
tograd_p = np.array(f["/PartType0/InternalEnergy"])

#rad_p = np.linalg.norm(r_p,axis=1)
#tograd_p = 1.-rad_p

ip = 50000
r0 = r_p[ip]
tograd0 = tograd_p[ip]
h0 = h_p[ip]

dr3_p = r_p-r0
dr_p = np.linalg.norm(dr3_p,axis=1)
w_p = vec_kernel(dr_p,h_p)

# only loop over closest particles

close_p = (dr_p<h_p)
print(np.sum(close_p),np.sum(w_p>0.))

# find gradient at r0

#make matrix to invert
E_mat = np.zeros((3,3))

for ix in range(3):
    for iy in range(3):
        E_mat[ix,iy] = np.sum((r_p[close_p,ix]-r0[ix])*(r_p[close_p,iy]-r0[iy])*w_p[close_p])

#print(E_mat)

B_mat = np.linalg.inv(E_mat)

gradvec = np.zeros(3)

for ix in range(3):
    for iy in range(3):
        gradvec[ix]+=np.sum(B_mat[ix,iy]*(tograd_p[close_p]-tograd0)*(r_p[close_p,iy]-r0[iy])*w_p[close_p])

gradvec/=1.e3 # convert to /pc
print(gradvec)
print(r0)

print(np.sum(gradvec*r0)/np.linalg.norm(r0))
