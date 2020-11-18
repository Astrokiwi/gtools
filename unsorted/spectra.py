# Import libraries to do our magic
import numpy as np
import h5py

from sys import path
path.append("../visualisation/")

from joblib import Parallel, delayed
from sph_plotter import sph_plotter

import matplotlib as mpl
mpl.use('Agg')

import pylab as P



molecular_mass = 4./(1.+3.*.76)
proton_mass_cgs = 1.6726e-24
gamma_minus_one = 5./3.-1.
boltzmann_cgs = 1.38066e-16

# def parallel_prof(

def vel_absorb(*args):
    return sph_plotter.sph_vel_absorb(*args)

if __name__ == '__main__':

    #infile = "/export/1/djw/gizmos/1012/fataradtest/snapshot_070.hdf5"
    infile = "/export/1/djw/gizmos/1020/aniso/snapshot_103.hdf5"
    nprocs = 40


    f = h5py.File(infile,"r")

    xyz_p = np.array(f["/PartType0/Coordinates"]) # kpc
    xyz_p*=1.e3 # to pc

    m_p = np.array(f["/PartType0/Masses"]) # 10^10 msun
    m_p*=1.e10 # 10^10 solar masses to solar masses

    n = m_p.size

    h_p = np.array(f["/PartType0/SmoothingLength"]) # kpc
    h_p*=1.e3 # to pc

    # for doppler broadening
    u_p = np.array(f["/PartType0/InternalEnergy"]) # 1e10 erg/g
    u_p*=1.e10 # to erg/g
    #TK_p = (gamma_minus_one/boltzmann_cgs*(molecular_mass*proton_mass_cgs)*u_p)
    
    opac_p = np.zeros(n)
    opac = 65. # cm^2/g
    opac *= 0.000208908219 # to pc**2/Msun
    opac_p += opac

    vel_p = np.array(f["/PartType0/Velocities"]) # in km/s

    theta = np.pi/8.
    phi = 0.
    dv = 1000.
    #vsteps = 10**4
    vsteps = 10**3

    prof_sum = np.zeros(vsteps)

    #gx = np.linspace(-1.,1.,9)
    #gy = np.linspace(-1.,1.,9)

    #gx = [3.]
    #gy = [3.]
    thetas = np.linspace(0.,2.*np.pi,9)
    cxy = []
    for theta in thetas:
        cxy.append([np.sin(theta)*5.,np.cos(theta)*5.])
    


    P.figure() # subplots maybe?
    x = np.linspace(-dv,dv,vsteps)
    P.xlabel(r"$v$ (km/s)")

    outp = [x]

#     cxy = []
#     for cx in gx:
#         for cy in gy:
#             cxy.append([cx,cy])

    #for phi in np.linspace(0.,np.pi/2.,5):
    #for phi in np.linspace(0.,np.pi/2.,3):
    for phi in [np.pi/4.]:
        print(phi)
        prof = np.mean(Parallel(n_jobs=nprocs)(delayed(vel_absorb)(xyz_p[:,0],xyz_p[:,1],xyz_p[:,2],m_p,h_p,u_p,opac_p,vel_p[:,0],vel_p[:,1],vel_p[:,2],t_cxy[0],t_cxy[1],theta,phi,dv,vsteps,n) for t_cxy in cxy),axis=0)
        #prof = prof/np.max(prof)
        P.plot(x,prof,label=r"$\phi={}^\circ$".format(phi*180./np.pi))
        outp.append(prof)

    P.legend()
    P.savefig("pics/proftest.png")
    P.close('all')

    np.savetxt("data/proftest.dat",np.array(outp).T)


























