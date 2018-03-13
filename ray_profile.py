print("Importing")

# Import libraries to do our magic
import numpy as np
import h5py
import sys

from sys import path
path.append("visualisation/")
from sph_plotter import sph_plotter

import gizmo_tools

sigma_convert = gizmo_tools.msun_g/(gizmo_tools.molecular_mass*gizmo_tools.proton_mass_cgs)/(3.086e+18)**2


if __name__ == '__main__':
    run_id = "2014"
    output_dir = "q2redo"
    snap_str = "200"

    rcut=80.

    gizmoDir = gizmo_tools.getGizmoDir()
    
    f = h5py.File(gizmoDir+"/"+run_id+"/"+output_dir+"/snapshot_"+snap_str+".hdf5","r")

    xyz = np.array(f["/PartType0/Coordinates"])
    xyz*=1.e3 # to pc
    rad2d_p = np.sqrt(xyz[:,0]**2+xyz[:,1]**2)
    rad3d_p = np.sqrt(xyz[:,0]**2+xyz[:,1]**2++xyz[:,2]**2)

    m_p = np.array(f["/PartType0/Masses"])
    m_p*=1.e10 # 10^10 solar masses to Msun

    h_p = np.array(f["/PartType0/SmoothingLength"]) # kpc
    h_p*=1.e3 # to pc

    n = m_p.size

    u_p = np.array(f["/PartType0/InternalEnergy"]) # 1e10 erg/g
    u_p*=1.e10 # to erg/g

    vel_p = np.array(f["/PartType0/Velocities"])
    vrad = np.sum(vel_p*xyz,axis=1)/rad3d_p

    TK_p = (gizmo_tools.gamma_minus_one/gizmo_tools.boltzmann_cgs*(gizmo_tools.molecular_mass*gizmo_tools.proton_mass_cgs)*u_p)


    
#     xyzray = np.array([[0.,1.,0.],[1.,0.,0.],[0.,0.,1.],[0.,-1.,0.],[-1.,0.,0.],[0.,0.,-1.]])
#     nray = 6

    nray = 20

    phis = np.linspace(0.,np.pi/4.,nray)
    thetas = np.linspace(0.,0.,nray)

    xyzray = np.zeros((nray,3))
    xyzray[:,0] = np.cos(phis)*np.cos(thetas)
    xyzray[:,1] = np.cos(phis)*np.sin(thetas)
    xyzray[:,2] = np.sin(phis)

    # flatten things, for better res
    xyz[:,0] = rad2d_p
    xyz[:,1] = 0.
    
    #normalisation
    m_p*=h_p/(2.*np.pi*rad3d_p)

    if False:
        TK_p = np.log10(TK_p)

        Tmin = 1.
        Tmax = 4.
        nbins = 100
        mask = np.full(n,True,dtype=np.bool)

        prof = sph_plotter.sph_ray_histogram(xyz,m_p,h_p,TK_p,Tmin,Tmax,xyzray,mask,nbins,nray,n)
        prof*=sigma_convert
    
        Tbins = np.arange(Tmin,Tmax,(Tmax-Tmin)/nbins)
    
        outtable = np.vstack((Tbins.T,prof.T)).T
    
        np.savetxt("data/proftest.dat",outtable)
    
    if True:
        for f_base,Tlowcut,Thighcut in [("vprof_all",0.,1.e10),("vprof_cold",0.,100.),("vprof_allwarm",500.,10000.),("vprof_warm",900.,1100.),("vprof_hot",5000.,15000.)]:
    
            vmin = -50
            vmax = 250.
            nbins = 100
        
#             Tlowcut = 900.
#             Thighcut = 1100.
            mask = (TK_p>=Tlowcut) & (TK_p<=Thighcut) & (rad2d_p<rcut)

            prof = sph_plotter.sph_ray_histogram(xyz,m_p,h_p,vrad,vmin,vmax,xyzray,mask,nbins,nray,n)
            prof*=sigma_convert
    
            vbins = np.arange(vmin,vmax,(vmax-vmin)/nbins)
    
            outtable = np.vstack((vbins.T,prof.T)).T
    
            np.savetxt("data/{}.dat".format(f_base),outtable)
        