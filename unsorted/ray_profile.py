print("Importing")

# Import libraries to do our magic
import numpy as np
import h5py
import sys

from sys import path
path.append("../visualisation/")
from sph_plotter import sph_plotter

import gizmo_tools
path.append("../src/")
import tab_interp

sigma_convert = gizmo_tools.msun_g/(gizmo_tools.molecular_mass*gizmo_tools.proton_mass_cgs)/(3.086e+18)**2


if __name__ == '__main__':
#     run_id = "2014"
#     output_dir = "q2redo"
#     snap_str = "200"

    run_id = "3001"
    output_dir = "a2_e01"
    snap_str = "020"
#     snap_str = "100"

    rcut=80.

    gizmoDir = gizmo_tools.getGizmoDir(run_id)
    
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
#     vrad = np.sum(vel_p*xyz,axis=1)/rad3d_p

    TK_p = (gizmo_tools.gamma_minus_one/gizmo_tools.boltzmann_cgs*(gizmo_tools.molecular_mass*gizmo_tools.proton_mass_cgs)*u_p)

    # table of line fluxes
    rho_p = np.array(f["/PartType0/Density"])
    rho_p*=6.77e-22 # to g/cm**3 
    nH_p = rho_p/(gizmo_tools.molecular_mass*gizmo_tools.proton_mass_cgs)

    tau_p = np.array(f["/PartType0/AGNDepth"])
    flux_p = np.array(f["/PartType0/AGNIntensity"]) # energy per surface area per time
    flux_p*=1.989e+53/(3.086e21)**2/(3.08568e+16)

    print("Loading table")
    chTab = tab_interp.CoolHeatTab( ("coolheat_tab_marta/shrunk_table_labels_161118tau.dat"),
                                    ("coolheat_tab_marta/shrunk_table_161118_m0.0001_hsmooth_tau.dat"),
                                    ("coolheat_tab_marta/shrunk_table_labels_161118taunodust.dat"),
                                    ("coolheat_tab_marta/shrunk_table_161118_m0.0001_hsmooth_taunodust.dat"),
                                    ("coolheat_tab_marta/shrunk_table_labels_161118taudense.dat"),
                                    ("coolheat_tab_marta/shrunk_table_161118_m0.0001_hsmooth_taudense.dat")
                                    )
    interpTabVec = np.vectorize(chTab.interpTab)
    print("Interpolating from table")
    tabStructs = interpTabVec(nH_p.astype(np.float64),TK_p.astype(np.float64),flux_p.astype(np.float64),tau_p.astype(np.float64))
    line_emission = dict()
    line_names = ["co1","co2","hcn1","hcn2"]
    for lineVal in line_names:
        line_emission[lineVal] = np.array(list(map(lambda y: y.__getattr__("line_"+lineVal), tabStructs)))
    line_emission["all"] = np.zeros(m_p.shape)+1.
    
#     xyzray = np.array([[0.,1.,0.],[1.,0.,0.],[0.,0.,1.],[0.,-1.,0.],[-1.,0.,0.],[0.,0.,-1.]])
#     nray = 6

#     nray = 20
#     nray = 3

    phi_degrees = [0.,45.,89.]

    nray = len(phi_degrees)
    
#     phis = np.linspace(0.,np.pi/2.,nray)
#     phis = np.linspace(0.,np.pi/4.,nray)
    phis = np.array(phi_degrees)/360.*2.*np.pi
    thetas = np.linspace(0.,0.,nray)

    xyzray = np.zeros((nray,3))
    xyzray[:,0] = np.cos(phis)*np.cos(thetas)
    xyzray[:,1] = np.cos(phis)*np.sin(thetas)
    xyzray[:,2] = np.sin(phis)


    
    broaden = np.empty([n]) # need to update numpy to get np.full...
#     broaden[:] = -1. # no broadening
    broaden[:] = 2. # fixed broadening
    
    # broaden by sound speed
#     broaden = np.sqrt(10.*u_p/9.)/1.e5 # cm/s to km/s

#     # flatten things, for better res
#     xyz[:,0] = rad2d_p
#     xyz[:,1] = 0.
#     
#     #normalisation
#     m_p*=h_p/(2.*np.pi*rad3d_p)

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
#         for f_base,Tlowcut,Thighcut in [("vprof_all",0.,1.e10),("vprof_cold",0.,100.),("vprof_allwarm",500.,10000.),("vprof_warm",900.,1100.),("vprof_hot",5000.,15000.)]:
#         for f_base,Tlowcut,Thighcut in [("vprof_all",0.,1.e10)]:#,("vprof_cold",0.,100.),("vprof_allwarm",500.,10000.),("vprof_warm",900.,1100.),("vprof_hot",5000.,15000.)]:
        for line_name in line_names + ["all"]:
            vmin = -250
            vmax = 250.
            nbins = 100
        
#             Tlowcut = 900.
#             Thighcut = 1100.
#             mask = (TK_p>=Tlowcut) & (TK_p<=Thighcut)# & (rad2d_p<rcut)
            mask = TK_p>=-1.

            vbins = np.arange(vmin,vmax,(vmax-vmin)/nbins)            
            outtable = [vbins.T]
            
            line_flux_p = line_emission[line_name]*m_p*gizmo_tools.msun_g
            print([f(line_flux_p) for f in (np.min,np.median,np.mean,np.max)])
#             line_flux_p = m_p

            for iray in range(nray):
                print("ray:",iray+1,"/",nray)
#             if nray!=1: raise Exception("Only one ray at a time!")
                raynorm = np.sqrt(np.sum(xyzray[iray,:]**2))
                vrad = np.sum(xyzray[iray,:]*vel_p,axis=1)/raynorm # velocity along line of sight
#                 prof = sph_plotter.sph_ray_histogram(xyz,m_p,h_p,vrad,vmin,vmax,xyzray,mask,broaden,nbins,nray,n)
                thisray = np.empty([1,3])
                thisray[0,:] = xyzray[iray,:]
                prof = sph_plotter.sph_ray_histogram(xyz,line_flux_p,h_p,vrad,vmin,vmax,thisray,mask,broaden,nbins,True,1,n)
                prof*=sigma_convert
                print(np.max(prof))
                outtable+=[prof.T]
    
    
#             outtable = np.vstack((vbins.T,prof.T)).T

            outtable = np.vstack(outtable).T    
            np.savetxt("data/line{}{}.dat".format(line_name,snap_str),outtable)
        