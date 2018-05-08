print("Importing")

# Import libraries to do our magic
import numpy as np
import h5py
import sys
# import matplotlib as mpl
# mpl.use('Agg')
# 
# import pylab as P

from sys import path
path.append("visualisation/")
from sph_plotter import sph_plotter

import itertools as it
import gizmo_tools


molecular_mass = 4./(1.+3.*.76)
proton_mass_cgs = 1.6726e-24
gamma_minus_one = 5./3.-1.
boltzmann_cgs = 1.38066e-16
msun_g = 1.989e33
sigma_convert = msun_g/(molecular_mass*proton_mass_cgs)/(3.086e+18)**2

# P.figure()

def calc_and_dump_sigma_angle(run_ids,output_dirs,snap_strs,outpFile):
    nrun = len(run_ids)
    if not (len(run_ids)==len(output_dirs)==len(snap_strs)):
        raise Exception("All arrays must be same length")

    print(run_ids)
    print(output_dirs)
    print(snap_strs)
    print(outpFile)


    nray_phi = 40
    #nray_phi = 1
    nray_theta = 40
    #nray_theta = 4
    nray = nray_phi*nray_theta
    xyzray = np.zeros((nray,3))

    phis_basis = np.linspace(0.,np.pi/2.,nray_phi)
    thetas_basis = np.linspace(0.,2.*np.pi,nray_theta,endpoint=False)

    phithetas = np.array([[x,y] for x,y in it.product(phis_basis,thetas_basis)]) 
    phis = phithetas[:,0]
    thetas = phithetas[:,1]

    xyzray[:,0] = np.cos(phis)*np.cos(thetas)
    xyzray[:,1] = np.cos(phis)*np.sin(thetas)
    xyzray[:,2] = np.sin(phis)

    outp = np.zeros((nray_phi,nruns*4+1))
    outp[:,0] = phis_basis

    gizmoDir = gizmo_tools.getGizmoDir()

    for irun in range(nruns):
        run_id = run_ids[irun]
        output_dir = output_dirs[irun]
        snap_str = snap_strs[irun]
#         run_id = sys.argv[1+irun*3]
#         output_dir = sys.argv[2+irun*3]
#         snap_str = sys.argv[3+irun*3]
    
        print("Run:",run_id+output_dir+snap_str)

        f = h5py.File(gizmoDir+"/"+run_id+"/"+output_dir+"/snapshot_"+snap_str+".hdf5","r")

        xyz = np.array(f["/PartType0/Coordinates"])
        xyz*=1.e3 # to pc

        m_p = np.array(f["/PartType0/Masses"])
        #m_p*=1.989e+43 # 10^10 solar masses to g
        m_p*=1.e10 # 10^10 solar masses to Msun

        h_p = np.array(f["/PartType0/SmoothingLength"]) # kpc
        h_p*=1.e3 # to pc

        opac_p = np.array(f["/PartType0/AGNOpacity"]) # internal units: kpc**2/1e10 solar mass
        opac_p = 1.e-4 # to pc**2/solar mass


    #     u_p = np.array(f["/PartType0/InternalEnergy"]) # 1e10 erg/g
    #     u_p*=1.e10 # to erg/g
    #     TK_p = (gamma_minus_one/boltzmann_cgs*(molecular_mass*proton_mass_cgs)*u_p)
    #     m_p[TK_p>9000.]=0.

        n = m_p.size

        print("Raytracing")
        sigma_all = sph_plotter.sph_ray_integrate(xyz,m_p,h_p,xyzray,nray,n)
        sigma = np.mean(sigma_all.reshape(-1,nray_theta),axis=1)
        sigma_low = np.min(sigma_all.reshape(-1,nray_theta),axis=1)
        sigma_high = np.max(sigma_all.reshape(-1,nray_theta),axis=1)

        tau_all = sph_plotter.sph_ray_integrate(xyz,m_p*opac_p,h_p,xyzray,nray,n)
        tau = np.mean(tau_all.reshape(-1,nray_theta),axis=1)

    #    phideg = phis_basis/(2.*np.pi)*360

    #     molecular_mass = 4./(1.+3.*.76)
    #     proton_mass_cgs = 1.6726e-24
    #     msun_g = 1.989e33

        sigma*=sigma_convert
        sigma_low*=sigma_convert
        sigma_high*=sigma_convert
    
        outp[:,irun*4+1] = sigma
        outp[:,irun*4+2] = sigma_low
        outp[:,irun*4+3] = sigma_high
        outp[:,irun*4+4] = tau

    #     print("Plotting")
    #     P.plot(phideg,sigma,label=run_id+output_dir+snap_str)
    # 
    # P.legend()
    # P.yscale('log')
    # P.xlabel(r"$\phi$ ($^\degree$)")
    # #P.ylabel(r"$\Sigma$ ($M_\odot$/pc$^2$)")
    # P.ylabel(r"$\Sigma$ (cm$^{-2}$)")
    # P.savefig("../figures/sigma_angle_"+run_id+"_"+output_dir+"_"+snap_str+".png",dpi=200)
    # P.close()

    np.savetxt(outpFile,outp)


if __name__ == '__main__':

    print("Running")

    if len(sys.argv)>1:
        nruns = len(sys.argv)-2
        if ( nruns%3!=0 ):
            raise Exception("python sigma_angle.py folder1 subdir1 dump1 folder2 subdir2 dump2 [etc] outputfile")
        nruns = nruns//3

        outpFile = sys.argv[-1]
    
        run_ids = sys.argv[1:-1:3]
        output_dirs = sys.argv[2:-1:3]
        snap_strs = sys.argv[3:-1:3]
        
        calc_and_dump_sigma_angle(run_ids,output_dirs,snap_strs,outpFile)
    else:
        for snap in ["100","200","500","1000"]:
            outpFile = "data/2014q2redo_sigangle_"+snap+".dat"
            output_dirs=["q2edd05redo","q2edd10_aniso1redo","q2edd10_aniso3redo","q2edd10redo","q2edd20redo","q2redo"]
            nruns = len(output_dirs)
            run_ids = ["2014"]*nruns
            snap_strs = [snap]*nruns

            calc_and_dump_sigma_angle(run_ids,output_dirs,snap_strs,outpFile)

