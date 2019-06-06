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

    outp = np.zeros((nray_phi,nruns*6+1))
    outp[:,0] = phis_basis


    for irun in range(nruns):
        run_id = run_ids[irun]
        gizmoDir = gizmo_tools.getGizmoDir(run_id)
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
        opac_p*= 1.e-4 # to pc**2/solar mass


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
        tau_low = np.min(tau_all.reshape(-1,nray_theta),axis=1)
        tau_high = np.max(tau_all.reshape(-1,nray_theta),axis=1)

    #    phideg = phis_basis/(2.*np.pi)*360

    #     molecular_mass = 4./(1.+3.*.76)
    #     proton_mass_cgs = 1.6726e-24
    #     msun_g = 1.989e33

        sigma*=sigma_convert
        sigma_low*=sigma_convert
        sigma_high*=sigma_convert
    
        outp[:,irun*6+1] = sigma
        outp[:,irun*6+2] = sigma_low
        outp[:,irun*6+3] = sigma_high
        outp[:,irun*6+4] = tau
        outp[:,irun*6+5] = tau_low
        outp[:,irun*6+6] = tau_high

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
#         for snap in ["100","200","500","1000"]:
#         run_groups = [["run_a0_e1","run_a1_e1","run_a2_e1","run_a3_e1"],["run_a2_e01_T1000","run_a2_e01_T300","run_a2_e01_T30","run_a2_e01"],["run_a2_e01","run_a2_e05","run_a2_e1","run_a2_e2"],["run_a2_e01","run_a2_e01_SN100","run_a2_e01_SN1000"]]
#         run_groups = [["a0_e1","a1_e1","a2_e1","a3_e1"],["a2_e01_T1000","a2_e01_T300","a2_e01_T30","a2_e01"],["a2_e01","a2_e05","a2_e1","a2_e2"],["a2_e01","a2_e01_SN100","a2_e01_SN1000"]]
#         group_names = ["aniso","mintemp","edd","SN"]
#         run_groups = [["run_a2_e01","run_a2_e01_SN100","run_a2_e01_SN1000"]]
#         group_names = ["SN"]
#         run_groups = [["run_a2_e01","run_a2_e02","run_a2_e05","run_a2_e1","run_a2_e2"]]
#         group_names = ["edd"]
#         run_groups = [["blobrot","blobclose"]]
#         group_names = ["blob"]

#         run_groups = [["SF_test_high_rho_floor_30_thinner_Q2_cut"]]
#         group_names = ["Q2_basis"]
#         snaps = ["025","050","075","100"]

        run_groups = [  ["longrun_weakflow_settled_defaultaniso","longrun_weakflow_vesc_defaultaniso","longrun_weakflow_rapid_defaultaniso"],
                        ["longrun_weakflow_settled_defaultaniso_polar","longrun_weakflow_vesc_defaultaniso_polar","longrun_weakflow_rapid_defaultaniso_polar"]]
        group_names = ["equatorial","polar"]
        
        
        
        snaps = ["100"]
        
        for igroup,group_name in enumerate(group_names):

            output_dirs = run_groups[igroup]
            
            for snap in snaps:
#             for snap in ["100"]:
                outpFile = "data/prodrun_"+group_name+snap+".dat"
#                 output_dirs = ["run_a0_e1","run_a1_e1","run_a2_e01","run_a2_e01_T1000","run_a2_e01_T300","run_a2_e01_T30","run_a2_e05","run_a2_e1","run_a2_e2","run_a3_e1"]
                nruns = len(output_dirs)
#                 run_ids = ["2022"]*nruns
#                 run_ids = ["3001"]*nruns
#                 run_ids = ["2030"]*nruns
                run_ids = ["3032"]*nruns
                snap_strs = [snap]*nruns

                calc_and_dump_sigma_angle(run_ids,output_dirs,snap_strs,outpFile)


#             output_dirs=["q2edd05redo","q2edd10_aniso1redo","q2edd10_aniso3redo","q2edd10redo","q2edd20redo","q2redo"]
#             nruns = len(output_dirs)
#             run_ids = ["2014"]*nruns



