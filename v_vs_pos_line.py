import numpy as np

import gizmo_tools

from sys import path,exit
path.append("src/")
import tab_interp

import sys
sys.path.append("visualisation/")

from sph_plotter import sph_plotter

import matplotlib.pyplot as plt
import numpy as np

# import time

from multiprocessing import Pool
from functools import partial


gizmo_tools.derive_opacities()

run_id = "3032"
# run_names = [
# "longrun_medflow_settled_defaultaniso_polar",
# "longrun_medflow_vesc_defaultaniso",
# "longrun_medflow_vesc_defaultaniso_polar",
# "longrun_weakflow_rapid_defaultaniso",
# "longrun_weakflow_rapid_defaultaniso_polar",
# "longrun_weakflow_settled_defaultaniso",
# "longrun_weakflow_settled_defaultaniso_polar",
# "longrun_weakflow_vesc_defaultaniso",
# "longrun_weakflow_vesc_defaultaniso_polar",
# "newflow_settled_thin_up",
# "newflow_vesc_thin_45",
# "newflow_vesc_thin_side",
# "newflow_vesc_thin_up"]

run_names = [
"longrun_medflow_vesc_defaultaniso_polar",
"newflow_vesc_thin_45"
]
# run_names = ["longrun_medflow_vesc_defaultaniso_polar"]

nruns = len(run_names)
snap_str = "200"

nbins = 400
# broaden_value = 1.

lines = ["co1","co2","hcn1","hcn2","h2_1","h2_2","h2_3"]
# lines = ["co1"]
# lines = ["h2_1","co1","hcn1"]
line_codes = ["line_"+line for line in lines]


tableDate="060319"
tableRes="0.1"
cloudy_table = gizmo_tools.cloudy_table(tableDate,tableRes)


vmin=-200.
# vmin=-100.
vmax=-vmin


# rotate
# rotate_phi = 20.



# nray_strip = 41
# nray = nray_strip*2
# ray_dir_z = np.array([0.,0.,1.])
# ray_dir_y = np.array([0.,1.,0.])
# ray_steps = np.linspace(-20.,20.,nray_strip)
# 
# ray_xyz_z = np.zeros((nray_strip,3))
# ray_xyz_z[:,1] = ray_steps
# ray_xyz_y = np.zeros((nray_strip,3))
# ray_xyz_y[:,2] = ray_steps
# ray_dirs_z = np.broadcast_to(ray_dir_z,(nray_strip,3))
# ray_dirs_y = np.broadcast_to(ray_dir_y,(nray_strip,3))
# 


# ray_dirs = np.vstack((ray_dirs_z,ray_dirs_y))

def v_vs_pos_rays(run_id,run_name,snap_str,rotate_phi_deg):

    rotate_phi = rotate_phi_deg/360.*2.*np.pi
    print(f"Running with angle {rotate_phi_deg}")

    nray_strip = 41
    nray = nray_strip*2
#     ray_dir_z = np.array([0.,np.sin(rotate_phi),np.cos(rotate_phi)])
    ray_dir_y = np.array([0.,np.cos(rotate_phi),np.sin(rotate_phi)])
    ray_dir_z = ray_dir_y
    ray_steps = np.linspace(-20.,20.,nray_strip)
#     ray_steps = np.linspace(-.5,.5,nray_strip)

    ray_xyz_z = np.zeros((nray_strip,3))
    ray_xyz_z[:,1] = -ray_steps*np.sin(rotate_phi)
    ray_xyz_z[:,2] = ray_steps*np.cos(rotate_phi)
    ray_xyz_y = np.zeros((nray_strip,3))
    ray_xyz_y[:,0] = ray_steps
    ray_dirs_z = np.broadcast_to(ray_dir_z,(nray_strip,3))
    ray_dirs_y = np.broadcast_to(ray_dir_y,(nray_strip,3))
    
    print(rotate_phi_deg,rotate_phi/(np.pi),ray_dir_z,ray_xyz_z[0],ray_dir_y,ray_xyz_y[0])

#     start = time.time()
    header,particles = gizmo_tools.load_gizmo_pandas(run_id,run_name,snap_str,
        ["Masses","Coordinates","Velocities","AGNIntensity","AGNDepth","Density","InternalEnergy","SmoothingLength"])
    gizmo_tools.load_calc_reqs(particles,["nH","temp","rad3d"])
    print("Getting dust etc information")
    cloudy_table.interp(particles)
#     duration = time.time()-start
#     print("Loading {} time:".format(run_name),duration)

    particles["SoundSpeed"] = np.sqrt(10./9.*particles["InternalEnergy"]) # in cm/s
    particles["SoundSpeed"]*=1.e-5 # to km/s
#     particles["SoundSpeed"] = 4.#np.sqrt(10./9.*particles["InternalEnergy"])
    xyz = np.array((particles['Coordinates_x'],particles['Coordinates_y'],particles['Coordinates_z'])).T
    vel = np.array((particles['Velocities_x'],particles['Velocities_y'],particles['Velocities_z'])).T
    n = len(particles)
    mask = np.ones(n,dtype=bool)


#     for line_code in ["line_hcn1"]:
#     for line_code in line_codes:
    for line in lines:
        line_code = "line_"+line
        line_emission = particles[line_code]*particles["Masses"]*gizmo_tools.msun_g
        line_cross_sections = gizmo_tools.line_opacities[line]*particles["Masses"]
        
        broaden = np.sqrt(particles["SoundSpeed"]**2+0**2) # add some extra broadening

    #     start = time.time()
    #     rayhists = sph_plotter.sph_ray_histogram_2(xyz,line_emission,particles["SmoothingLength"],particles['Velocities_z'],vmin,vmax,ray_dirs,ray_xyzs,mask,particles["SoundSpeed"],nbins,False,nray,n)
    #     duration = time.time()-start
    #     print("Method 2:",duration)
#         start = time.time()
#         rayhists_z = sph_plotter.sph_ray_histogram(xyz,line_emission,particles["SmoothingLength"],particles['Velocities_z'],vmin,vmax,ray_dirs_z,ray_xyz_z,mask,particles["SoundSpeed"],nbins,False,nray_strip,n)
#         rayhists_y = sph_plotter.sph_ray_histogram(xyz,line_emission,particles["SmoothingLength"],particles['Velocities_y'],vmin,vmax,ray_dirs_y,ray_xyz_y,mask,particles["SoundSpeed"],nbins,False,nray_strip,n)
#         rayhists = np.hstack([rayhists_z,rayhists_y])
#         duration = time.time()-start
#         print("Integration time:",duration)
# 
#         start = time.time()
#         np.savetxt("data/line_profs_{}_{}.dat".format(run_name,line_code),rayhists)
#         duration = time.time()-start
#         print("Saving time:",duration)
        
#         line_emission[:] = np.max(line_emission) # UNIFORM
#         print("UNIFORM: FOR TEST")
        print(ray_dirs_z.shape,vel.shape)
        vlos = np.sum(ray_dirs_y[0,:]*vel,axis=1)
        print(np.max(vlos),np.min(vlos))
        print(np.max(vel),np.min(vel))
        print(vlos)
        
        
        print(rotate_phi,"Starting integration with extinction")
#         start = time.time()
        rayhists_z = sph_plotter.sph_ray_histogram_opacity(xyz,line_emission,particles["SmoothingLength"],vel,vmin,vmax,line_cross_sections,ray_dirs_z,ray_xyz_z,mask,broaden,nbins,nray_strip,n)
        rayhists_y = sph_plotter.sph_ray_histogram_opacity(xyz,line_emission,particles["SmoothingLength"],vel,vmin,vmax,line_cross_sections,ray_dirs_y,ray_xyz_y,mask,broaden,nbins,nray_strip,n)
        rayhists = np.hstack([rayhists_z,rayhists_y])
#         duration = time.time()-start
#         print("Ext Integration time:",duration)

#         start = time.time()
        print(rotate_phi,"Saving with extinction")
        np.savetxt("data/line_profs_ext_{}_{}_{}deg_{}.dat".format(run_name,line_code,rotate_phi_deg,snap_str),rayhists)

#         duration = time.time()-start
#         print("Saving time:",duration)
#         sys.exit()

#         start = time.time()
        print(rotate_phi,"Starting integration without extinction")
        rayhists_z = sph_plotter.sph_ray_histogram(xyz,line_emission,particles["SmoothingLength"],vel,vmin,vmax,ray_dirs_z,ray_xyz_z,mask,broaden,nbins,False,nray_strip,n)
        rayhists_y = sph_plotter.sph_ray_histogram(xyz,line_emission,particles["SmoothingLength"],vel,vmin,vmax,ray_dirs_y,ray_xyz_y,mask,broaden,nbins,False,nray_strip,n)
        rayhists = np.hstack([rayhists_z,rayhists_y])
#         duration = time.time()-start
#         print("Integration time:",duration)

#         start = time.time()
        print(rotate_phi,"Saving without extinction")
        np.savetxt("data/line_profs_{}_{}_{}deg_{}.dat".format(run_name,line_code,rotate_phi_deg,snap_str),rayhists)
#         duration = time.time()-start
#         print("Saving time:",duration)



# angles = np.arange(0.,90.,5.)
# angles = np.arange(0.,25.,5.)
# angles = [0.,10.,20.]
angles = [10.]
for run_name in run_names:
    v_vs_pos_rays_for_pool = partial(v_vs_pos_rays,run_id,run_name,snap_str)

    # serial, for test (or like OMP_NUM_THREADS)
    for angle in angles:
        v_vs_pos_rays_for_pool(angle)
        
#     with Pool(processes=3) as pool:
#         pool.map(v_vs_pos_rays_for_pool,angles)
    




#     
#     print("Method 1:",duration)
#     
#     
#     v_values = np.arange(vmin,vmax,(vmax-vmin)/nbins)
#     plt.figure()
#     for iray in range(nray):
#         plt.plot(v_values,rayhists[:,iray])
#     plt.yscale('log')
#     plt.ylim([np.percentile(rayhists[rayhists>0.],30.),None])
# #     plt.legend()
#     plt.savefig("pics/rayhist_test.pdf")
#     plt.close('all')
#     