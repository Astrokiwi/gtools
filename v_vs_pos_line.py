import numpy as np

import gizmo_tools

from sys import path,exit
path.append("src/")
import tab_interp

from sys import path
path.append("visualisation/")
from sph_plotter import sph_plotter

import matplotlib.pyplot as plt


import time

run_id = "3032"
run_names = [
"longrun_medflow_settled_defaultaniso_polar",
"longrun_medflow_vesc_defaultaniso",
"longrun_medflow_vesc_defaultaniso_polar",
"longrun_weakflow_rapid_defaultaniso",
"longrun_weakflow_rapid_defaultaniso_polar",
"longrun_weakflow_settled_defaultaniso",
"longrun_weakflow_settled_defaultaniso_polar",
"longrun_weakflow_vesc_defaultaniso",
"longrun_weakflow_vesc_defaultaniso_polar",
"newflow_settled_thin_up",
"newflow_vesc_thin_45",
"newflow_vesc_thin_side",
"newflow_vesc_thin_up"]

nruns = len(run_names)
snap_str = "100"

nbins = 400
broaden_value = 1.

lines = ["co1","co2","hcn1","hcn2","h2_1","h2_2","h2_3"]
line_codes = ["line_"+line for line in lines]


tableDate="060319"
tableRes="0.1"
cloudy_table = gizmo_tools.cloudy_table(tableDate,tableRes)


nray_strip = 41
nray = nray_strip*2
ray_dir_z = np.array([0.,0.,1.])
ray_dir_y = np.array([0.,1.,0.])
#     ray_z = np.linspace(-10.,10.,nray)
ray_steps = np.linspace(-20.,20.,nray_strip)
# ray_y = np.linspace(-2.,2.,nray)
# ray_xyzs = np.zeros((nray,3))
# # ray_xyzs[:,1] = ray_y
# ray_xyzs[0:nray_strip,1] = ray_steps
# ray_xyzs[nray_strip:nray,2] = ray_steps

ray_xyz_z = np.zeros((nray_strip,3))
ray_xyz_z[:,1] = ray_steps
ray_xyz_y = np.zeros((nray_strip,3))
ray_xyz_y[:,2] = ray_steps
ray_dirs_z = np.broadcast_to(ray_dir_z,(nray_strip,3))
ray_dirs_y = np.broadcast_to(ray_dir_y,(nray_strip,3))

vmin=-1000.
vmax=-vmin

# ray_dirs = np.vstack((ray_dirs_z,ray_dirs_y))

# for run_name in ["longrun_weakflow_vesc_defaultaniso_polar"]:
for run_name in run_names:
    start = time.time()
    header,particles = gizmo_tools.load_gizmo_pandas(run_id,run_name,snap_str,
        ["Masses","Coordinates","Velocities","AGNIntensity","AGNDepth","Density","InternalEnergy","SmoothingLength"])
    gizmo_tools.load_calc_reqs(particles,["nH","temp","rad3d"])
    cloudy_table.interp(particles)
    duration = time.time()-start
    print("Loading {} time:".format(run_name),duration)

    particles["SoundSpeed"] = np.sqrt(10./9.*particles["InternalEnergy"])
#     particles["SoundSpeed"] = 4.#np.sqrt(10./9.*particles["InternalEnergy"])
    xyz = np.array((particles['Coordinates_x'],particles['Coordinates_y'],particles['Coordinates_z'])).T
    n = len(particles)
    mask = np.ones(n,dtype=bool)


#     for line_code in ["line_hcn1"]:
    for line_code in line_codes:
        line_emission = particles[line_code]*particles["Masses"]*gizmo_tools.msun_g                                                                                   

    #     start = time.time()
    #     rayhists = sph_plotter.sph_ray_histogram_2(xyz,line_emission,particles["SmoothingLength"],particles['Velocities_z'],vmin,vmax,ray_dirs,ray_xyzs,mask,particles["SoundSpeed"],nbins,False,nray,n)
    #     duration = time.time()-start
    #     print("Method 2:",duration)

        start = time.time()
        rayhists_z = sph_plotter.sph_ray_histogram(xyz,line_emission,particles["SmoothingLength"],particles['Velocities_z'],vmin,vmax,ray_dirs_z,ray_xyz_z,mask,particles["SoundSpeed"],nbins,False,nray_strip,n)
        rayhists_y = sph_plotter.sph_ray_histogram(xyz,line_emission,particles["SmoothingLength"],particles['Velocities_y'],vmin,vmax,ray_dirs_y,ray_xyz_y,mask,particles["SoundSpeed"],nbins,False,nray_strip,n)
        rayhists = np.hstack([rayhists_z,rayhists_y])
        duration = time.time()-start
        print("Integration time:",duration)

        start = time.time()
        np.savetxt("data/line_profs_{}_{}.dat".format(run_name,line_code),rayhists)
        duration = time.time()-start
        print("Saving time:",duration)
    
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