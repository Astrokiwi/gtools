import numpy as np

import gizmo_tools

from sys import path,exit
path.append("src/")
import tab_interp

from sys import path
path.append("visualisation/")
from sph_plotter import sph_plotter

import matplotlib.pyplot as plt

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

nbins = 200
broaden_value = 1.

lines = ["co1","co2","hcn1","hcn2","h2_1","h2_2","h2_3"]
line_codes = ["line_"+line for line in lines]


tableDate="060319"
tableRes="0.1"
cloudy_table = gizmo_tools.cloudy_table(tableDate,tableRes)

for run_name in ["longrun_weakflow_vesc_defaultaniso_polar"]:
    header,particles = gizmo_tools.load_gizmo_pandas(run_id,run_name,snap_str,
        ["Masses","Coordinates","Velocities","AGNIntensity","AGNDepth","Density","InternalEnergy","SmoothingLength"])
    gizmo_tools.load_calc_reqs(particles,["nH","temp","rad3d"])
    cloudy_table.interp(particles)
    
    
    nray = 3
    ray_dir = np.array([0.,0.,1.])
#     ray_z = np.linspace(-10.,10.,nray)
    ray_y = np.linspace(-2.,2.,nray)
    ray_xyzs = np.zeros((nray,3))
    ray_xyzs[:,1] = ray_y
    
    vmin=-420
    vmax=-vmin
    
    ray_dirs = np.broadcast_to(ray_dir,ray_xyzs.shape)
    
    particles["SoundSpeed"] = np.sqrt(10./9.*particles["InternalEnergy"])
    
    xyz = np.array((particles['Coordinates_x'],particles['Coordinates_y'],particles['Coordinates_z'])).T
    line_emission = particles["line_hcn1"]*particles["Masses"]*gizmo_tools.msun_g                                                                                   
    n = len(particles)
    mask = np.ones(n,dtype=bool)
#     broaden = np.full(n,broaden_value)
    rayhists = sph_plotter.sph_ray_histogram(xyz,line_emission,particles["SmoothingLength"],particles['Velocities_z'],vmin,vmax,ray_dirs,ray_xyzs,mask,particles["SoundSpeed"],nbins,False,nray,n)
    
    v_values = np.arange(vmin,vmax,(vmax-vmin)/nbins)
    plt.figure()
    for iray in range(nray):
        plt.plot(v_values,rayhists[:,iray],label=ray_y[iray])
    plt.yscale('log')
    plt.ylim([np.percentile(rayhists[rayhists>0.],10.),None])
    plt.legend()
    plt.savefig("pics/rayhist_test.pdf")
    plt.close('all')
    