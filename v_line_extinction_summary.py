import numpy as np

import gizmo_tools

import sys
sys.path.append("src/")
sys.path.append("visualisation/")
import tab_interp
from sph_plotter import sph_plotter

# sph_plotter.set_parallel()

from multiprocessing import Pool


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
nprocs = len(run_names)

# run_names = [
# "newflow_vesc_thin_45"]

nruns = len(run_names)
snap_str = "100"

lines = ["co1","co2","hcn1","hcn2","h2_1","h2_2","h2_3"]
# lines = ["co1","h2_1"]
# lines = ["h2_1"]
line_codes = ["line_"+line for line in lines]

lums = lines+["dust"]

grain_size = 1. # microns
grain_density = 3 # g/cm**3
grain_specific_emission_cross_section = 3./(4.*grain_size*1.e-6*grain_density) # cm**2/g, = 250000 !
grain_specific_emission_cross_section_astronomical = grain_specific_emission_cross_section * 0.000208908219 # to pc**2/Msun
IR_opacity = 1. * 0.000208908219 # 1 cm**2/g, somewhat arbitrary

tableDate="060319"
tableRes="0.1"
cloudy_table = gizmo_tools.cloudy_table(tableDate,tableRes)

gizmo_tools.derive_opacities()

taugrid_L=1024
# taugrid_L=1
taugrid_width = 200. # pc
taugrid_corner = [-taugrid_width/2,-taugrid_width/2]

def weighted_mean(values,weights=1):
    return np.sum(values*weights)/np.sum(weights)

def calc_dump_v_line_extinction(run_name):
# for irun,run_name in enumerate(["newflow_vesc_thin_45"]):
    print("Loading",run_name)
    header,particles = gizmo_tools.load_gizmo_pandas(run_id,run_name,snap_str,
        ["Masses","Coordinates","Velocities","AGNIntensity","AGNDepth","Density","InternalEnergy","SmoothingLength"])
    gizmo_tools.load_calc_reqs(particles,["nH","temp","v_circ","v_rad"])
    
    print("n=",len(particles))

    print("tabulating values")
    cloudy_table.interp(particles)
    particles["brightness"] = 5.67e-5 * particles["dustT"]**4. * particles["dg"]/np.nanmax(particles["dg"]) # erg/s/cm^2
    particles["emitArea"] = np.minimum(particles["dg"]*particles["Masses"]*grain_specific_emission_cross_section_astronomical,
                                            np.pi*particles["SmoothingLength"]**2)
    particles["lum_dust"] = particles["emitArea"]*particles["brightness"]/particles["Masses"]
    particles["opac_dust"] = IR_opacity

    print("calculating optical depths")
    flagarray = np.ones(len(particles),dtype=np.bool)
#   doing face-on LOS right now
#     particles["vlos"] = particles["Velocities_z"]

    for line in lines:
        particles["opac_"+line] = gizmo_tools.line_opacities[line]
        particles["lum_"+line] = particles["line_"+line]*particles["Masses"]

    for lum_code in lums:
        particles["tau_"+lum_code]=np.zeros(len(particles))
        
#         print("mass grid")
#         g = sph_plotter.sph_dense(particles["Coordinates_x"],particles["Coordinates_y"],particles["Masses"],particles["SmoothingLength"],taugrid_L,taugrid_corner,taugrid_width,flagarray)

        for (direction_name,directions) in [("face",["x","y","z"]),("edge",["x","z","y"])]:
            print(f"{run_name} opacity grid - {lum_code}, {direction_name}-on")
            opac_grid = sph_plotter.sph_mean_optical_depth(
                            particles[f"Coordinates_{directions[0]}"],particles[f"Coordinates_{directions[1]}"],particles[f"Coordinates_{directions[2]}"]
                            ,particles["Masses"],particles["SmoothingLength"]
                            ,particles["opac_"+lum_code]
                            ,taugrid_L,taugrid_corner,taugrid_width,flagarray,particles["tau_"+lum_code])
            particles["lumext_"+direction_name+"_"+lum_code]=particles["lum_"+lum_code]*np.exp(-particles["tau_"+lum_code])
                
    
#     particles.to_csv("data/line_ext_v.txt",sep=' ')
    
# 
    summary_outp = []
    for brightness_key in   ["lum_"+lum_code for lum_code in lums]\
                            +["lumext_"+direction_name+"_"+lum_code for lum_code in lums for direction_name in ["face","edge"]]:

        mean_v = [weighted_mean(particles[f"Velocities_{direction}"],weights=particles[brightness_key]) for direction in ["y","z"]]
        summary_outp.append([weighted_mean(particles[v_key],weights=particles[brightness_key])
                        for v_key in ("v_circ","v_rad")] + \
                        [np.sqrt(weighted_mean((particles[f"Velocities_{direction}"]-mean_vlos)**2,weights=particles[brightness_key])) for\
                                direction,mean_vlos in zip(["y","z"],mean_v)]
                        )
    np.savetxt(f"data/v_summary_{run_name}.dat",summary_outp)

with Pool(nprocs) as pool:
    pool.map(calc_dump_v_line_extinction,run_names)

# for irun,run_name in enumerate(run_names):
#     calc_dump_v_line_extinction(run_name)

    
