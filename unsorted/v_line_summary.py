import numpy as np

import gizmo_tools

from sys import path,exit
path.append("../src/")
import tab_interp

run_id = "3032"
run_names = [
# "longrun_medflow_settled_defaultaniso",
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
max_rad = 100.
min_rad = 1.5



lines = ["co1","co2","hcn1","hcn2","h2_1","h2_2","h2_3"]
line_codes = ["line_"+line for line in lines]

grain_size = 1. # microns
grain_density = 3 # g/cm**3
grain_specific_emission_cross_section = 3./(4.*grain_size*1.e-6*grain_density) # cm**2/g, = 250000 !
grain_specific_emission_cross_section_astronomical = grain_specific_emission_cross_section * 0.000208908219 # to pc**2/Msun

tableDate="060319"
tableRes="0.1"
cloudy_table = gizmo_tools.cloudy_table(tableDate,tableRes)



for irun,run_name in enumerate(run_names):
# for irun,run_name in enumerate(["newflow_vesc_thin_45"]):
    print(irun,run_name)
    header,particles = gizmo_tools.load_gizmo_pandas(run_id,run_name,snap_str,
        ["Masses","Coordinates","Velocities","AGNIntensity","AGNDepth","Density","InternalEnergy","SmoothingLength"])
    gizmo_tools.load_calc_reqs(particles,["nH","temp","v_circ","v_rad"])

    print("tabulating values")
    cloudy_table.interp(particles)
    particles["brightness"] = 5.67e-5 * particles["dustT"]**4. * particles["dg"]/np.nanmax(particles["dg"]) # erg/s/cm^2
    particles["emitArea"] = np.minimum(particles["dg"]*particles["Masses"]*grain_specific_emission_cross_section_astronomical,
                                            np.pi*particles["SmoothingLength"]**2)
    particles["luminosity"] = particles["emitArea"]*particles["brightness"]
    for line in lines:
        particles["lum_"+line] = particles["line_"+line]*particles["Masses"] # n.b. units irrelevant, we are normalising later
        
    summary_outp = []
    for brightness_key in ["luminosity"]+["lum_"+line for line in lines]:
        summary_outp.append([np.average(particles[v_key],weights=particles[brightness_key])
                        for v_key in
                        ("v_circ","v_rad")])
    np.savetxt(f"data/v_summary_{run_name}.dat",summary_outp)
    
