import pandas as pd
import numpy as np

base_dir = "/srv/djw1g16/gizmos/3032"

runs = [
"longrun_medflow_settled_defaultaniso",
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

parameters_tofind = ["outflowRate","outflowThetaCentre","outflowThetaWidth","outflowVel"]

run_parameters = {}

for run in runs:
    run_prams = dict(zip(parameters_tofind,[0]*len(parameters_tofind)))
    pram_file_name = "{}/{}.param".format(base_dir,run)
    with open(pram_file_name) as f:
        for line in f:
            for param in parameters_tofind:
                if param in line:
                    v = float(line.split()[1])
                    run_prams[param] = v
    run_parameters[run] = run_prams
#     print(pram_file_name,run_prams)


# demonstration
for run in runs:
    run_prams = run_parameters[run]
    theta_min = run_prams["outflowThetaCentre"]-run_prams["outflowThetaWidth"]/2.
    theta_max = run_prams["outflowThetaCentre"]+run_prams["outflowThetaWidth"]/2.
    surface_covered = (np.cos(np.radians(theta_min))-np.cos(np.radians(theta_max)))/2. # as a fraction of total surface
    steradians_covered = surface_covered*4.*np.pi
    run_prams["solid_angle"]=steradians_covered
    # outflow_per_steradian = run_prams["outflowRate"]/steradians_covered
    print("{:48}{:6.1f}{:6.1f}{:6.3f}{:6.3f}{:6.3f}".format(run,theta_min,theta_max,surface_covered,run_prams["outflowRate"],outflow_per_steradian))

