import numpy as np

import gizmo_tools

from sys import path,exit

snap_str = "100"

nbins = 200

for run_name in ["longrun_weakflow_vesc_defaultaniso_polar"]:
    header,particles = gizmo_tools.load_gizmo_pandas("3032",run_name,snap_str,
        ["Masses","Coordinates","Velocities","AGNIntensity","AGNDepth","InternalEnergy","SmoothingLength"])

    outp_strings = ["Masses","Coordinates_x","Coordinates_y","Coordinates_z","Velocities_x","Velocities_y","Velocities_z","SmoothingLength"]
    outp = [particles[outp_string].values for outp_string in outp_strings]
    outp = np.array(outp).T
    
    np.savetxt("../quickdump.dat",outp)