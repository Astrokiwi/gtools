print("Importing")

# Import libraries to do our magic
import numpy as np
import sys
import os
import h5py

import re

import gizmo_tools

if __name__ == '__main__':
    print("Running")

    run_id = sys.argv[1]
    output_dir = sys.argv[2]
    i_track = int(sys.argv[3])

    if len(sys.argv)>4:
        snapf = int(sys.argv[4])
    else:
        snapf = 0

    snapi = 0

    if snapf==0:
        snapf = gizmo_tools.lastConsecutiveSnapshot(run_id,output_dir)


    all_outp = []

    track_vals = ["x"
                 ,"y"
                 ,"z"
                 ,"r"
                 ,"vr"
                 ,"AGNDepth"
                 ,"AGNDepth_2"
                 ,"AGNDepth_Effective"
                 ,"AGNIntensity"
                 ,"AGNIntensity_2"
                 ,"AGNIntensity_Effective"
                 ,"AGNOpacity"
                 ,"RadiativeAcceleration_r"
                 ,"AGNHeat"
                 ,"nH"
                 ,"TimeStep"
                    ]

    for snapx in range(snapi,snapf+1):
        print(snapx)
        header,snap = gizmo_tools.load_gizmo_nbody(run_id,output_dir,f"{snapx:03d}",load_vals=["temp","nH","RadiativeAcceleration_r"])
        time = header["time"]

        this_p = np.argwhere(snap["iord"]==i_track)[0][0]
    
        outp = [time]
        for key in track_vals:
            outp+=[snap[key][this_p]]
            
#         print(outp)
    
        all_outp+=[outp]
        
#     print(all_outp)

    all_outp = np.array(all_outp)

    np.savetxt("data/trackp"+run_id+output_dir+"_"+str(i_track)+".dat",all_outp)





