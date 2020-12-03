#!/usr/bin/env python3

import sys
import numpy as np
from . import gizmo_tools

file_header = """"# Gas particles
# SKIRT 9 import format for a medium source using M_dust = f_dust x Z x M_gas
#
# Column 1: x-coordinate (pc)
# Column 2: y-coordinate (pc)
# Column 3: z-coordinate (pc)
# Column 4: smoothing length (pc)
# Column 5: gas mass (Msun)
# Column 6: temperature (K)
#"""

# Column 5: x-velocity (km/s)
# Column 6: y-velocity (km/s)
# Column 7: z-velocity (km/s)
# Column 6: metallicity (1)


if __name__ == '__main__' :
    run_dir = sys.argv[1]
    run_id = sys.argv[2]
    snap_str = sys.argv[3]

    try:
        outfile = sys.argv[4]
    except:
        outfile = None # print to stream

    header,snap = gizmo_tools.load_gizmo_nbody(run_dir,run_id,snap_str)

    gizmo_tools.nbody_calc_val(snap,'temp')

    out_table = np.vstack([
                            snap.gas['x'].in_units('pc')
                            ,snap.gas['y'].in_units('pc')
                            ,snap.gas['z'].in_units('pc')
                            ,snap.gas['smooth'].in_units('pc')
                            # ,snap.gas['vel'][:,0].in_units('km s**-1')
                            # ,snap.gas['vel'][:,1].in_units('km s**-1')
                            # ,snap.gas['vel'][:,2].in_units('km s**-1')
                            ,snap.gas['mass'].in_units('Msol')
                            # ,np.ones(len(snap.gas)) # dust
                            ,snap.gas['temp'].in_units('K')
        ]).T

    np.savetxt(outfile,out_table,header=file_header)
