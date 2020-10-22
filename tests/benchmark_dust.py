import numpy as np

import gizmo_tools

from sys import path
path.append("src/")
import tab_interp

import time

import ctypes

run_id = "3032"
run_name = "longrun_medflow_settled_defaultaniso_polar"
snap_str = "100"

tableDate="281118"
tableRes="0.0001"
chTab = tab_interp.CoolHeatTab( ("coolheat_tab_marta/shrunk_table_labels_"+tableDate+"tau.dat"),
                                ("coolheat_tab_marta/shrunk_table_"+tableDate+"_m"+tableRes+"_hsmooth_tau.dat"),
                                ("coolheat_tab_marta/shrunk_table_labels_"+tableDate+"taunodust.dat"),
                                ("coolheat_tab_marta/shrunk_table_"+tableDate+"_m"+tableRes+"_hsmooth_taunodust.dat"),
                                ("coolheat_tab_marta/shrunk_table_labels_"+tableDate+"taudense.dat"),
                                ("coolheat_tab_marta/shrunk_table_"+tableDate+"_m"+tableRes+"_hsmooth_taudense.dat")
                                )


interpTabVec = np.vectorize(chTab.interpTab)


header,particles = gizmo_tools.load_gizmo_pandas(run_id,run_name,snap_str,
    ["Masses","Coordinates","Velocities","AGNIntensity","AGNDepth","Density","InternalEnergy","SmoothingLength"])
gizmo_tools.load_calc_reqs(particles,["nH","temp","rad3d"])


def old_version(particles):
    # this part is slow
    start = time.time()
    tabStructs = interpTabVec(  particles["nH"].astype(np.float64),
                                particles["temp"].astype(np.float64),
                                particles["AGNIntensity"].astype(np.float64),
                                particles["AGNDepth"].astype(np.float64))
    end = time.time()
    print("old table lookup:",end-start)

    start = time.time()
    dustTemp = np.array(list(map(lambda y: y.dustT, tabStructs)))
    dg = np.array(list(map(lambda y: y.dg, tabStructs)))
    end = time.time()
    print("old assignment:",end-start)

    return dustTemp,dg

def double_pointer_to_array(x,n):
    ptr = int(x)
    type_size = ctypes.c_double * n
    base_array = type_size.from_address(ptr)
    numpy_array = np.array(base_array)
    return numpy_array

def new_version(particles):
    #tabStructs = 
    start = time.time()
    n = len(particles)
    tabStructs=chTab.interpTabArray( 
                           particles["nH"].astype(np.float64),
                           particles["temp"].astype(np.float64),
                           particles["AGNIntensity"].astype(np.float64),
                           particles["AGNDepth"].astype(np.float64)
                           )
    end = time.time()
    print("new table lookup:",end-start)

    start = time.time()
#     tabStructs = tab_interp.coolHeatDustArray.frompointer(x)
#     dustTemp = np.empty((n))
#     dg = np.empty((n))
#     for i in range(n):
#         dg[i] = tabStructs[i].dg
#         dustTemp[i] = tabStructs[i].dustT
    dustTemp = double_pointer_to_array(tabStructs.dustT,n)                                                                                                
    dg =  double_pointer_to_array(tabStructs.dg,n)
    end = time.time()
    print("new assignment:",end-start)

    return dustTemp,dg

start = time.time()
T_new,dg_new = new_version(particles)
end = time.time()
print("new:",end-start)


start = time.time()
T_old,dg_old = old_version(particles)
end = time.time()
print("old:",end-start)


