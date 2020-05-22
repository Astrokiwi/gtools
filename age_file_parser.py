import gizmo_tools
import pandas as pd
import numpy as np
from multiprocessing import Pool

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

run_id="3032"


def parse_particle_ages(run):
    snapf = gizmo_tools.lastConsecutiveSnapshot(run_id,run,False)
    
    age_file = "data/particle_ages_{}_{}.dat".format(run_id,run)
    print(age_file)
    id_table,age_table = np.loadtxt(age_file,unpack=True)
    id_table=id_table.astype(np.int)
    sorted_indices = np.argsort(id_table)
    id_table=id_table[sorted_indices]
    age_table=age_table[sorted_indices]
    print((id_table == np.sort(id_table)).all())

    for isnap in range(snapf+1):
        print(run,isnap,snapf)
        snap_str = f"{isnap:03d}"
        header,snap = gizmo_tools.load_gizmo_pandas(run_id,run,snap_str,["ParticleIDs"])
        id_indices = np.searchsorted(id_table,snap["ParticleIDs"])
        print(np.min(id_indices),np.max(id_indices),age_table.size)
        ages = age_table[id_indices]
        np.savetxt("data/age_{}_{}_{}.dat".format(run_id,run,isnap),ages)


if __name__ == '__main__':
    with Pool(processes=16) as pool:
        pool.map(parse_particle_ages,run_names)

#     for run in run_names:
#     run="longrun_weakflow_vesc_defaultaniso"
#         parse_particle_ages(run_id,run)