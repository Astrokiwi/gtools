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


def dump_particle_ages(run):
    snapf = gizmo_tools.lastConsecutiveSnapshot(run_id,run,False)
    # birth_table = pd.DataFrame(columns=["id","birth"])
    birth_ids = None
    birth_dates = None

    for isnap in range(snapf+1):
        print(isnap)
        snap_str = f"{isnap:03d}"
        header,snap = gizmo_tools.load_gizmo_pandas(run_id,run,snap_str,["ParticleIDs"])
        if birth_ids is None:
            birth_ids = snap["ParticleIDs"].values
            birth_dates = np.full(len(snap),header["time"])
        else:
            new_ids = ~np.isin(snap["ParticleIDs"],birth_ids)
            new_id_snap = snap[new_ids]
            birth_ids=np.append(birth_ids,new_id_snap["ParticleIDs"].values)
            birth_dates=np.append(birth_dates,np.full(len(new_id_snap),header["time"]))
    np.savetxt("data/particle_ages_{}_{}.dat".format(run_id,run),np.array([birth_ids,birth_dates]).T)
            # birth_table = birth_table.append(   {"id":new_id_snap["ParticleIDs"].values
    #                                             ,"birth":np.full(len(new_id_snap),header["time"])}
    #                                         ,ignore_index=True)
    #         birth_table = birth_table.append(   [new_id_snap["ParticleIDs"].values
    #                                             ,np.full(len(new_id_snap),header["time"])]
    #                                         ,ignore_index=True)


if __name__ == '__main__':
    with Pool(processes=16) as pool:
        pool.map(dump_particle_ages,run_names)
#     for run in run_names:
# #     run="longrun_weakflow_vesc_defaultaniso"
# 
#         dump_particle_ages(run_id,run)
#         