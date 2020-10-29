import numpy as np
import gizmo_tools
from multiprocessing import Pool
from functools import partial

run_id="3032"

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

T_threshold = 1.e2

rmax = 20.

def summarise_dump(run_id,run_name,snap_str):
    header,particles = gizmo_tools.load_gizmo_pandas(run_id,run_name,snap_str,
        ["Masses","Coordinates","Velocities","InternalEnergy"],False)
    print(header["time"])
    gizmo_tools.load_calc_reqs(particles,["temp","rad3d","v_circ","v_rad"])

    warm_particles = particles[(particles.temp>T_threshold) & (particles.rad3d<rmax)]
    cool_particles = particles[(particles.temp<=T_threshold) & (particles.rad3d<rmax)]
    
    
    out_line = [header["time"]]
    for p in [cool_particles,warm_particles]:
        for val in ["rad3d","v_circ","v_rad"]:
            out_line.append(p[val].mean())
    return out_line

def summarise_evolution(run_id,run_name,nprocs=64):
    snapf = gizmo_tools.lastConsecutiveSnapshot(run_id,run_name,dumpsOrdered=False)
    snap_strs = ["{:03d}".format(x) for x in range(snapf)]
    
    with Pool(processes=nprocs) as pool:
        one_run_summarise_dump = partial(summarise_dump,run_id,run_name)
        out_tab = pool.map(one_run_summarise_dump,snap_strs)

    np.savetxt(
                f"data/components_summary_{run_id}_{run_name}.dat"
                ,out_tab
                )

if __name__=='__main__':
#     with Pool(processes=3) as pool:
#         summarise_evolution_this_id = partial(summarise_evolution,run_id)
#         pool.map(summarise_evolution_this_id,run_names)    

    for run_name in run_names:
        summarise_evolution(run_id,run_name,nprocs=176)