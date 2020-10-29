#import pynbody
import gizmo_tools
import matplotlib.pyplot as plt
import numpy as np
import pynbody
import functools
from multiprocessing import Pool
import itertools

m_smbh = 1.e6
smbh_soft = 0.01
m_hernquist = 1.e9
a_hernquist = 250.
G_msun_pcyr2_pc2 = 4.5e-15
G_msun_km2s2_pc = 0.0043022682


def accel_smbh_hernquist(m_smbh,smbh_soft,m_hernquist,a_hernquist,r):
    """accel_smbh_hernquist

    Calculates acceleration due to a Plummer softend SMBH plus a Hernquist bulge.

    Mass units are in Msun, length units are in pc, time units are in years
    """
#     G_msun_km2s2_pc = 0.004302052
    x = r/a_hernquist
    r2 = r**2
    m =   m_smbh*r2/(r2+smbh_soft**2) + m_hernquist * (x/(1.+x))**2
    accel = -(G_msun_pcyr2_pc2 * m/r2)
    return accel

# grav_accel = functools.partial(accel_smbh_hernquist,m_smbh,smbh_soft,m_hernquist,a_hernquist)

def t_dyn(r):
#     return np.sqrt(r/-grav_accel(r))
    return np.sqrt(r/-accel_smbh_hernquist(m_smbh,smbh_soft,m_hernquist,a_hernquist,r))

def v_esc(r):
    return np.sqrt(
                2*G_msun_km2s2_pc*(
                    m_hernquist/(a_hernquist+r) +
                    m_smbh/np.sqrt(r**2+smbh_soft**2)
                )
            )
    

run_id = "3032"
# run_name = "longrun_medflow_vesc_defaultaniso_polar_highedd"

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

run_names+= ["longrun_medflow_vesc_defaultaniso_polar_highedd"]

t_factors = np.array([1,2,4,10])

def t_dyn_run(isnap,run_name,verbose=False):
    snap_id = "{:03d}".format(isnap)
    header,snap = gizmo_tools.load_gizmo_nbody(run_id,run_name,snap_id)
    snap['t_dyn'] = pynbody.array.SimArray(t_dyn(snap['r']*.1e3),'yr')
    snap['vesc'] = pynbody.array.SimArray(v_esc(snap['r']*.1e3),'km s**-1')
    time = header['time']*1.e6
    t_scales = time/t_factors
    outp_line = [time]
    outp_line+=[np.sum(snap['t_dyn']<t)/len(snap) for t in t_scales]
    outp_line+=[np.median(time/snap['t_dyn'])]
    outp_line+=[np.sum((snap['t_dyn']<t)|(snap['vr']>snap['vesc']))/len(snap) for t in t_scales]
    if verbose:
        print(outp_line)
        sample_snap = snap[:10]
        print(sample_snap['vr'],sample_snap['vesc'],sample_snap['r'])
    return outp_line

def full_run():
    for run_name in run_names:
        print(run_name)
        snapf = gizmo_tools.lastConsecutiveSnapshot(run_id,run_name,dumpsOrdered=False)

        outp = []


        # outp = map(t_dyn_run,range(snapf))
        with Pool(172) as pool:
            outp = pool.starmap(t_dyn_run,zip(range(snapf),itertools.repeat(run_name)))
        # outp = list(map(t_dyn_run,range(32)))

        # for isnap in range(snapf):
        #     snap_id = "{:03d}".format(isnap)
        #     header,snap = gizmo_tools.load_gizmo_nbody(run_id,run_name,snap_id)
        #     snap['t_dyn'] = pynbody.array.SimArray(t_dyn(snap['r']*.1e3),'yr')
        #     time = header['time']*1.e6
        #     outp+=[[time]+[np.sum(snap['t_dyn']<t)/len(snap) for t in [time,time/2.,time/4.,time/10.]]]
        #     print(outp[-1])

        np.savetxt(f"data/t_dyn_frac_{run_id}_{run_name}.dat",outp)

def test_run():
    run_name = run_names[-1]
    snapf = 10
    for isnap in range(snapf):
        t_dyn_run(isnap,run_name,verbose=True)

# test_run()

full_run()

# print(time,snap['t_dyn'].min(),snap['t_dyn'].max())
# 
# fig,sp = plt.subplots()
# hist,xs, patches = sp.hist(snap['t_dyn'],bins=100,density=True)
# sp.axvline(x=time*1.e6)
# fig.savefig(f"../figures/t_dyn_{run_id}_{run_name}_{snap_id}.png")
