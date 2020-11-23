# Import libraries to do our magic

import os
this_dir, this_filename = os.path.split(__file__)


import numpy as np

from ..tools import gizmo_tools
from multiprocessing import Pool

from scipy import linalg
import functools

import pynbody

import matplotlib as mpl

mpl.use('Agg')

import matplotlib.pyplot as plt


# G_pc3_yr_2_msun_1 = 4.5e-15

# runs = [ ["binary_ecc0_norm_fast","binary_ecc0"]
#         ,["binary_ecc0_HPM_MM_radoff","binary_ecc0"]
#         ,["binary_ecc0_HPM_radoff_long","binary_ecc0"] ]

# gizmoDir = "/export/2/lb1g19/data"

# runs = [ ["classic","test_two_BH"]]
# 
# gizmoDir = "/srv/djw1g16/gizmos/bh_merge/"

# runs = [ ["devel_norad","open_binary_test"]]
# runs = [ ["devel_norad","open_binary_inslice_cont"]]
# runs = [ ["rad_cont","open_binary_rad_cont_inslice_evolved"]]
# runs = [ ["devel_independent","test_0m0001"]]
# runs = [ ["devel_norad","open_binary_test"], ["devel_norad","open_binary_inslice_cont"],
#           ["rad_cont","open_binary_rad_cont_inslice_full"]]
#runs = [["devel_norad", "scaled_binary_stable"]]
runs = [["norad_tracktorque", "scaled_binary_stable_test"]]

gizmoDir = "/srv/djw1g16/gizmos/"


def load_gadget(run_id, output_dir, snap_str):
    return gizmo_tools.load_gizmo_nbody(run_id, output_dir, snap_str=snap_str, gizmoDir=gizmoDir
                                        , load_binary_headers=True, only_header=True)


def extract_bh_data(run_id, output_dir, gizmodir=gizmoDir, snap0=0, maxsnapf=-1, nprocs=64):
    snapi = 0
    if snap0 > 0:
        snapi = snap0
    snapf = gizmo_tools.lastConsecutiveSnapshot(run_id, output_dir, False, gizmoDir=gizmodir)

    # max run
    if -1 < maxsnapf < snapf:
        print("Forcing snapf from {} to {}".format(snapf, maxsnapf))
        snapf = maxsnapf

    print("nfiles:", snapf - snapi + 1)

    isnaps = range(snapi, snapf + 1)

    snap_strs = [f"{i:03d}" for i in isnaps]
    extract_header_func = functools.partial(load_gadget, run_id, output_dir)

    with Pool(processes=nprocs) as pool:
        bh_binary_alltimes = pool.map(extract_header_func, snap_strs)

    # convert list of dicts of 1D arrays to dict of 2D arrays
    d = {k: np.array([dic[k] for dic in bh_binary_alltimes]) for k in bh_binary_alltimes[0]}

    # load force manually
    bh_force = np.loadtxt(os.path.join(this_dir, f"data/bh_forces_{run_id}_{output_dir}.dat"))
    # Units: km**2 kpc**-1 s**-2, i.e. actually acceleration
    # Target unit: pc yr**-2, then multiply by mass in Msun to get Msun pc yr**-2
    bh_force *= pynbody.units.Unit("km**2 kpc**-1 s**-2").in_units("pc yr**-2")

    d['BH_force_1'] = bh_force[isnaps, 0:3] * d['BH_mass_1'][:, np.newaxis]
    d['BH_force_2'] = bh_force[isnaps, 3:6] * d['BH_mass_2'][:, np.newaxis]

    return d


def calc_torque_etc(d):
    # fix BH vels
    dtime = np.gradient(d['time'])
    # assume fixed dt
    # dt = dtime[0]
    d['BH_vel_1'] = np.gradient(d['BH_pos_1'], axis=0) / dtime[:, np.newaxis] / 1.e6
    d['BH_vel_2'] = np.gradient(d['BH_pos_2'], axis=0) / dtime[:, np.newaxis] / 1.e6

    d['dist'] = linalg.norm(d['BH_pos_1'] - d['BH_pos_2'], axis=1)
    d['angmom'] = d['BH_mass_1'][:, np.newaxis] * np.cross(d['BH_pos_1'], d['BH_vel_1']) \
                  + d['BH_mass_2'][:, np.newaxis] * np.cross(d['BH_pos_2'], d['BH_vel_2'])

    d['torque'] = np.cross(d['BH_pos_1'], d['BH_force_1']) + np.cross(d['BH_pos_2'], d['BH_force_2'])
    d['cum_torque'] = np.cumsum(d['torque'] * dtime[:,np.newaxis], axis=0)

    #     d['radial_force_1'] = (d['BH_pos_1']-d['BH_pos_2'])*d['BH_force_1']/linalg.norm(d['BH_pos_1']-d['BH_pos_2'],axis=1)[:,np.newaxis]
    #     d['radial_force_2'] = (d['BH_pos_2']-d['BH_pos_1'])*d['BH_force_2']/linalg.norm(d['BH_pos_2']-d['BH_pos_1'],axis=1)[:,np.newaxis]
    d['radial_force_1'] = d['BH_pos_1'] * d['BH_force_1'] / linalg.norm(d['BH_pos_1'], axis=1)[:, np.newaxis]
    d['radial_force_2'] = d['BH_pos_2'] * d['BH_force_2'] / linalg.norm(d['BH_pos_2'], axis=1)[:, np.newaxis]

    #     d['radial_acc_mom_1'] = d['BH_pos_1']*d['BH_acc_mom_1']/linalg.norm(d['BH_pos_1'],axis=1)[:,np.newaxis]
    #     d['radial_acc_mom_2'] = d['BH_pos_2']*d['BH_acc_mom_2']/linalg.norm(d['BH_pos_2'],axis=1)[:,np.newaxis]

    #     d['acc_cum_torque'] = np.cross(d['BH_pos_1'],d['BH_acc_mom_1']) + np.cross(d['BH_pos_2'],d['BH_acc_mom_2'])
    d['acc_force_1'] = np.gradient(d['BH_acc_mom_1'], axis=0) / dtime[:,np.newaxis]
    d['acc_force_2'] = np.gradient(d['BH_acc_mom_2'], axis=0) / dtime[:,np.newaxis]
    d['acc_torque'] = np.cross(d['BH_pos_1'], d['acc_force_1']) + np.cross(d['BH_pos_2'], d['acc_force_2'])
    d['acc_cum_torque'] = np.cumsum(d['acc_torque'] * dtime[:,np.newaxis], axis=0)

    d['net_force'] = d['BH_force_1'] + d['BH_force_2']
    d['acc_net_force'] = d['acc_force_1'] + d['acc_force_2']

    #     d['BH_BH_accel_1'] = G_pc3_yr_2_msun_1*d['BH_mass_1'][:,np.newaxis]*d['BH_mass_2'][:,np.newaxis]*(d['BH_pos_2']-d['BH_pos_1'])/d['dist'][:,np.newaxis]**3
    #     d['BH_BH_accel_2'] =-d['BH_BH_accel_1']

    # changing sign to get positive torque, just because it looks tidier
    d['angmom'] *= -1
    d['torque'] *= -1


def next_axis(sp):
    irow = 0
    icol = 0

    while True:
        print(irow, icol)
        yield sp[irow, icol]
        icol += 1
        if icol >= sp.shape[1]:
            icol = 0
            irow += 1


def plot_torque(run_id, output_dir, d, tmax=0.05):
    time = d['time']
    dist = d['dist']
    angmom = d['angmom']
    torque = d['torque']

    plots = [
        #                'acc_force_1'
        #                 ,'grav_force_1'
        #                 ,'acc_force_2'
        #                 ,'grav_force_2'
        'gtorque'
        , 'acc_torque'
        , 'gcumtorque'
        , 'acc_cumtorque'

    ]

    nplots = len(plots)
    nrows = int(np.ceil(nplots / 2))
    ncols = 2

    #     nrows=7
    #     ncols=2

    scale = 3.
    fig, sp = plt.subplots(nrows, ncols, sharex=True, constrained_layout=True,
                           figsize=(ncols * scale * 2, nrows * scale), squeeze=False)

    print(plots, nplots)

    for iaxis, axis in enumerate(["x", "y", "z"]):
        na = next_axis(sp)

        if 'dist' in plots:
            ax = next(na)
            if iaxis == 0:
                ax.plot(time, dist)
                ax.set_ylabel(r'$d$ (pc)')

        if 'angmom' in plots:
            ax = next(na)
            ax.plot(time, angmom[:, iaxis], label=r"$L_{{{}}}$".format(axis))
            ax.set_ylabel('Angular momentum\n(M$_{\odot}$ pc$^2$ / yr)')

        if 'gtorque' in plots:
            ax = next(na)
            ax.plot(time, torque[:, iaxis], label=r"$\tau_{{{}}}$".format(axis))
            ax.set_ylabel('Grav Torque\n(M$_{\odot}$ pc$^2$ / yr$^{{-2}}$)')
        #         ax.plot(time,angmom[:,iaxis]/torque[:,iaxis],label=r"$'\tau_{{{}}}$".format(axis))
        #         ax.plot(time,d['net_force'][:,iaxis],label=r"$F_{{{}}}$".format(axis))
        if 'gcumtorque' in plots:
            ax = next(na)
            ax.plot(time, d['cum_torque'][:, iaxis], label=r"$\int\tau_{{{}}}dt$".format(axis))
            ax.set_ylabel('Integrated Grav Torque\n(M$_{\odot}$ pc$^2$ / yr)')

        if 'bh_pos_1' in plots:
            ax = next(na)
            ax.plot(time, d['BH_pos_1'][:, iaxis], label=r"${{{}}}_1$".format(axis))
            ax.set_ylabel('Pos 1')
        if 'bh_pos_2' in plots:
            ax = next(na)
            ax.plot(time, d['BH_pos_2'][:, iaxis], label=r"${{{}}}_2$".format(axis))
            ax.set_ylabel('Pos 2')

        if 'bh_vel_1' in plots:
            ax = next(na)
            ax.plot(time, d['BH_vel_1'][:, iaxis], label=r"$v_{{{},1}}$".format(axis))
            ax.set_ylabel('Vel 1')
        if 'bh_pos_2' in plots:
            ax = next(na)
            ax.plot(time, d['BH_vel_2'][:, iaxis], label=r"$v_{{{},2}}$".format(axis))
            ax.set_ylabel('Vel 2')

        if 'grav_force_1' in plots:
            ax = next(na)
            ax.plot(time, d['BH_force_1'][:, iaxis], label=r"$F_{{{},1}}$".format(axis))
            ax.set_ylabel('Grav Force 1')
        if 'grav_force_2' in plots:
            ax = next(na)
            ax.plot(time, d['BH_force_2'][:, iaxis], label=r"$F_{{{},2}}$".format(axis))
            ax.set_ylabel('Grav Force 2')

        if 'acc_force_1' in plots:
            ax = next(na)
            ax.plot(time, d['acc_force_1'][:, iaxis], label=r"$F_{{{},1}}$".format(axis))
            ax.set_ylabel('Acc Force 1')
        if 'acc_force_2' in plots:
            ax = next(na)
            ax.plot(time, d['acc_force_2'][:, iaxis], label=r"$F_{{{},2}}$".format(axis))
            ax.set_ylabel('Acc Force 2')

        if 'bh_circradforce_1' in plots:
            ax = next(na)
            ax.plot(time, (d['BH_force_1'][:, iaxis] - d['radial_force_1'][:, iaxis]),
                    label=r"$F_{{c,{},1}}$".format(axis))
            ax.plot(time, (d['radial_force_1'][:, iaxis]), label=r"$F_{{r,{},1}}$".format(axis))
            ax.set_ylabel('Radial/non-radial force 1')
        if 'bh_circradforce_2' in plots:
            ax = next(na)
            ax.plot(time, (d['BH_force_2'][:, iaxis] - d['radial_force_2'][:, iaxis]),
                    label=r"$F_{{c,{},2}}$".format(axis))
            ax.plot(time, (d['radial_force_2'][:, iaxis]), label=r"$F_{{r,{},2}}$".format(axis, axis))
            ax.set_ylabel('Radial/non-radial force 2')

        if 'bh_accmom' in plots:
            ax = next(na)
            ax.plot(time, d['BH_acc_mom_1'][:, iaxis], label=r"$P_{{a,{},1}}$".format(axis))
            ax.plot(time, d['BH_acc_mom_2'][:, iaxis], label=r"$P_{{a,{},2}}$".format(axis))
            ax.set_ylabel('Accreted momentum')

        if 'acc_cumtorque' in plots:
            ax = next(na)
            ax.plot(time, d['acc_torque'][:, iaxis], label=r"$\tau_{{{}}}$".format(axis))
            ax.set_ylabel('Accreted torque')

        if 'acc_cumtorque' in plots:
            ax = next(na)
            ax.plot(time, d['acc_cum_torque'][:, iaxis], label=r"$\int\tau_{{{}}}dt$".format(axis))
            ax.set_ylabel('Integrated Accretion Torque')

    #     sp[6,0].plot(time,(1.-d['BH_50_1']),label='1')
    #     sp[6,0].plot(time,(1.-d['BH_50_2']),label='2')

    #         sp[5,0].plot(time,d['BH_force_1'][:,iaxis]-d['BH_BH_accel_1'][:,iaxis],label=r"$\Delta F_{{{},1}}$".format(axis))
    #         sp[5,1].plot(time,d['BH_force_2'][:,iaxis]-d['BH_BH_accel_2'][:,iaxis],label=r"$\Delta F_{{{},2}}$".format(axis))
    #     sp[1,1].set_ylabel('Torque time-scale\n(yr)')
    #     sp[1,1].set_ylabel('Force sum')

    #     sp[6,0].set_ylabel('BH force 50 centile')
    #     sp[6,0].set_yscale('log')

    for ix in range(ncols):
        sp[-1, ix].set_xlabel(r'Time / Myr')
        #         for iy in [4,5]:
        #             sp[iy,ix].set_yscale('symlog',linthresh=1.e-1)
        #             sp[iy,ix].set_ylim(-0.06,0.06)
        for iy in range(nrows):
            sp[iy, ix].legend()
            sp[iy, ix].set_xlim([0, tmax])

    fig.savefig(os.path.join(this_dir,f"../figures/torque_{run_id}_{output_dir}.pdf"))
    plt.close('all')


def extract_plot_torque(run_id, output_dir, **kwargs):
    d = extract_bh_data(run_id, output_dir, **kwargs)
    calc_torque_etc(d)
    plot_torque(run_id, output_dir, d)


def plot_all_torques():
    for run in runs:
        extract_plot_torque(*run)


if __name__ == '__main__':
    plot_all_torques()
#     extract_plot_torque(*runs[0],patch_accretion=True)    
#     extract_plot_torque(*runs[0])    
#     extract_plot_torque(*runs[1])#,maxsnapf=20)
