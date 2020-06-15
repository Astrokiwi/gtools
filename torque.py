# Import libraries to do our magic

import math

import numpy as np
import matplotlib as mpl
mpl.use('Agg')

from matplotlib import colors,rc
import matplotlib.pyplot as plt

import gizmo_tools
from multiprocessing import Pool

import pandas as pd
from scipy import linalg
import functools

import pynbody
import sys

# G_pc3_yr_2_msun_1 = 4.5e-15

# runs = [ ["binary_ecc0_norm_fast","binary_ecc0"]
#         ,["binary_ecc0_HPM_MM_radoff","binary_ecc0"]
#         ,["binary_ecc0_HPM_radoff_long","binary_ecc0"] ]

# gizmoDir = "/export/2/lb1g19/data"

runs = [ ["classic","test_two_BH"]]

gizmoDir = "/srv/djw1g16/gizmos/bh_merge/"


def load_gadget(run_id,output_dir,snap_str):
    return gizmo_tools.load_gizmo_nbody(run_id,output_dir,snap_str=snap_str,gizmoDir=gizmoDir
                                        ,load_binary_headers=True,only_header=True)

def extract_bh_data(run_id,output_dir,gizmoDir=gizmoDir,snap0=0,maxsnapf=-1,nprocs=64):
    snapi = 0
    if snap0>0:
        snapi=snap0
    snapf = gizmo_tools.lastConsecutiveSnapshot(run_id,output_dir,False,gizmoDir=gizmoDir)
    
# max run
    if ( maxsnapf>-1 and snapf>maxsnapf ):
        print("Forcing snapf from {} to {}".format(snapf,maxsnapf))
        snapf = maxsnapf

    print("nfiles:",snapf-snapi+1)
    
    isnaps = range(snapi,snapf+1)
    
    nfiles = len(isnaps)

    angmom = np.zeros((nfiles,3),dtype=float)
    torque = np.zeros((nfiles,3),dtype=float)
    angmom_zyl = np.zeros((nfiles,3),dtype=float)
    torque_zyl = np.zeros((nfiles,3),dtype=float)

    snap_strs = [f"{i:03d}" for i in isnaps]
    extract_header_func = functools.partial(load_gadget,run_id,output_dir)
    
    with Pool(processes=nprocs) as pool:
#       BH_Binary_alltimes = pool.map(load_gadget,range(snapf+1-snapi))
      BH_Binary_alltimes = pool.map(extract_header_func,snap_strs)

    
    # convert list of dicts of 1D arrays to dict of 2D arrays
    d={k: np.array([dic[k] for dic in BH_Binary_alltimes]) for k in BH_Binary_alltimes[0]} 
    
    # load force manually
    BH_force = np.loadtxt(f"data/bh_forces_{run_id}_{output_dir}.dat")
    # Units: km**2 kpc**-1 s**-2, i.e. actually acceleration
    # Target unit: pc yr**-2, then multiply by mass in Msun to get Msun pc yr**-2
    BH_force *= pynbody.units.Unit("km**2 kpc**-1 s**-2").in_units("pc yr**-2")

    d['BH_force_1']=BH_force[:,0:3]*d['BH_mass_1'][:,np.newaxis]
    d['BH_force_2']=BH_force[:,3:6]*d['BH_mass_2'][:,np.newaxis]

    BH_50centiles = np.loadtxt(f"data/bh_50cent_{run_id}_{output_dir}.dat")
    d['BH_50_1'] = BH_50centiles[:,0]
    d['BH_50_2'] = BH_50centiles[:,1]

    return d

def calc_torque_etc(d):
    # fix BH vels
    dtime = np.gradient(d['time'])
    # assume fixed dt
    dt = dtime[0]
    d['BH_vel_1']=np.gradient(d['BH_pos_1'],axis=0)/dtime[:,np.newaxis]/1.e6
    d['BH_vel_2']=np.gradient(d['BH_pos_2'],axis=0)/dtime[:,np.newaxis]/1.e6

    d['dist'] = linalg.norm(d['BH_pos_1']-d['BH_pos_2'],axis=1)
    d['angmom'] = d['BH_mass_1'][:,np.newaxis]*np.cross(d['BH_pos_1'],d['BH_vel_1'])  \
             + d['BH_mass_2'][:,np.newaxis]*np.cross(d['BH_pos_2'],d['BH_vel_2'])
    d['torque'] = np.cross(d['BH_pos_1'],d['BH_force_1']) + np.cross(d['BH_pos_2'],d['BH_force_2'])
    d['cum_torque'] = np.cumsum(d['torque']*dt,axis=0)
    
#     d['radial_force_1'] = (d['BH_pos_1']-d['BH_pos_2'])*d['BH_force_1']/linalg.norm(d['BH_pos_1']-d['BH_pos_2'],axis=1)[:,np.newaxis]
#     d['radial_force_2'] = (d['BH_pos_2']-d['BH_pos_1'])*d['BH_force_2']/linalg.norm(d['BH_pos_2']-d['BH_pos_1'],axis=1)[:,np.newaxis]
    d['radial_force_1'] = d['BH_pos_1']*d['BH_force_1']/linalg.norm(d['BH_pos_1'],axis=1)[:,np.newaxis]
    d['radial_force_2'] = d['BH_pos_1']*d['BH_force_2']/linalg.norm(d['BH_pos_2'],axis=1)[:,np.newaxis]
    
    d['net_force'] = d['BH_force_1']+d['BH_force_2']
        
#     d['BH_BH_accel_1'] = G_pc3_yr_2_msun_1*d['BH_mass_1'][:,np.newaxis]*d['BH_mass_2'][:,np.newaxis]*(d['BH_pos_2']-d['BH_pos_1'])/d['dist'][:,np.newaxis]**3
#     d['BH_BH_accel_2'] =-d['BH_BH_accel_1']

    # changing sign to get positive torque, just because it looks tidier
    d['angmom']*=-1
    d['torque']*=-1

def plot_torque(run_id,output_dir,d):
    time = d['time']
    dist = d['dist']
    angmom = d['angmom']
    torque = d['torque']

    nrows=7
    ncols=2

    scale = 3.
    fig,sp = plt.subplots(nrows,ncols,sharex=True,constrained_layout=True,figsize=(ncols*scale*2,nrows*scale))

    sp[0,0].plot(time,dist)
    sp[0,0].set_ylabel(r'$d$ (pc)')
    
    time_offset = (time[:-1]+time[1:])/2
    dtime = np.gradient(time)

    for iaxis,axis in enumerate(["x","y","z"]):
        sp[0,1].plot(time,angmom[:,iaxis],label=r"$L_{{{}}}$".format(axis))

        sp[1,0].plot(time,torque[:,iaxis],label=r"$\tau_{{{}}}$".format(axis))
#         sp[1,1].plot(time,angmom[:,iaxis]/torque[:,iaxis],label=r"$'\tau_{{{}}}$".format(axis))
#         sp[1,1].plot(time,d['net_force'][:,iaxis],label=r"$F_{{{}}}$".format(axis))
        sp[1,1].plot(time,d['cum_torque'][:,iaxis],label=r"$\int\tau_{{{}}}dt$".format(axis))

        sp[2,0].plot(time,d['BH_pos_1'][:,iaxis],label=r"${{{}}}_1$".format(axis))
        sp[2,1].plot(time,d['BH_pos_2'][:,iaxis],label=r"${{{}}}_2$".format(axis))

        sp[3,0].plot(time,d['BH_vel_1'][:,iaxis],label=r"$v_{{{},1}}$".format(axis))
        sp[3,1].plot(time,d['BH_vel_2'][:,iaxis],label=r"$v_{{{},2}}$".format(axis))

        sp[4,0].plot(time,d['BH_force_1'][:,iaxis],label=r"$F_{{{},1}}$".format(axis))
        sp[4,1].plot(time,d['BH_force_2'][:,iaxis],label=r"$F_{{{},2}}$".format(axis))

        sp[5,0].plot(time,(d['BH_force_1'][:,iaxis]-d['radial_force_1'][:,iaxis]),label=r"$F_{{c,{},1}}$".format(axis))
        sp[5,0].plot(time,(d['radial_force_1'][:,iaxis]),label=r"$F_{{r,{},1}}$".format(axis))
        sp[5,1].plot(time,(d['BH_force_2'][:,iaxis]-d['radial_force_2'][:,iaxis]),label=r"$F_{{c,{},2}}$".format(axis))
        sp[5,1].plot(time,(d['radial_force_2'][:,iaxis]),label=r"$F_{{r,{},2}}$".format(axis,axis))
        
    sp[6,0].plot(time,(1.-d['BH_50_1']),label='1')
    sp[6,0].plot(time,(1.-d['BH_50_2']),label='2')

#         sp[5,0].plot(time,d['BH_force_1'][:,iaxis]-d['BH_BH_accel_1'][:,iaxis],label=r"$\Delta F_{{{},1}}$".format(axis))
#         sp[5,1].plot(time,d['BH_force_2'][:,iaxis]-d['BH_BH_accel_2'][:,iaxis],label=r"$\Delta F_{{{},2}}$".format(axis))
    sp[0,1].set_ylabel('Angular momentum\n(M$_{\odot}$ pc$^2$ / yr)')
    sp[1,0].set_ylabel('Torque\n(M$_{\odot}$ pc$^2$ / yr$^{{-2}}$)')
#     sp[1,1].set_ylabel('Torque time-scale\n(yr)')
#     sp[1,1].set_ylabel('Force sum')
    sp[1,1].set_ylabel('Integraed Torque\n(M$_{\odot}$ pc$^2$ / yr)')
    sp[2,0].set_ylabel('Pos 1')
    sp[2,1].set_ylabel('Pos 2')
    sp[3,0].set_ylabel('Vel 1')
    sp[3,1].set_ylabel('Vel 2')

    sp[4,0].set_ylabel('Force 1')
    sp[4,1].set_ylabel('Force 2')
    sp[5,0].set_ylabel('Radial/non-radial force 1')
    sp[5,1].set_ylabel('Radial/non-radial force 2')

    sp[6,0].set_ylabel('BH force 50 centile')
    sp[6,0].set_yscale('log')
    
    for ix in range(ncols):
        sp[-1,ix].set_xlabel(r'Time / Myr')
#         for iy in [4,5]:
#             sp[iy,ix].set_yscale('symlog',linthresh=1.e-1)
#             sp[iy,ix].set_ylim(-0.06,0.06)
        for iy in range(nrows):
            sp[iy,ix].legend()


    fig.savefig(f"../figures/torque_{run_id}_{output_dir}.pdf")
    plt.close('all')
    

def extract_plot_torque(run_id,output_dir,**kwargs):
    d = extract_bh_data(run_id,output_dir,**kwargs)
    calc_torque_etc(d)
    plot_torque(run_id,output_dir,d)


def plot_all_torques():
    for run in runs:
        extract_plot_torque(*run)

if __name__ == '__main__':
    plot_all_torques()      
#     extract_plot_torque(*runs[0])      
