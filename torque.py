# Import libraries to do our magic
import tblib.pickling_support
tblib.pickling_support.install()

import math

import numpy as np
import h5py
import matplotlib as mpl
mpl.use('Agg')

from matplotlib import colors,rc
import matplotlib.pyplot as plt
# rc('text', usetex=True)
# rc('font',family='serif')

from sys import path, version, exit
path.append("src/")
import tab_interp

import sys

import itertools

import gizmo_tools
import argparse
from multiprocessing import Pool

import pandas as pd
from scipy import linalg
import functools


runs = [ ["binary_ecc0_norm_fast","binary_ecc0"]
        ,["binary_ecc0_HPM_MM_radoff","binary_ecc0"]
        ,["binary_ecc0_HPM_radoff_long","binary_ecc0"] ]


def load_gadget(run_id,output_dir,snap_str):
    return gizmo_tools.load_gizmo_nbody(run_id,output_dir,snap_str=snap_str,gizmoDir="/export/2/lb1g19/data"
                                        ,load_binary_headers=True,only_header=True)

def extract_bh_data(run_id,output_dir,gizmoDir="/export/2/lb1g19/data/",snap0=0,maxsnapf=-1,nprocs=64):
    snapi = 0
    if snap0>0:
        snapi=snap0
    snapf = gizmo_tools.lastConsecutiveSnapshot(run_id,output_dir,False,gizmoDir=gizmoDir)

#     gizmoDir = gizmo_tools.getGizmoDir(run_id)
    gizmoDir = "/export/2/lb1g19/data"
    movieDir = gizmo_tools.getMovieDir()
    fullDir = gizmoDir+"/"+run_id+"/"+output_dir
    
# max run
#     if ( maxsnapf>-1 and snapf>maxsnapf ):
    if maxsnapf>-1:
        print("Forcing snapf from {} to {}".format(snapf,maxsnapf))
        snapf = maxsnapf

    print("nfiles:",snapf-snapi+1)
    
    isnaps = range(snapi,snapf+1)
    
    infiles = [fullDir+"/snapshot_"+("%03d" % snapx)+".hdf5" for snapx in isnaps]
    
    nfiles = len(infiles)

    angmom = np.zeros((nfiles,3),dtype=float)
    torque = np.zeros((nfiles,3),dtype=float)
    angmom_zyl = np.zeros((nfiles,3),dtype=float)
    torque_zyl = np.zeros((nfiles,3),dtype=float)

    snap_strs = [f"{i:03d}" for i in isnaps]
    extract_header_func = functools.partial(load_gadget,run_id,output_dir)
    
    with Pool(processes=nprocs) as pool:
#       BH_Binary_alltimes = pool.map(load_gadget,range(snapf+1-snapi))
      BH_Binary_alltimes = pool.map(extract_header_func,snap_strs)

    d={k: np.array([dic[k] for dic in BH_Binary_alltimes]) for k in BH_Binary_alltimes[0]} 
    return d

def calc_torque_etc(d):
    # fix BH vels
    dtime = np.gradient(d['time'])
    d['BH_vel_1']=np.gradient(d['BH_pos_1'],axis=0)/dtime[:,np.newaxis]/1.e6
    d['BH_vel_2']=np.gradient(d['BH_pos_2'],axis=0)/dtime[:,np.newaxis]/1.e6

    d['dist'] = linalg.norm(d['BH_pos_1']-d['BH_pos_2'],axis=1)
    d['angmom'] = d['BH_mass_1'][:,np.newaxis]*np.cross(d['BH_pos_1'],d['BH_vel_1'])  \
             + d['BH_mass_2'][:,np.newaxis]*np.cross(d['BH_pos_2'],d['BH_vel_2'])
    d['torque'] = d['BH_mass_1'][:,np.newaxis]*np.cross(d['BH_pos_1'],d['BH_force_1'])  \
             + d['BH_mass_2'][:,np.newaxis]*np.cross(d['BH_pos_2'],d['BH_force_2'])

    # changing sign to get positive torque, just because it looks tidier
    d['angmom']*=-1
    d['torque']*=-1

def plot_torque(run_id,output_dir,d):
    time = d['time']
    dist = d['dist']
    angmom = d['angmom']
    torque = d['torque']

    nrows=5
    ncols=2

    scale = 2.
    fig,sp = plt.subplots(nrows,ncols,sharex=True,constrained_layout=True,figsize=(ncols*scale*2,nrows*scale))

    sp[0,0].plot(time,dist)
    sp[0,0].set_ylabel(r'$d$ (pc)')
    
    time_offset = (time[:-1]+time[1:])/2
    dtime = np.gradient(time)

    for iaxis,axis in enumerate(["x","y","z"]):
        sp[0,1].plot(time,angmom[:,iaxis],label=r"$L_{{{}}}$".format(axis))

        sp[1,0].plot(time,torque[:,iaxis],label=r"$\tau_{{{}}}$".format(axis))
        sp[1,1].plot(time,angmom[:,iaxis]/torque[:,iaxis],label=r"$t_{{{}}}$".format(axis))

        sp[2,0].plot(time,d['BH_pos_1'][:,iaxis],label=r"${{{}}}_1$".format(axis))
        sp[2,1].plot(time,d['BH_pos_2'][:,iaxis],label=r"${{{}}}_2$".format(axis))

        sp[3,0].plot(time,d['BH_vel_1'][:,iaxis],label=r"$v_{{{},1}}$".format(axis))
        sp[3,1].plot(time,d['BH_vel_2'][:,iaxis],label=r"$v_{{{},2}}$".format(axis))

        sp[4,0].plot(time,d['BH_force_1'][:,iaxis],label=r"$F_{{{},1}}$".format(axis))
        sp[4,1].plot(time,d['BH_force_2'][:,iaxis],label=r"$F_{{{},2}}$".format(axis))
    sp[0,1].set_ylabel('Angular momentum\n(M$_{\odot}$ pc$^2$ / yr)')
    sp[1,0].set_ylabel('Torque\n(M$_{\odot}$ pc$^2$ / yr$^{{-2}}$)')
    sp[1,1].set_ylabel('Torque time-scale\n(yr)')
    sp[2,0].set_ylabel('Pos 1')
    sp[2,1].set_ylabel('Pos 2')
    sp[3,0].set_ylabel('Vel 1')
    sp[3,1].set_ylabel('Vel 2')
    sp[4,0].set_ylabel('Force 1')
    sp[4,1].set_ylabel('Force 2')
    
    for ix in range(ncols):
        sp[-1,ix].set_xlabel(r'Time / Myr')
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
#     plot_all_torques()      
    extract_plot_torque(*runs[0])      
