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
import itertools

import pynbody

G_pc3_yr_2_msun_1 = 4.5e-15

runs = [ ["binary_ecc0_norm_fast","binary_ecc0"]
        ,["binary_ecc0_HPM_MM_radoff","binary_ecc0"]
        ,["binary_ecc0_HPM_radoff_long","binary_ecc0"] ]



def load_gadget(run_id,output_dir,snap_str,gizmoDir=None):
    print(snap_str)
    return gizmo_tools.load_gizmo_nbody(run_id,output_dir,snap_str=snap_str,gizmoDir=gizmoDir
                                        ,load_binary_headers=True)

def nbody_norm(pos):
    return (pos**2).sum(axis=1)**(3,2)


def bh_force(snap,bh_pos,bh_mas):
    dr = snap['pos']-bh_pos
    all_forces = dr*snap['mass'][:,np.newaxis]/nbody_norm(dr)[:,np.newaxis]
    force = pynbody.units.G*np.sum(all_forces,axis=0)
    force = force.in_original_units() # km**2 kpc**-1 s**-2
    
    force_mag = np.cumsum(np.sort(linalg.norm(all_forces,axis=1)))
    force_mag/=force_mag[-1]
    half_force_loc = np.searchsorted(force_mag,.5)/force_mag.size
    
    
    print("mass:",snap['mass'].in_units("Msol").sum())

    return force,half_force_loc

def bh_force_mag(snap,bh_pos,bh_mas):
    dr = snap['pos']-bh_pos
    all_forces = dr*snap['mass'][:,np.newaxis]/nbody_norm(dr)[:,np.newaxis]
    force_mag = np.sort(linalg.norm(all_forces,axis=1))
    return force_mag

def bh_force_mag(force):
    dr = snap['pos']-bh_pos
    all_forces = dr*snap['mass'][:,np.newaxis]/nbody_norm(dr)[:,np.newaxis]
    force_mag = np.sort(linalg.norm(all_forces,axis=1))
    return force_mag

def bh_forces(header,snap):
    bh_positions = [header['BH_pos_1'],header['BH_pos_2']]
    bh_masses = [header['BH_mass_1'],header['BH_mass_2']]
#     output_grid = [bh_force(snap,bh_pos,bh_mass) for bh_pos,bh_mass in zip(bh_positions,bh_masses)]

    return [bh_force(snap,bh_pos,bh_mass) for bh_pos,bh_mass in zip(bh_positions,bh_masses)]

def bh_force_load_calc(*args,**kwargs):
    return bh_forces(
                        *load_gadget(*args,**kwargs)
                    )

def calc_smbh_force(run_id,output_dir,gizmoDir="/export/2/lb1g19/data/",snap0=0,maxsnapf=-1,nprocs=64):

    snapi = 0
    if snap0>0:
        snapi=snap0
    snapf = gizmo_tools.lastConsecutiveSnapshot(run_id,output_dir,False,gizmoDir=gizmoDir)

    if ( maxsnapf>-1 and snapf>maxsnapf ):
        print("Forcing snapf from {} to {}".format(snapf,maxsnapf))
        snapf = maxsnapf

    snap_strs = [f"{i:03d}" for i in range(snapi,snapf+1)]
    
    
    
    with Pool(nprocs) as pool:
        snap_forces = pool.starmap(bh_force_load_calc,zip(
                                                            itertools.repeat(run_id)
                                                            ,itertools.repeat(output_dir)
                                                            ,snap_strs
                                                            ,itertools.repeat(gizmoDir)
                                                            ))
    just_forces = [[x[0] for x in y] for y in snap_forces]
    np.savetxt(f"data/bh_forces_{run_id}_{output_dir}.dat",np.array(just_forces).reshape(len(just_forces),6)) 

    force_50centiles = [[x[1] for x in y] for y in snap_forces]
    np.savetxt(f"data/bh_50cent_{run_id}_{output_dir}.dat",np.array(force_50centiles).reshape(len(force_50centiles),2)) 

    return snap_forces
#     for snap_str in ["050"]:
#         forces = bh_force_load_calc(run_id,output_dir,snap_str,gizmoDir)
#         print(forces)

def calc_smbh_forces(runs):
    if type(runs[0]) is list:
        output = []
        for run in runs:
            output+=[calc_smbh_force(*run,maxsnapf=1)]
        return output
    else:
        return calc_smbh_force(*runs)

if __name__ == '__main__':
     snap_forces=calc_smbh_forces(runs)
