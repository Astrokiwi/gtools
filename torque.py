# Import libraries to do our magic
import tblib.pickling_support
tblib.pickling_support.install()

import math

import numpy as np
import h5py
import matplotlib as mpl
#mpl.use('Agg')

from matplotlib import colors,rc
import matplotlib.pyplot as plt
rc('text', usetex=True)
rc('font',family='serif')

from sys import path, version, exit
path.append("src/")
import tab_interp

import sys

import itertools

import gizmo_tools
import argparse
from multiprocessing import Pool

def load_gadget(i):
    f = h5py.File(infiles[i],"r")
    
    header = f["/Header"]
    BH_data= f["/BH_binary"]
#    print(BH_data.attrs.get("BH_pos_1"))
#    print(list(header.attrs.keys()))
#    print(header)

    BH_binary = dict()
    BH_binary['time'] = header.attrs.get("Time")
    BH_binary['time']*= 0.9778e9 # to yr
    BH_binary['time']/=1.e6 # to Myr

    BH_binary['BH_pos_1'] = np.array(BH_data.attrs.get("Binary_pos_1"))
    BH_binary['BH_pos_2'] = np.array(BH_data.attrs.get("Binary_pos_2"))
    BH_binary['BH_vel_1'] = np.array(BH_data.attrs.get("Binary_vel_1"))
    BH_binary['BH_vel_2'] = np.array(BH_data.attrs.get("Binary_vel_2"))
    BH_binary['BH_force_1'] = np.array(BH_data.attrs.get("Binary_force_1"))
    BH_binary['BH_force_2'] = np.array(BH_data.attrs.get("Binary_force_2"))
    BH_binary['BH_mass_1'] = BH_data.attrs.get("Binary_mass_1")
    BH_binary['BH_mass_2'] = BH_data.attrs.get("Binary_mass_2")

    BH_binary['BH_pos_1'] *= 1000.0           # pc
    BH_binary['BH_pos_2'] *= 1000.0           # pc
    BH_binary['BH_vel_1'] *= 1000.0/0.9778e9  # pc / yr
    BH_binary['BH_vel_2'] *= 1000.0/0.9778e9  # pc / yr
    BH_binary['BH_force_1'] *= 1e10 * 1000.0/(0.9778e9)**2  # msun * pc / yr**2
    BH_binary['BH_force_2'] *= 1e10 * 1000.0/(0.9778e9)**2  # msun * pc / yr**2
    BH_binary['BH_mass_1'] *= 1e10            # msun
    BH_binary['BH_mass_2'] *= 1e10            # msun

    return BH_binary


if __name__ == '__main__':
    default_values = dict()
    default_values["nprocs"]=8
    default_values["maxsnapf"]=-1
    default_values["snap0"]=-1
    default_values["L"]=400
    default_values["data_ranges"]=""
 
    parsevals = ["data_ranges","nprocs","maxsnapf","run_id","output_dir","snap0"]

    parser = argparse.ArgumentParser()
    parser.add_argument('run_id',help="name of superdirectory for runs")
    parser.add_argument('output_dir',help="name of subdirectory for run")
    parser.add_argument('--nprocs',type=int,help="processors to run on (default {})".format(default_values["nprocs"]))
    parser.add_argument('--maxsnapf',type=int,help="snapshot to end on (default=-1=do all snapshots)")
    parser.add_argument('--snap0',type=int,help="snapshot to start on (default=-1=do all snapshots)")
    parser.add_argument('--L',type=int,help="size of plot area in pixels")
    parser.add_argument('--data_ranges',type=str,help="explicit data bounds for plot - format: --data_ranges=min1,max1+min2,max2+min3,max3 etc; the equals sign may be necessary!")
    args = parser.parse_args()
    
    # TODO: don't do this!
    for parseval in parsevals:
        if ( vars(args)[parseval] ):
            vars()[parseval] = vars(args)[parseval]
            print("setting {} to {}, current value:{}".format(parseval,vars(args)[parseval],vars()[parseval]))
        else:
            if ( parseval in default_values ):
                vars()[parseval] = default_values[parseval]
                print("setting {} to default {}, current value:{}".format(parseval,default_values[parseval],vars()[parseval]))
            else:
                raise Exception("No default value for {} - it must be specified!".format(parseval))


    if len(data_ranges)<=0:
        data_ranges = None
    else:
        data_ranges=[[float(y) for y in x.split(",")] for x in data_ranges.split("+")]
    
    snapi = 0
    if snap0>0:
        snapi=snap0
    snapf = gizmo_tools.lastConsecutiveSnapshot(run_id,output_dir,False)

    gizmoDir = gizmo_tools.getGizmoDir(run_id)
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

    time = np.zeros(nfiles,dtype=float)
    dist = np.zeros(nfiles,dtype=float)
    angmom = np.zeros((nfiles,3),dtype=float)
    torque = np.zeros((nfiles,3),dtype=float)
    angmom_zyl = np.zeros((nfiles,3),dtype=float)
    torque_zyl = np.zeros((nfiles,3),dtype=float)

    
    with Pool(processes=nprocs) as pool:
      BH_Binary_alltimes = pool.map(load_gadget,range(snapf+1-snapi))


    for i in range(nfiles):
      time[i] = maps[i]['time']
      dist[i] = math.sqrt((maps[i]['BH_pos_1'][0]-maps[i]['BH_pos_2'][0])**2 + (maps[i]['BH_pos_1'][1]-maps[i]['BH_pos_2'][1])**2 + (maps[i]['BH_pos_1'][2]-maps[i]['BH_pos_2'][2])**2)

      angmom[i][0] = maps[i]['BH_mass_1']*(maps[i]['BH_pos_1'][1]*maps[i]['BH_vel_1'][2]-maps[i]['BH_pos_1'][2]*maps[i]['BH_vel_1'][1]) + maps[i]['BH_mass_2']*(maps[i]['BH_pos_2'][1]*maps[i]['BH_vel_2'][2]-maps[i]['BH_pos_2'][2]*maps[i]['BH_vel_2'][1])
      angmom[i][1] = maps[i]['BH_mass_1']*(maps[i]['BH_pos_1'][2]*maps[i]['BH_vel_1'][0]-maps[i]['BH_pos_1'][0]*maps[i]['BH_vel_1'][2]) + maps[i]['BH_mass_2']*(maps[i]['BH_pos_2'][2]*maps[i]['BH_vel_2'][0]-maps[i]['BH_pos_2'][0]*maps[i]['BH_vel_2'][2])
      angmom[i][2] = maps[i]['BH_mass_1']*(maps[i]['BH_pos_1'][0]*maps[i]['BH_vel_1'][1]-maps[i]['BH_pos_1'][1]*maps[i]['BH_vel_1'][0]) + maps[i]['BH_mass_2']*(maps[i]['BH_pos_2'][0]*maps[i]['BH_vel_2'][1]-maps[i]['BH_pos_2'][1]*maps[i]['BH_vel_2'][0])


      torque[i][0] = (maps[i]['BH_pos_1'][1]*maps[i]['BH_force_1'][2]-maps[i]['BH_pos_1'][2]*maps[i]['BH_force_1'][1]) + (maps[i]['BH_pos_2'][1]*maps[i]['BH_force_2'][2]-maps[i]['BH_pos_2'][2]*maps[i]['BH_force_2'][1])
      torque[i][1] = (maps[i]['BH_pos_1'][2]*maps[i]['BH_force_1'][0]-maps[i]['BH_pos_1'][0]*maps[i]['BH_force_1'][2]) + (maps[i]['BH_pos_2'][2]*maps[i]['BH_force_2'][0]-maps[i]['BH_pos_2'][0]*maps[i]['BH_force_2'][2])
      torque[i][2] = (maps[i]['BH_pos_1'][0]*maps[i]['BH_force_1'][1]-maps[i]['BH_pos_1'][1]*maps[i]['BH_force_1'][0]) + (maps[i]['BH_pos_2'][0]*maps[i]['BH_force_2'][1]-maps[i]['BH_pos_2'][1]*maps[i]['BH_force_2'][0])


      rho1     = math.sqrt(maps[i]['BH_pos_1'][0]**2 + maps[i]['BH_pos_1'][1]**2)
      rho1_dot = (maps[i]['BH_pos_1'][0]*maps[i]['BH_vel_1'][0] + maps[i]['BH_pos_1'][1]*maps[i]['BH_vel_1'][1])/rho1
      phi1     = math.atan2(maps[i]['BH_pos_1'][1],maps[i]['BH_pos_1'][0])
      phi1_dot = (maps[i]['BH_pos_1'][0]*maps[i]['BH_vel_1'][1] - maps[i]['BH_pos_1'][1]*maps[i]['BH_vel_1'][0])/rho1**2
      z1       = maps[i]['BH_pos_1'][2]
      z1_dot   = maps[i]['BH_vel_1'][2]
      rho2     = math.sqrt(maps[i]['BH_pos_2'][0]**2 + maps[i]['BH_pos_2'][1]**2)
      rho2_dot = (maps[i]['BH_pos_2'][0]*maps[i]['BH_vel_2'][0] + maps[i]['BH_pos_2'][1]*maps[i]['BH_vel_2'][1])/rho2
      phi2     = math.atan2(maps[i]['BH_pos_2'][1],maps[i]['BH_pos_2'][0])
      phi2_dot = (maps[i]['BH_pos_2'][0]*maps[i]['BH_vel_2'][1] - maps[i]['BH_pos_2'][1]*maps[i]['BH_vel_2'][0])/rho2**2
      z2       = maps[i]['BH_pos_2'][2]
      z2_dot   = maps[i]['BH_vel_2'][2]

      angmom_zyl[0] = -maps[i]['BH_mass_1']*rho1*phi1_dot*z1 - maps[i]['BH_mass_2']*rho2*phi2_dot*z2
      angmom_zyl[1] = maps[i]['BH_mass_1']*(rho1_dot*z1-rho1*z1_dot) + maps[i]['BH_mass_2']*(rho2_dot*z2-rho2*z2_dot)
      angmom_zyl[1] = maps[i]['BH_mass_1']*rho1**2*phi1_dot + maps[i]['BH_mass_2']*rho2**2*phi2_dot

      torque_zyl[0] = -maps[i]['BH_mass_1']*rho1*phi1_dot*z1 - maps[i]['BH_mass_2']*rho2*phi2_dot*z2
      torque_zyl[1] = maps[i]['BH_mass_1']*(rho1_dot*z1-rho1*z1_dot) + maps[i]['BH_mass_2']*(rho2_dot*z2-rho2*z2_dot)
      torque_zyl[1] = maps[i]['BH_mass_1']*rho1**2*phi1_dot + maps[i]['BH_mass_2']*rho2**2*phi2_dot

      print(maps[i]['BH_pos_1'][0]*maps[i]['BH_vel_1'][1],maps[i]['BH_pos_1'][1]*maps[i]['BH_vel_1'][0],maps[i]['BH_pos_2'][0]*maps[i]['BH_vel_2'][1],maps[i]['BH_pos_2'][1]*maps[i]['BH_vel_2'][0],maps[i]['BH_pos_1'][0]*maps[i]['BH_force_1'][1],maps[i]['BH_pos_1'][1]*maps[i]['BH_force_1'][0],maps[i]['BH_pos_2'][0]*maps[i]['BH_force_2'][1],maps[i]['BH_pos_2'][1]*maps[i]['BH_force_2'][0])

ax0 = plt.subplot(221)
ax1 = plt.subplot(222)
ax12=ax1.twinx()
ax2 = plt.subplot(223)
ax22=ax2.twinx()
ax3 = plt.subplot(224)
ax32=ax3.twinx()

ax0.plot(time,dist)
ax0.set_xlabel(r'Time / Myr')
ax0.set_ylabel(r'distance / pc')

ax1.plot(time,angmom[:,0],'b')
ax1.set_xlabel(r'Time / Myr')
ax1.set_ylabel(r'$L_x / M_{\odot} pc / yr$')
ax12.plot(time,torque[:,0],'r')
ax12.set_ylabel(r'$M_x / M_{\odot} pc^2 / yr^2$')
ax12.yaxis.label.set_color('red')

ax2.plot(time,angmom[:,1],'b')
ax2.set_xlabel(r'Time / Myr')
ax2.set_ylabel(r'$L_y / M_{\odot} pc / yr$')
ax22.plot(time,torque[:,1],'r')
ax22.set_ylabel(r'$M_y / M_{\odot} pc^2 / yr^2$')
ax22.yaxis.label.set_color('red')

ax3.plot(time,angmom[:,2],'b')
ax3.set_xlabel(r'Time / Myr')
ax3.set_ylabel(r'$L_z / M_{\odot} pc / yr$')
ax32.plot(time,torque[:,2],'r')
ax32.set_ylabel(r'$M_z / M_{\odot} pc^2 / yr^2$')
ax32.yaxis.label.set_color('red')

plt.show()

