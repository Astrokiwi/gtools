print("Importing")

# Import libraries to do our magic
import numpy as np
import sys
import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as plt

import gizmo_tools

import pandas as pd

from visualisation.sph_plotter import sph_plotter

# def vec_kernel(r,h):
#     x = r/h
#     k = np.zeros(x.size)
#     k[x<.5] = 1.+6.*(x[x<.5]-1.)*x[x<.5]**2
#     k[(x<1.)&(x>=.5)] = 2.*(1.-x[(x<1.)&(x>=.5)])**3
# 
#     # normalise
#     k*=8./np.pi/h**3
#     return k

if __name__=='__main__':
    run_id = "3032"
    run_name = "newflow_vesc_thin_45"
    snap_str="100"

    print("Loading")
    t,snap=gizmo_tools.load_gizmo_nbody(run_id,run_name,snap_str)
    
#     fig,sp = plt.subplots()
#     ax = sp
#     snap["mask"] = True
#     sph_plotter.set_parallel()
#     map = sph_plotter.sph_dense(snap["x"],
#                                 snap["z"],
#                                 snap["mass"],
#                                 snap["smooth"],
#                                 512,
#                                 [-0.01,-0.01],
#                                 0.02,
#                                 snap["mask"]
#                                 )
# 
#     ax.imshow(np.log10(map.T))
#     fig.savefig("pics/sph_test.png",dpi=400)
#     plt.close('all')

    print("Reducing")

    Nclouds = 10000
    Mcloud = snap["mass"].sum()/Nclouds
    
    Nparticles = len(snap)
    
    clouds = pd.DataFrame(  index=range(Nclouds),
                            columns=["x","y","z","vx","vy","vz","mass","rad"])

    clouds["mass"] = Mcloud

    cloud_ids = np.random.randint(0,len(snap),Nclouds)
    for coord in "x","y","z","vx","vy","vz":
        clouds[coord] = snap[coord][cloud_ids]
    clouds["rad"] = (Mcloud/snap["rho"][cloud_ids]*3./4./np.pi)**(1,3)
    
    for key in "x","y","z","rad":
        clouds[key]*=1.e3 # kpc to pc
    clouds["mass"]*=1.e10 # to Msol
    #n.b. v already in km/s
    
    clouds.to_csv(f"data/clouds_{run_name}_{snap_str}.txt",sep=' ',index=False)
    
    print("Plotting")
    fig,sp = plt.subplots()
    ax = sp
    lower = np.min([np.min(clouds["x"]),
                    np.min(clouds["z"])])
    upper = np.min([np.max(clouds["x"]),
                    np.max(clouds["z"])])
    ax.set_xlim(lower,upper)
    ax.set_ylim(lower,upper)
    for index,cloud in clouds.iterrows():
        ax.add_artist(
                        plt.Circle(
                                    (cloud["x"],cloud["z"]),
                                    cloud["rad"]
                                )
                    )
    
    fig.savefig("pics/cloud_test.png",dpi=400)
    plt.close('all')
    
    
    
    
    
    
    
    
        
