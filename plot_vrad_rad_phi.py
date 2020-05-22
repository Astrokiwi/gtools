import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import parameter_mangler
import h5py

import gizmo_tools

run_id = "2030"
snap_str = "010"

plot_r = [0,6]

mode = "v"

mdots = [0.0378,0.378]
v_inits = [10,100,500]

mdot_labels = [r"$\dot{{M}}={{{}}}$ M$_\odot$ yr$^{{-1}}$".format(mdot) for mdot in mdots]
v_labels = [r"$v={{{}}}$ km s$^{{-1}}$".format(v) for v in v_inits]

# scale = .7
scale = .6

if __name__=='__main__':
    if mode!="v" and mode!="p":
        raise Exception("mode must be v (velocity) or p (momentum)")
    
    file_name = "data/{}rad_rad_phi_{}.hdf5".format(mode,run_id)
    df = pd.read_hdf(file_name,"vrad_vrad_phi")
    with h5py.File(file_name,'r') as f:
        valid_rundirs = np.array(f["valid_rundirs"])                                                                                                                                             

    manglia = parameter_mangler.load_manglia("prams/pram_options_testflows.dat")
    key_combinations = parameter_mangler.mangle_combinations(manglia)
    run_dirs = [parameter_mangler.filebase_from_combination("testflows",x) for x in key_combinations]
#     run_dirs = np.array(run_dirs,dtype=bytes)
    key_bases = [list(x.mangle_texts.keys()) for x in manglia]
    key_bases[2][1]="thinaniso" # hack for output
    mangle_n = [len(x) for x in key_bases]
    n_mangles = len(mangle_n)

    iaxis = 1
    nx = mangle_n[iaxis]
    jaxis = 0
    ny = mangle_n[jaxis]
    ny = 2 # cut out "heavy" outflow, too slow
    color_axis = 2
    colors = gizmo_tools.colors

    r_bin_edges = np.linspace(.5,7.5,8)
    n_r_bins = 7

#     ytick_values = [-50,-20,0,20,50,70,100,200,500,700]    
    ytick_values = [-50,-20,0,20,50,100,200,500]    

    nfig = len(plot_r)
    
    fig = np.empty((nfig),dtype=object)
    sp = np.empty((nfig,ny,nx),dtype=object)
    
    
    for ifig in range(nfig):
         fig[ifig],sp[ifig,:,:] = plt.subplots(ny,nx,figsize=(3*nx*scale,3*ny*scale),sharex=True,sharey=True,squeeze=False,constrained_layout=True)
    
    for irun,run_dir in enumerate(run_dirs):
        if run_dir.encode() not in valid_rundirs:
            continue
    
        coords = np.zeros((n_mangles),dtype=int)
        for ii in range(n_mangles):
            if ii==n_mangles-1:
                coords[ii] = irun%mangle_n[ii]
            coords[ii] = (irun//np.prod(mangle_n[ii+1:]))%mangle_n[ii]
        ix = coords[iaxis]
        iy = coords[jaxis]
        ic = coords[color_axis]
        color_label = key_bases[color_axis][ic]

        for ifig,ir in enumerate(plot_r):
            print(ifig,ix,iy,sp.shape)
            sp[ifig,iy,ix].plot(df["phi"],df["{}_{}".format(run_dir,ir)],label=color_label,color=colors[ic])

    x = np.linspace(0.,90.,100)
    y = np.zeros(x.shape)
    for ifig,ir in enumerate(plot_r):
        for ix in range(nx):
            sp[ifig,0,ix].set_title(v_labels[ix])
        for iy in range(ny):
            sp[ifig,iy,nx-1].yaxis.set_label_position("right")
            sp[ifig,iy,nx-1].set_ylabel(mdot_labels[iy])
#                 sp[ifig,iy,ix].set_title(key_bases[jaxis][iy]+"_"+key_bases[iaxis][ix])
#                 sp[ifig,iy,ix].set_title(mdot_labels[iy]+", "+v_labels[ix])
        sp[ifig,-1,-1].legend(loc='best',fontsize='xx-small')
#         sp[ifig,-1,-1].legend(loc='best',fontsize='x-small')
#         sp[ifig,0,0].legend(loc='best',fontsize='xx-small')
#         sp[ifig,0,0].legend(loc='best',fontsize='x-small')
        for sp_row in sp[ifig,:,:]:
            for p in sp_row:
#                 p.legend(loc='best',fontsize='xx-small')
                p.plot(x,y,linestyle=':')
#                 p.set_xticks(np.linspace(0,90,10))
                p.set_xticks( np.linspace(0,90,7) )
                if mode=='p':
                    p.set_ylim(-1.e5,1.e5)
                    p.set_yscale("symlog",linthreshy=100.)
                elif mode=='v':
                    p.set_ylim(ytick_values[0],ytick_values[-1])
                    p.set_yscale("symlog",linthreshy=50.)
                    p.set_yticks(ytick_values)
                    p.set_ylim([-50.,700.])
                    p.set_yticklabels(ytick_values)
        fig[ifig].suptitle(r"$%g\mathrm{pc}\leq r<%g\mathrm{pc}$" % (r_bin_edges[ir],r_bin_edges[ir+1]) )
#         fig[ifig].text(0.5, 0.04, r"$\phi$ ($^\circ$)", ha='center')
        sp[ifig,ny-1,0].set_xlabel(r"$\phi$ ($^\circ$)")
        if mode=='p':
#             fig[ifig].text(0.04, 0.5, r"$\Sigma m v_r$ (M$_\odot$ km s$^{-1}$)", va='center', rotation='vertical')
            sp[ifig,ny-1,0].set_ylabel(r"$\Sigma m v_r$ (M$_\odot$ km s$^{-1}$)")
        elif mode=='v':
#             fig[ifig].text(0.04, 0.5, r"$v_r$ (km s$^{-1}$)", va='center', rotation='vertical')
            sp[ifig,ny-1,0].set_ylabel(r"$v_r$ (km s$^{-1}$)")
#         fig[ifig].savefig("../figures/vrad_rad_phi_testflows_r{}.png".format(ifig),dpi=150)
        fig[ifig].savefig("../figures/{}rad_rad_phi_testflows_r{}.pdf".format(mode,ifig),dpi=150)
    plt.close('all')
        