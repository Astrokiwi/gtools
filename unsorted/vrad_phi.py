import numpy as np
import matplotlib.pyplot as plt

from unsorted import parameter_mangler
import gizmo_tools
    
if __name__=='__main__':
    run_id = "2030"

#     snap_str = "010"
#     for snap_str in "001","004","010"
#     snap_str = "001"
        
    manglia = parameter_mangler.load_manglia("prams/pram_options_testflows.dat")
    key_combinations = parameter_mangler.mangle_combinations(manglia)
#     print(key_combinations)
    run_dirs = [parameter_mangler.filebase_from_combination("testflows", x) for x in key_combinations]
#     print(run_dirs)
    key_bases = [list(x.mangle_texts.keys()) for x in manglia]
    mangle_n = [len(x) for x in key_bases]
    n_mangles = len(mangle_n)
    
    iaxis = 1
    nx = mangle_n[iaxis]
    jaxis = 0
    ny = mangle_n[jaxis]
    color_axis = 2
    colors = ['r','g','b','k']
#     linestyles = ['dashed','dashdot','solid']

    for snap_str in "001","004","010":
        fig,sp = plt.subplots(ny,nx,figsize=(4*nx,4*ny),sharex=True,sharey=True,squeeze=False)
        phi_bin_edges = np.linspace(0.,90.,19)
        phi_centres = (phi_bin_edges[:-1]+phi_bin_edges[1:])/2.
    #     for irun,run_dir in enumerate(run_dirs):
        for irun,run_dir in enumerate(run_dirs[:34]):
            try:
                coords = np.zeros((n_mangles),dtype=int)
                for ii in range(n_mangles):
                    if ii==n_mangles-1:
                        coords[ii] = irun%mangle_n[ii]
                    coords[ii] = (irun//np.prod(mangle_n[ii+1:]))%mangle_n[ii]
    #             iplot = coords[plot_axis]
                ix = coords[iaxis]
                iy = coords[jaxis]
                ic = coords[color_axis]
                color_label = key_bases[color_axis][ic]
    #             color_label = run_dir
                print(coords,iaxis,jaxis,ix,iy,sp.shape)
                t,data=gizmo_tools.load_gizmo_pandas(run_id,run_dir,snap_str,["Masses","Coordinates","Velocities","AGNIntensity"])
                print(run_dir,"found")
                gizmo_tools.calculate_phi(data)
                gizmo_tools.calculate_vrad(data)
                bin_N,x = np.histogram(data["phi"],phi_bin_edges)
                bin_vrad_sum,x = np.histogram(data["phi"],phi_bin_edges,weights=data["vrad"])
    #             bin_vrad_sum,x = np.histogram(data["phi"],phi_bin_edges,weights=data["AGNIntensity"]*data["rad3d"]**2)
                bin_vrad_mean = bin_vrad_sum/bin_N
                bin_vrad_mean*=np.cos(np.radians(phi_centres))
                bin_vrad_mean=np.cumsum(bin_vrad_mean)
                print(ix,iy,sp.shape)
                sp[ix,iy].plot(phi_centres,bin_vrad_mean,label=color_label,color=colors[ic])
            except OSError:
                print(run_dir," not found")
                continue
        x = np.linspace(0.,90.,100)
        y = np.zeros(x.shape)
        for ix in range(nx):
            for iy in range(ny):
                sp[ix,iy].set_title(key_bases[jaxis][iy]+"_"+key_bases[iaxis][ix])
        sp[0,0].legend(loc='best',fontsize='xx-small')
        ytick_values = [-100,-70,-50,-20,0,20,50,70,100,200,500,700,1000]
        for sp_row in sp:
            for p in sp_row:
                p.legend(loc='best',fontsize='xx-small')
                p.plot(x,y,linestyle=':')
                p.set_xticks(phi_bin_edges)
                p.set_yscale("symlog",linthreshy=100.)
                p.set_yticks(ytick_values)
                p.set_yticklabels(ytick_values)
        plt.savefig("../figures/vrad_phi_testflows{}.png".format(snap_str),dpi=150)
        plt.close()
        