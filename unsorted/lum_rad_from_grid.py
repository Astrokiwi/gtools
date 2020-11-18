print("Importing")

import numpy as np

import gizmo_tools

import matplotlib.pyplot as plt

from scipy import interpolate

import multiprocessing
import functools
import subprocess
import sys

# extinguished_line_codes = ["view"+line_basis for line_basis in line_bases]
# unextinguished_line_codes = [line_basis+"m" for line_basis in line_bases]

rad_lines = ["IRdust"]+gizmo_tools.line_bases

# everything, for analysis
output_dirs = [
"longrun_medflow_settled_defaultaniso",
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

# samples for some plots
# output_dirs = ["longrun_medflow_vesc_defaultaniso_polar","newflow_vesc_thin_45"]

run_id = "3032"
lum_factors = np.array([0.5,0.9])


resolution = 2048

def extinguished_formatter(line_basis):
    return "view"+line_basis

def unextinguished_formatter(line_basis):
    return line_basis+"m"

def render_all(idump,nprocs=64):
    print("Rendering")
    snap_str = "{:03d}".format(idump)
    commands = []
    for output_dir in output_dirs:
# no extinction, include IR
        commands.append("cd visualisation; python sph_anim.py 3032 {run_name} --rad 10.,10.,10.,10.,10.,10.,10.,10.,100.,100.,100.,100.,100.,100.,100.,100. --savemap --view side --noring --snap0 {snap} --maxsnapf {snap} --plot IRdustm,co1m,co2m,hcn1m,hcn2m,h2_1m,h2_2m,h2_3m,IRdustm,co1m,co2m,hcn1m,hcn2m,h2_1m,h2_2m,h2_3m --suffix grid_unextinguished --L {resolution} &".format(run_name=output_dir,snap=snap_str,resolution=resolution))
#extinction, include IR
        commands.append("cd visualisation; python sph_anim.py 3032 {run_name} --rad 10.,10.,10.,10.,10.,10.,10.,10.,100.,100.,100.,100.,100.,100.,100.,100. --savemap --view side --noring --snap0 {snap} --maxsnapf {snap} --plot viewIRdust,viewco1,viewco2,viewhcn1,viewhcn2,viewh2_1,viewh2_2,viewh2_3,viewIRdust,viewco1,viewco2,viewhcn1,viewhcn2,viewh2_1,viewh2_2,viewh2_3  --suffix grid_extinguished --L {resolution}  &".format(run_name=output_dir,snap=snap_str,resolution=resolution))

    with multiprocessing.Pool(nprocs) as pool:
        pool.map(functools.partial(subprocess.call,shell=True),commands)


def grids_to_radially_cumulative(run_name,snap="100",suffix="001",code_formatter=extinguished_formatter,radius=100.,fig_suffix="extinguished",doplot=False):
    try:
        filestart = f"data/smoothsum_giz_{run_id}_{run_name}grid_{fig_suffix}_{snap}_"
        fileend = f"_{suffix}.dat"
    
        if doplot:
            fig,sp = plt.subplots(figsize=(8,8))
    
            ax = sp
        
        summary_outp=[]

        for line_basis in rad_lines:
            map = np.loadtxt(filestart+code_formatter(line_basis)+fileend)
        
            flat_map = map.flatten()
        
            flat_map_radii = np.linalg.norm(
                                np.meshgrid(np.linspace(-radius,radius,map.shape[0]),
                                            np.linspace(-radius,radius,map.shape[1])),
                                            axis=0)\
                                .flatten()
        
            sort_indices = np.argsort(flat_map_radii)

            nice_slice = np.isfinite(flat_map[sort_indices])
            sort_indices = sort_indices[nice_slice]
        
            cum_lum = np.cumsum(10**flat_map[sort_indices])
            cum_lum/=np.max(cum_lum)
        
            radlum_interp = interpolate.interp1d(cum_lum,flat_map_radii[sort_indices])
        
            r_lums = radlum_interp(lum_factors)
            print(r_lums,map.size)
            
            if map.size!=resolution**2:
                raise ValueError("INCORRECT MAP SIZE {} {} {} {} ".format(map.shape,map.size,resolution,resolution**2))
        
    #         ax.scatter(flat_map_radii[sort_indices],flat_map[sort_indices],label=line_basis)
            if doplot:
                ax.plot(flat_map_radii[sort_indices],cum_lum,label=line_basis)
            summary_outp+=[r_lums]
        np.savetxt("data/summary_lumrads_{}_{}.dat".format(run_name,fig_suffix),summary_outp)
        
        if doplot:
            ax.legend()
            ax.set_xscale('log')
            ax.set_ylim([1.e-4,1.1])
            ax.set_yscale('log')
        #     fig.savefig("../figures/lum_rad_from_grid_{}.png".format(fig_suffix),dpi=200)
            fig.savefig("../figures/lum_rad_from_grid_{}.pdf".format(fig_suffix),dpi=200)
            plt.close('all')
    except OSError as e:
        print("File not found "+str(e))
        
print("loading and plotting")

if __name__=='__main__':
    print("Running")
    idump = 100
    if len(sys.argv)>=2:
        render_all(idump)
    else:
        nprocs=172
        with multiprocessing.Pool(nprocs) as pool:
            pool.map(functools.partial(grids_to_radially_cumulative,snap=idump),output_dirs)
            pool.map(functools.partial(grids_to_radially_cumulative,code_formatter=unextinguished_formatter,fig_suffix="unextinguished",snap=idump),output_dirs)
