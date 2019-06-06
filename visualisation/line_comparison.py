    # Data produced with:
    # python sph_anim.py 3032 $i --rad 50. --savemap --view side --noring --snap0 100 --maxsnapf 100 --plot co1m,co2m,hcn1m,hcn2m,h2_1m,h2_2m,h2_3m

    # newer runs:
    # python sph_anim.py 3032 $i --rad 100. --savemap --view side --noring --snap0 172 --maxsnapf 172 --plot co1m,co2m,hcn1m,hcn2m,h2_1m,h2_2m,h2_3m

    # redone:
    # python sph_anim.py 3032 $i --rad 100. --savemap --view side --noring --snap0 100 --maxsnapf 100 --plot co1m,co2m,hcn1m,hcn2m,h2_1m,h2_2m,h2_3m,view

    # ok, actual one:
    #for i in longrun_weakflow_rapid_defaultaniso longrun_weakflow_rapid_defaultaniso_polar longrun_weakflow_settled_defaultaniso longrun_weakflow_settled_defaultaniso_polar longrun_weakflow_vesc_defaultaniso longrun_weakflow_vesc_defaultaniso_polar
#just long #for i in longrun_weakflow_rapid_defaultaniso_polar longrun_weakflow_settled_defaultaniso_polar longrun_weakflow_vesc_defaultaniso_polar longrun_weakflow_vesc_defaultaniso
#IR    #do python sph_anim.py 3032 $i --rad 10.,100.,10.,100. --savemap --view side --noring --snap0 100 --maxsnapf 100 --plot view,view,view,view --gaussian=-1.,-1.,1.,10.
#lines    #do python sph_anim.py 3032 $i --rad 20. --savemap --view side --noring --snap0 100 --maxsnapf 100 --plot co1m,co2m,hcn1m,hcn2m,h2_1m,h2_2m,h2_3m,co1m,co2m,hcn1m,hcn2m,h2_1m,h2_2m,h2_3m --gaussian=-1.,-1.,-1.,-1.,-1.,-1.,-1.,4.,4.,4.,4.,4.,4.,4.
    #done

print("Importing")

# Import libraries to do our magic
import numpy as np
import sys
import os
#from joblib import Parallel, delayed

#from sys import path
#path.append("../src/")
import sph_frame

from sys import path
path.append("../")
import gizmo_tools

import matplotlib.pyplot as plt

import rgb_image

print("Running")

output_dirs = [ "longrun_weakflow_rapid_defaultaniso",
                "longrun_weakflow_rapid_defaultaniso_polar",
                "longrun_weakflow_settled_defaultaniso",
                "longrun_weakflow_settled_defaultaniso_polar",
                "longrun_weakflow_vesc_defaultaniso",
                "longrun_weakflow_vesc_defaultaniso_polar"]
rads = [20.]*14+[10.,100.,10.,100.]

idump = 100

# output_dirs = [ "longrun_weakflow_rapid_defaultaniso_polar",
#                 "longrun_weakflow_settled_defaultaniso_polar",
#                 "longrun_weakflow_vesc_defaultaniso_polar"]
# rad = 100.

line_codes = ["co1m","co2m","hcn1m","hcn2m","h2_1m","h2_2m","h2_3m"]*2+["view"]*4


run_id = "3032"


nruns = len(output_dirs)
nlines = len(line_codes)# len(rgb_image.labels)

cmap = "inferno"
# cmap = "bwr"

nx = 2
ny = nruns//2
if ny*2<nruns: ny+=1

vranges = { "view":[0.,5.],
            "co1m":[-14.,-2.],
            "co2m":[-14.,-2.],
            "hcn1m":[-20.,-2.],
            "hcn2m":[-20.,-2.],
            "h2_1m":[-8.,0.],
            "h2_2m":[-20.,-10.],
            "h2_3m":[-8.,-0.]
            }

# for iline,line_label in enumerate(rgb_image.labels):
# for iline,line_label in [[7,"view"]]:
for iline,line_code in enumerate(line_codes):
#     if line_code!="view": continue
    line_label = rgb_image.line_label_dict[line_code]
    print("plotting {}".format(line_label))
    fig,sp = plt.subplots(ny,1+nx,figsize=(4.*nx,3.*ny), gridspec_kw = {'width_ratios':[16]*nx+[1],'height_ratios':[16]*ny})

    for irow in range(ny-1):
        for icol in [0,1]:
            sp[irow,icol].set_xticklabels([])
            sp[irow,icol].set_xticklabels([])
        for icol in [1]:
            sp[irow,icol].set_yticklabels([])
        for icol in [0]:
            sp[irow,icol].get_shared_x_axes().join(sp[irow,icol], sp[irow+1,icol])
            sp[irow,icol].get_shared_x_axes().join(sp[irow,icol], sp[irow,icol+1])

    for irow in range(ny):
        for icol in [0,1]:
            sp[irow,icol].set(adjustable='box-forced', aspect='equal')

    plot_str = line_code
    if line_codes.count(plot_str)>0:
        plot_str+="_%03d"%line_codes[:iline].count(line_code)

    if line_code in vranges:
        vrange = vranges[line_code]
    else:
        vrange = None

    for irow,output_dir in enumerate(output_dirs):
        ix = irow%2
        iy = irow//2
#         basefile = "../data/smoothsum_co1mco2mhcn1mhcn2mh2_1mh2_2mh2_3mgiz_{}_{}_000_%03d.dat".format(run_id,output_dir)
#         basefile = "../data/smoothsum_co1mco2mhcn1mhcn2mh2_1mh2_2mh2_3mgiz_{}_{}_000_%03d.dat".format(run_id,output_dir)

        basefile = "../data/smoothsum_co1mco2mhcn1mhcn2mh2_1mh2_2mh2_3mviewgiz_{}_{}_000_%03d.dat".format(run_id,output_dir)
        infile = "../data/smoothsum_giz_{}_{}_{:03d}_{}.dat".format(run_id,output_dir,idump,plot_str)
        infile = [infile]

        label = [line_label]
        print("..subplot {}".format(infile[0]))
        
        ax = sp[iy,ix]
        
        mappable=rgb_image.produce_and_save_rgb_table_image(infile,ax=ax,rad=rads[iline],labels=label,vrange=vrange,cmap='inferno')
        ax.set_title(output_dir)
        if ix==nx-1:
            plt.colorbar(mappable,cax=sp[iy,nx])

    if nx*ny>nruns:
        sp[ny-1,nx].remove()
        sp[ny-1,nx-1].remove()

    # for irow in range(ncuts):
    #     sp[irow,0].yaxis.tick_left()
    #     sp[irow,0].yaxis.set_label_position("left")

    sp[ny-1,0].set_xlabel("pc")
    sp[ny-1,0].set_ylabel("pc")

    # fig.subplots_adjust(hspace=0., wspace=0.) 
#     fig.subplots_adjust(left=0.1,right=.9,bottom=0.05,top=.99,hspace=0.03,wspace=0.01) 
    fig.tight_layout()
    print(".finishing plot")

    plt.savefig("../../figures/line_montage_{}.png".format(plot_str),dpi=150)
