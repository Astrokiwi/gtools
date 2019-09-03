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
import config3032

config3032.setup("../")

# output_dirs = ["longrun_medflow_settled_defaultaniso_polar","longrun_medflow_vesc_defaultaniso","longrun_medflow_vesc_defaultaniso_polar","longrun_weakflow_rapid_defaultaniso","longrun_weakflow_rapid_defaultaniso_polar","longrun_weakflow_settled_defaultaniso","longrun_weakflow_settled_defaultaniso_polar","longrun_weakflow_vesc_defaultaniso","longrun_weakflow_vesc_defaultaniso_polar","newflow_settled_thin_up","newflow_vesc_thin_45","newflow_vesc_thin_side","newflow_vesc_thin_up"]
# output_dirs = ["longrun_medflow_settled_defaultaniso_polar","longrun_medflow_vesc_defaultaniso","longrun_medflow_vesc_defaultaniso_polar","longrun_weakflow_rapid_defaultaniso","longrun_weakflow_rapid_defaultaniso_polar","longrun_weakflow_settled_defaultaniso","longrun_weakflow_settled_defaultaniso_polar","longrun_weakflow_vesc_defaultaniso","longrun_weakflow_vesc_defaultaniso_polar","newflow_settled_thin_up","newflow_vesc_thin_45","newflow_vesc_thin_side","newflow_vesc_thin_up"]
# output_dirs = ["newflow_vesc_thin_45","longrun_medflow_vesc_defaultaniso_polar"] # for presentation
# output_dirs = [ "longrun_weakflow_rapid_defaultaniso",
#                 "longrun_weakflow_rapid_defaultaniso_polar",
#                 "longrun_weakflow_settled_defaultaniso",
#                 "longrun_weakflow_settled_defaultaniso_polar",
#                 "longrun_weakflow_vesc_defaultaniso",
#                 "longrun_weakflow_vesc_defaultaniso_polar",
#                 "longrun_medflow_vesc_defaultaniso"]

# everything, for analysis
# output_dirs = [
# "longrun_medflow_settled_defaultaniso",
# "longrun_medflow_settled_defaultaniso_polar",
# "longrun_medflow_vesc_defaultaniso",
# "longrun_medflow_vesc_defaultaniso_polar",
# "longrun_weakflow_rapid_defaultaniso",
# "longrun_weakflow_rapid_defaultaniso_polar",
# "longrun_weakflow_settled_defaultaniso",
# "longrun_weakflow_settled_defaultaniso_polar",
# "longrun_weakflow_vesc_defaultaniso",
# "longrun_weakflow_vesc_defaultaniso_polar",
# "newflow_settled_thin_up",
# "newflow_vesc_thin_45",
# "newflow_vesc_thin_side",
# "newflow_vesc_thin_up"]


# samples, for production
output_dirs = ["longrun_medflow_vesc_defaultaniso_polar","newflow_vesc_thin_45"]


# rads = [20.]*7*2+[10.,100.,10.,100.]
# rads = [20.]*3*2+[10.,100.,10.,100.]

# idump = 100

# output_dirs = [ "longrun_weakflow_rapid_defaultaniso_polar",
#                 "longrun_weakflow_settled_defaultaniso_polar",
#                 "longrun_weakflow_vesc_defaultaniso_polar"]
# rad = 100.

run_id = "3032"

# line_codes = ["co1m","co2m","hcn1m","hcn2m","h2_1m","h2_2m","h2_3m"]*2+["view"]*4

# line_codes = ["hcn2m_001","co2m_001","h2_1m_001","view_002","view_003"] # for presentation
# rads = [20.]*3+[10.,100.]

line_codes = ["co1m_000","co2m_000","hcn1m_000","hcn2m_000","h2_1m_000","h2_2m_000","h2_3m_000","view_000","co1m_001","co2m_001","hcn1m_001","hcn2m_001","h2_1m_001","h2_2m_001","h2_3m_001","view_001"]
# rads = [20.]*7+[10.,100.]
rads = [10.]*8+[100.]*8

# line_labels = [r"HCN(8-7), $422.796$ $\mu$m",r"CO(6-5), $433.438$ $\mu$m",r"H$_2$ (1-0) S(1), $2.121$ $\mu$m",r"$F_{IR}$",r"$F_{IR}$"] # for presentation
line_name_lookup = gizmo_tools.lineNamesFull.copy()
line_name_lookup["vie"] = r"$F_{IR}$"
line_labels = [line_name_lookup[line[:-5]] for line in line_codes]

# clear
# vranges = { "view":[0.,5.],
#             "co1m":[-14.,-2.],
#             "co2m":[-14.,-2.],
#             "hcn1m":[-20.,-2.],
#             "hcn2m":[-20.,-2.],
#             "h2_1m":[-8.,0.],
#             "h2_2m":[-20.,-10.],
#             "h2_3m":[-8.,-0.]
#             }

# fair
vranges = { "view":[0.,5.],
            "co1m":[-20.,0.],
            "co2m":[-20.,0.],
            "hcn1m":[-20.,0.],
            "hcn2m":[-20.,0.],
            "h2_1m":[-20.,0.],
            "h2_2m":[-20.,0.],
            "h2_3m":[-20.,0.]
            }

nruns = len(output_dirs)
nlines = len(line_codes)# len(rgb_image.labels)

cmap = "inferno"
# cmap = "bwr"

def render_all(idump):
    print("Rendering")
    snap_str = "{:03d}".format(idump)
    for output_dir in output_dirs:
        commands = []
#         commands.append("python sph_anim.py 3032 {run_name} --rad 10.,100.,10.,100. --savemap --view side --noring --snap0 {snap} --maxsnapf {snap} --plot view,view --gaussian=1.,10. &")
# #         commands.append("python sph_anim.py 3032 {run_name} --rad 20. --savemap --view side --noring --snap0 {snap} --maxsnapf {snap} --plot co1m,co2m,hcn1m,hcn2m,h2_1m,h2_2m,h2_3m,co1m,co2m,hcn1m,hcn2m,h2_1m,h2_2m,h2_3m --gaussian=-1.,-1.,-1.,-1.,-1.,-1.,-1.,4.,4.,4.,4.,4.,4.,4. &")
#         commands.append("python sph_anim.py 3032 {run_name} --rad 10.,10.,10.,10.,10.,10.,10.,100.,100.,100.,100.,100.,100.,100. --savemap --view side --noring --snap0 {snap} --maxsnapf {snap} --plot co1m,co2m,hcn1m,hcn2m,h2_1m,h2_2m,h2_3m,co1m,co2m,hcn1m,hcn2m,h2_1m,h2_2m,h2_3m --gaussian=1.,1.,1.,1.,1.,1.,1.,10.,10.,10.,10.,10.,10.,10. &")

        commands.append("python sph_anim.py 3032 {run_name} --rad 10.,100.,10.,100. --savemap --view side --noring --snap0 {snap} --maxsnapf {snap} --plot view,view &")
        commands.append("python sph_anim.py 3032 {run_name} --rad 10.,10.,10.,10.,10.,10.,10.,100.,100.,100.,100.,100.,100.,100. --savemap --view side --noring --snap0 {snap} --maxsnapf {snap} --plot co1m,co2m,hcn1m,hcn2m,h2_1m,h2_2m,h2_3m,co1m,co2m,hcn1m,hcn2m,h2_1m,h2_2m,h2_3m &")

        for command in commands:
            formatted_command = command.format(run_name=output_dir,snap=snap_str)
            print(formatted_command)
            os.system(formatted_command)

def plot_lines_separately():
    nx = 2
    ny = nruns//nx
    if ny*nx<nruns: ny+=1

    # for iline,line_label in enumerate(rgb_image.labels):
    # for iline,line_label in [[7,"view"]]:
    for iline,line_code in enumerate(line_codes):
    #     if line_code!="view": continue
        line_label = rgb_image.line_label_dict[line_code]
        print("plotting {}".format(line_label))
        fig,sp = plt.subplots(ny,1+nx,figsize=(4.*nx+0.25*nx,3.*ny), gridspec_kw = {'width_ratios':[16]*nx+[1],'height_ratios':[16]*ny},squeeze=False)

        for irow in range(ny-1):
            for icol in range(nx):
                sp[irow,icol].set_xticklabels([])
            for icol in range(1,nx):
                sp[irow,icol].set_yticklabels([])
            for icol in range(nx-1):
                sp[irow,icol].get_shared_x_axes().join(sp[irow,icol], sp[irow+1,icol])
                sp[irow,icol].get_shared_x_axes().join(sp[irow,icol], sp[irow,icol+1])

        for irow in range(ny):
            for icol in range(nx):
                sp[irow,icol].set(adjustable='box-forced', aspect='equal')

        plot_str = line_code
        if line_codes.count(plot_str)>0:
            plot_str+="_%03d"%line_codes[:iline].count(line_code)

        if line_code in vranges:
            vrange = vranges[line_code]
        else:
            vrange = None

        for irow,output_dir in enumerate(output_dirs):
            ix = irow%nx
            iy = irow//nx
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
            ax.set_xlabel("")
            ax.set_ylabel("")
            if ix==0:
                plt.colorbar(mappable,cax=sp[iy,nx])

        if nx*ny>nruns:
            for iplot in range(nruns,nx*ny):
                ix = iplot%nx
                iy = iplot//nx
#                 print(ix,nx,iy,ny,nruns,iplot,nx*ny)
                sp[iy,ix].remove()

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

def which_files_exist(run_id,output_dirs,idump):
    output_dirs_which_exist = []
    for output_dir in output_dirs:
        file_exists = True
        for line_code in line_codes:
            infile = "../data/smoothsum_giz_{}_{}_{:03d}_{}.dat".format(run_id,output_dir,idump,line_code)
            if not os.path.isfile(infile):
                file_exists = False
#                 print(infile," n'existe pas")
                break
        if file_exists:
            output_dirs_which_exist.append(output_dir)
    return output_dirs_which_exist

def plot_lines_together(run_id,output_dirs,idump,split=False):
    ny = len(output_dirs)
    nx = len(line_codes)
    if split:
        nx//=2
        ny*=2
    fig,sp = plt.subplots(ny+1,nx,figsize=(3.*nx,3.*ny+3./16.), gridspec_kw = {'width_ratios':[16]*nx,'height_ratios':[16]*ny+[1]},squeeze=False,constrained_layout=True)
#     for irow in range(ny-1):
#         for icol in range(nx-1):
#             sp[irow,icol].get_shared_x_axes().join(sp[irow,icol], sp[irow+1,icol])
#             sp[irow,icol].get_shared_x_axes().join(sp[irow,icol], sp[irow,icol+1])
#             sp[irow,icol].get_shared_y_axes().join(sp[irow,icol], sp[irow+1,icol])
#             sp[irow,icol].get_shared_y_axes().join(sp[irow,icol], sp[irow,icol+1])

    for irow in range(ny):
        for icol in range(nx):
            sp[irow,icol].set(adjustable='box-forced', aspect='equal')


    for iline,line_code in enumerate(line_codes):
        for idir,output_dir in enumerate(output_dirs):
#             basefile = "../data/smoothsum_co1mco2mhcn1mhcn2mh2_1mh2_2mh2_3mviewgiz_{}_{}_000_%03d.dat".format(run_id,output_dir)
            infile = "../data/smoothsum_giz_{}_{}_{:03d}_{}.dat".format(run_id,output_dir,idump,line_code)
            infile = [infile]

            print("..subplot {}".format(infile[0]))
            if split:
                ix = iline%(len(line_codes)//2)
                iy = iline//(len(line_codes)//2)+idir*2
            else:
                ix = iline
                iy = idir
        
            ax = sp[iy,ix]
            if line_code in vranges:
                vrange = vranges[line_code]
            elif line_code[:-4] in vranges:
                vrange = vranges[line_code[:-4]]
            else:
                vrange = None
            
            mappable=rgb_image.produce_and_save_rgb_table_image(infile,ax=ax,rad=rads[iline],vrange=vrange,cmap='inferno')
            ax.set_title("")
            if iy==ny-1 or split and iy%2==1:
                ax.set_xlabel("pc")
            else:
                ax.set_xlabel("")
                ax.xaxis.set_visible(False)
#             if ix==0:
#                 ax.set_ylabel("pc")
#             else:
#                 ax.set_ylabel("")
#                 ax.yaxis.set_visible(False)
            if iy==0:
#                 sp[iy,ix].set_title(line_labels[ix])
                ax.set_title(line_labels[iline])
                plt.colorbar(mappable,cax=sp[ny,ix],orientation='horizontal')
            if ix==0:
#                 ax.set_ylabel(output_dir,fontsize=6)
                ax.set_ylabel(config3032.run_parameters[output_dir]["name"])
            else:
                ax.yaxis.set_visible(False)
                ax.set_ylabel("")

#     plt.savefig("../../figures/lines_together_ewass.png",dpi=150)
#     plt.savefig("../../figures/lines_together_analysis_{}.png".format(idump),dpi=150)
    plt.savefig("../../figures/lines_together_summary_{}.png".format(idump),dpi=150)

if __name__=='__main__':
    print("Running")
#     render_all(100)
#     render_all(200)
#     render_all(300)
#     idumps = [100,200,300]
    idumps = [100]


    if len(sys.argv)>=2:
        for idump in idumps:
            render_all(idump)
    else:
        for idump in idumps:
            plot_lines_together(run_id,which_files_exist(run_id,output_dirs,idump),idump,split=True)
#         print(which_files_exist(run_id,output_dirs,100))
#         plot_lines_separately()
#         plot_lines_together(idump)
