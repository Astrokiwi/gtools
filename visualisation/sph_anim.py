print("Importing")

import tblib.pickling_support
tblib.pickling_support.install()

# Import libraries to do our magic
import numpy as np
import sys
import os
# from joblib import Parallel, delayed
from multiprocessing import Pool
import functools

import sph_frame
import re

from sys import path
path.append("../")
import gizmo_tools

import argparse


# def sort_nicely( l ):
#     """ Sort the given list in the way that humans expect.
#     """
#     convert = lambda text: int(text) if text.isdigit() else text
#     alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
#     l.sort( key=alphanum_key )

def dump_rad0(infile):
    import h5py
    with h5py.File(infile,"r") as f:
        xyz_p = np.array(f["/PartType0/Coordinates"]) # kpc
        id_p = np.array(f["/PartType0/ParticleIDs"]).astype(int)
    id_p-=1
    rad_p = np.sqrt(xyz_p[:,0]**2+xyz_p[:,1]**2+xyz_p[:,2]**2)*1.e3
    rad_out = np.zeros_like(rad_p)
    rad_out[id_p] = rad_p
    np.save("rad0.npy",rad_out)
#     sys.exit()
    
def string_to_list_or_float(s):
    if s is None:
        return None
    if "," in s:
        return [float(x) for x in s.split(",")]
    return float(s)

# joblib doesn't dump exceptions well, just print them
def plotter_parallel_exception_wrapper(*args,**kwargs):
    try:
        return sph_frame.makesph_trhoz_frame(*args,**kwargs)
    except Exception as e:
        print(e) # should give some output at least
        raise(e) # might not work

if __name__ == '__main__':
    default_values = dict()
    default_values["nprocs"]=8
    default_values["maxsnapf"]=-1
    default_values["snap0"]=-1
    default_values["rad"]="15"
    #default_values["rads"]=""
    default_values["L"]=400
    default_values["plot"]="dens,temp"
    default_values["views"]="face,side"
    default_values["cmap"]="viridis"
    default_values["slice"]=False
    default_values["phi"]=0.
    default_values["theta"]=0.
    default_values["noaxes"]=False
    default_values["centredens"]=False
    default_values["centrecom"]=False
    default_values["suffix"]=""
    default_values["dotmode"]=""
    default_values["absurd"]=False
    default_values["pixsize"]=1
    default_values["noring"]=False
    default_values["data_ranges"]=""
    default_values["savemap"]=False
    default_values["gaussian"]=None
    
    parsevals = ["data_ranges","nprocs","maxsnapf","run_id","output_dir","plot","cmap","rad","L","slice","views","phi","theta","noaxes","centredens","centrecom","suffix","dotmode","absurd","pixsize","noring","snap0","savemap","gaussian"]

    parser = argparse.ArgumentParser()
    parser.add_argument('run_id',help="name of superdirectory for runs")
    parser.add_argument('output_dir',help="name of subdirectory for run")
    parser.add_argument('--nprocs',type=int,help="processors to run on (default {})".format(default_values["nprocs"]))
    parser.add_argument('--maxsnapf',type=int,help="snapshot to end on (default=-1=do all snapshots)")
    parser.add_argument('--snap0',type=int,help="snapshot to start on (default=-1=do all snapshots)")
    parser.add_argument('--rad',type=str,help="radius of all plots in parsecs - same value for all plots, or separated by commas")
    #parser.add_argument('--rads',type=float,help="radius of each plot in parsecs, separated by commas")
    parser.add_argument('--L',type=int,help="size of plot area in pixels")
    parser.add_argument('--plot',type=str,help="values to plot, separated by commas")
    parser.add_argument('--views',type=str,help="face and/or side view, separated by commas")
    parser.add_argument('--cmap',type=str,help="colourmap palette")
    parser.add_argument('--phi',type=float,help="phi (latitudinal) rotation angle in degrees")
    parser.add_argument('--theta',type=float,help="theta (azimuthal) rotation angle in degrees")
    parser.add_argument('--slice',help="option to be slice plot",action='store_true')
    parser.add_argument('--noaxes',help="don't show axes - fill frame with plot",action='store_true')
    parser.add_argument('--centredens',help="centre on densest particle",action='store_true')
    parser.add_argument('--centrecom',help="centre on centre of mass",action='store_true')
    parser.add_argument('--suffix',help="suffix on anim filename")
    parser.add_argument('--dotmode',help="max, min, or nothing - dot plot mode")
    parser.add_argument('--pixsize',type=int,help="pixelation level")
    parser.add_argument('--data_ranges',type=str,help="explicit data bounds for plot - format: --data_ranges=min1,max1+min2,max2+min3,max3 etc; the equals sign may be necessary!")
    parser.add_argument('--absurd',help="I'll try spinning, that's a good trick",action='store_true')
    parser.add_argument('--noring',help="Edge on ring, not wrapped",action='store_true')
    parser.add_argument('--savemap',help="Save maps as data file",action='store_true')
    parser.add_argument('--gaussian',type=str,help="Gaussian filter size - same value for all plots, or separated by commas")
    args = parser.parse_args()
    
#     run_id = args.run_id
#     output_dir = args.output_dir
    
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

    flatPlot = not slice
    ringPlot = flatPlot
    if noring:
        ringPlot = False

    if dotmode:
        smooth_str = dotmode
    elif flatPlot:
        smooth_str = "smooth"
    else:
        smooth_str = "slice"
    
    rad = string_to_list_or_float(rad)
    gaussian = string_to_list_or_float(gaussian)
    
    if len(data_ranges)<=0:
        data_ranges = None
    else:
        data_ranges=[[float(y) for y in x.split(",")] for x in data_ranges.split("+")]
    
    toplot = plot.split(",")
    outp_plot = "".join(toplot)

    toview = views.split(",")
    outp_views = "".join(toview)
    
    phi*=2.*np.pi/360.
    theta*=2.*np.pi/360.
    
    visibleAxes=not noaxes
    print("Running")

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

    os.system("rm ../pics/sphplot"+run_id+output_dir+"???.png")


    print("nfiles:",snapf-snapi+1)
    
    isnaps = range(snapi,snapf+1)
    
    infiles = [fullDir+"/snapshot_"+("%03d" % snapx)+".hdf5" for snapx in isnaps]
    outfiles = ["../pics/sphplot"+run_id+output_dir+"%03d.png"%snapx for snapx in isnaps]
    
    nfiles = len(infiles)
    if absurd:
        period2 = 11
        phis = np.linspace(0.,np.pi/period2*nfiles,nfiles)
        thetas = np.linspace(0.,np.pi/period2*nfiles,nfiles)
        period_alt = 10
        ungles = np.linspace(0.,np.pi/period_alt*nfiles,nfiles)
        rads = np.sin(ungles)*rad/2.+rad
    else:
        thetas = np.full(nfiles,theta)
        phis = np.full(nfiles,phi)
        if isinstance(rad,list):
            rads = [rad]*nfiles
        else: # if float
            rads = np.full(nfiles,rad)
    
    if "rad0" in toplot:
        dump_rad0(infiles[0])

    def frame_i(i):
        return sph_frame.makesph_trhoz_frame(infiles[i],outfiles[i],cmap=cmap,flat=flatPlot,ring=ringPlot,plot=toplot,L=L,scale=rads[i],views=toview,rot=[thetas[i],phis[i]],visibleAxes=visibleAxes,centredens=centredens,centrecom=centrecom,dotmode=dotmode,pixsize=pixsize,data_ranges=data_ranges,return_maps=savemap,gaussian=gaussian)
    
    with Pool(processes=nprocs) as pool:
        maps = pool.map(frame_i,range(snapf+1-snapi))
#     maps=Parallel(n_jobs=nprocs)(delayed(sph_frame.makesph_trhoz_frame)(infiles[i],outfiles[i],cmap=cmap,flat=flatPlot,ring=ringPlot,plot=toplot,L=L,scale=rads[i],views=toview,rot=[thetas[i],phis[i]],visibleAxes=visibleAxes,centredens=centredens,centrecom=centrecom,dotmode=dotmode,pixsize=pixsize,data_ranges=data_ranges,return_maps=savemap,gaussian=gaussian) for i in range(snapf+1-snapi))
#     if len(maps)==1:
#         maps = [maps]
    for result in maps:
        if isinstance(result, sph_frame.ExceptionWrapper):
            result.re_raise()
    if savemap:
        for itime,maps_timeslice in enumerate(maps):
            idump = isnaps[itime]
            for iplot,map in enumerate(maps_timeslice):
                plot_str = toplot[iplot]
                if toplot.count(plot_str)>0:
                    plot_str+="_%03d"%toplot[:iplot].count(toplot[iplot])
                outfile = "../data/"+smooth_str+"sum_giz_"+run_id+"_"+output_dir+suffix+"_%03d"%idump+"_"+plot_str+".dat"
                np.savetxt(outfile,map)
    #Parallel(n_jobs=nprocs)(delayed(sph_frame.makesph_trhoz_frame)(infiles[i],outfiles[i],cmap='plasma',flat=True,ring=True,plot=['dens'],L=400) for i in range(snapi,snapf+1))
    #Parallel(n_jobs=nprocs)(delayed(sph_frame.makesph_trhoz_frame)(infiles[i],outfiles[i],cmap='viridis',flat=True,ring=True,plot=['dt'],L=200,scale=15.,pixsize=2) for i in range(snapi,snapf+1))
    #Parallel(n_jobs=nprocs)(delayed(sph_frame.makesph_trhoz_frame)(infiles[i],outfiles[i],cmap='plasma',flat=True,ring=True,plot=['dens','temp'],L=400,scale=15.) for i in range(snapi,snapf+1))
    #Parallel(n_jobs=nprocs)(delayed(sph_frame.makesph_trhoz_frame)(infiles[i],outfiles[i],cmap='plasma',flat=True,ring=True,plot=['vels','dens'],L=400) for i in range(snapi,snapf+1))
    #Parallel(n_jobs=nprocs)(delayed(sph_frame.makesph_trhoz_frame)(infiles[i],outfiles[i],cmap='plasma',flat=True,ring=True,plot=['vels'],L=400,scale=10.) for i in range(snapi,snapf+1))
    #[sph_frame.makesph_trhoz_frame(infiles[i],outfiles[i],cmap='plasma',flat=True,ring=True,plot=['dens'],L=400) for i in range(snapi,snapf+1)]

#     for snapx in range(snapi,snapf+1):
#         infile = fullDir+"/snapshot_"+("%03d" % snapx)+".hdf5"
#         outfile = "../pics/sphplot"+run_id+output_dir+"%03d.png"%snapx
#         #sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',flat=True,ring=True)
#         #sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',flat=False,ring=False,plot=['dens','temp'],L=400)
#         #sph_frame.makesph_trhoz_frame(infile,outfile,cmap='Greys',flat=True,ring=True,plot=['dens'],L=400)
#         sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',flat=True,ring=True,plot=['dens'],L=400)
#         #sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',flat=True,ring=True,plot=['temp','dens'],L=400)
#         #sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',flat=False,ring=False,plot=['temp','tdust'],L=400)
#         #sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',flat=True,ring=True,plot=['temp','tdust'],L=400)
#         #sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',flat=False,ring=False,plot=['dens','temp'],L=400)
#         #sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',flat=True,ring=True,plot=['vels','dens'],L=400,subsample=10,pixsize=1)   
#         #sph_frame.makesph_trhoz_frame(infile,outfile,cmap='Greys',flat=True,ring=True,plot=['emit','temp','dens'],L=200)

    
    #for snapx in [37]:
    
    print("to mp4!")
    cmd = "ffmpeg -y -r 24 -i ../pics/sphplot"+run_id+output_dir+"%03d.png -c:v mpeg4 -q:v 1 "+movieDir+"/"+smooth_str+"sum_"+outp_plot+"giz_"+run_id+"_"+output_dir+suffix+".mp4"


    #cmd = "ffmpeg -y -r 24 -i ../pics/sphplot%03d.png -c:v mpeg4 -q:v 1 /export/1/djw/movies/smooth_rhotempgiz_"+run_id+"_"+output_dir+".mp4"
    #cmd = "ffmpeg -y -r 24 -i ../pics/sphplot%03d.png -c:v mpeg4 -q:v 1 /export/1/djw/movies/smooth_depthtempgiz_"+run_id+"_"+output_dir+".mp4"

    #cmd = "ffmpeg -y -r 24 -i ../pics/sphplot"+run_id+output_dir+"%03d.png -c:v mpeg4 -q:v 1 /export/1/djw/movies/smoothslice_rhodepthtempgiz_"+run_id+"_"+output_dir+".mp4"
    #cmd = "ffmpeg -y -r 24 -i ../pics/sphplot"+run_id+output_dir+"%03d.png -c:v mpeg4 -q:v 1 /export/1/djw/movies/smoothslice_velrhogiz_"+run_id+"_"+output_dir+".mp4"

    #cmd = "ffmpeg -y -r 24 -i ../pics/sphplot"+run_id+output_dir+"%03d.png -c:v mpeg4 -q:v 1 /export/1/djw/movies/smoothsum_TTdustgiz_"+run_id+"_"+output_dir+".mp4"
    #cmd = "ffmpeg -y -r 24 -i ../pics/sphplot"+run_id+output_dir+"%03d.png -c:v mpeg4 -q:v 1 /export/1/djw/movies/smoothslice_TTdustgiz_"+run_id+"_"+output_dir+".mp4"

    #cmd = "ffmpeg -y -r 24 -i ../pics/sphplot"+run_id+output_dir+"%03d.png -c:v mpeg4 -q:v 1 /export/1/djw/movies/smoothsum_rhogiz_"+run_id+"_"+output_dir+".mp4"

    #cmd = "ffmpeg -y -r 24 -i ../pics/sphplot"+run_id+output_dir+"%03d.png -c:v mpeg4 -q:v 1 /export/1/djw/movies/smoothmin_dtgiz_"+run_id+"_"+output_dir+".mp4"
    
    #cmd = "ffmpeg -y -r 24 -i ../pics/sphplot"+run_id+output_dir+"%03d.png -c:v mpeg4 -q:v 1 /export/1/djw/movies/smoothsum_rhoTgiz_"+run_id+"_"+output_dir+".mp4"
    #cmd = "ffmpeg -y -r 24 -i ../pics/sphplot"+run_id+output_dir+"%03d.png -c:v mpeg4 -q:v 1 /export/1/djw/movies/smoothslice_rhoTgiz_"+run_id+"_"+output_dir+".mp4"

    #cmd = "ffmpeg -y -r 24 -i ../pics/sphplot"+run_id+output_dir+"%03d.png -c:v mpeg4 -q:v 1 /export/1/djw/movies/smoothslice_rhozoomgiz_"+run_id+"_"+output_dir+".mp4"
    #cmd = "ffmpeg -y -r 24 -i ../pics/sphplot"+run_id+output_dir+"%03d.png -c:v mpeg4 -q:v 1 /export/1/djw/movies/smoothslice_velrhogiz_"+run_id+"_"+output_dir+".mp4"

    #cmd = "ffmpeg -y -r 24 -i ../pics/sphplot"+run_id+output_dir+"%03d.png -c:v mpeg4 -q:v 1 /export/1/djw/movies/smoothslice_emitrhoTgiz_"+run_id+"_"+output_dir+".mp4"

    #cmd = "ffmpeg -y -r 12 -i ../pics/sphplot"+run_id+output_dir+"%03d.png -c:v mpeg4 -q:v 1 /export/1/djw/movies/smoothsum_rhovelgiz_"+run_id+"_"+output_dir+".mp4"
    #cmd = "ffmpeg -y -r 24 -i ../pics/sphplot"+run_id+output_dir+"%03d.png -c:v mpeg4 -q:v 1 /export/1/djw/movies/smoothsum_velgiz_"+run_id+"_"+output_dir+".mp4"

    print(cmd)
    os.system(cmd)
