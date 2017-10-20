print("Importing")

# Import libraries to do our magic
import numpy as np
import sys
import os
from joblib import Parallel, delayed

import sph_frame
import re

from sys import path
path.append("../")
import gizmodatadir

import argparse

def doop(*args):
    print(args)

def sort_nicely( l ):
    """ Sort the given list in the way that humans expect.
    """
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    l.sort( key=alphanum_key )

if __name__ == '__main__':
    default_values = dict()
    default_values["nprocs"]=8
    default_values["maxsnapf"]=-1
    default_values["rad"]=15.
    default_values["L"]=400
    default_values["plot"]="dens,temp"
    default_values["views"]="face,side"
    default_values["cmap"]="viridis"
    default_values["slice"]=False
    
    parsevals = ["nprocs","maxsnapf","run_id","output_dir","plot","cmap","rad","L","slice","views"]

    parser = argparse.ArgumentParser()
    parser.add_argument('run_id',help="name of superdirectory for runs")
    parser.add_argument('output_dir',help="name of subdirectory for run")
    parser.add_argument('--nprocs',type=int,help="processors to run on (default {})".format(default_values["nprocs"]))
    parser.add_argument('--maxsnapf',type=int,help="snapshot to end on (default=-1=do all snapshots)")
    parser.add_argument('--rad',type=float,help="radius of plot in parsecs")
    parser.add_argument('--L',type=int,help="size of plot area in pixels")
    parser.add_argument('--plot',type=str,help="values to plot, separated by commas")
    parser.add_argument('--views',type=str,help="face and/or side view, separated by commas")
    parser.add_argument('--cmap',type=str,help="colourmap palette")
    parser.add_argument('--slice',help="option to be slice plot",action='store_true')
    args = parser.parse_args()
    
#     run_id = args.run_id
#     output_dir = args.output_dir
    
    # TODO: don't do this!
    for parseval in parsevals:
        if ( vars(args)[parseval] ):
            vars()[parseval] = vars(args)[parseval]
#             print("setting {} to {}, current value:{}".format(parseval,vars(args)[parseval],vars()[parseval]))
        else:
            if ( parseval in default_values ):
                vars()[parseval] = default_values[parseval]
#                 print("setting {} to default {}, current value:{}".format(parseval,default_values[parseval],vars()[parseval]))
            else:
                raise Exception("No default value for {} - it must be specified!".format(parseval))

    flatPlot = not slice

    if ( flatPlot ):
        smooth_str = "smooth"
    else:
        smooth_str = "slice"
    
    toplot = plot.split(",")
    outp_plot = "".join(toplot)

    toview = views.split(",")
    outp_views = "".join(toview)

    print("Running")

    #run_id = sys.argv[1]
    #output_dir = sys.argv[2]
    
#     if ( len(sys.argv)>3 ):
#         nprocs = int(sys.argv[3])
#     else:
#         nprocs = default_procs
# 
#     if ( len(sys.argv)>4 ):
#         maxsnapf = int(sys.argv[4])
#     else:
#         maxsnapf = -1

    #snapi = int(sys.argv[3])
    #snapf = int(sys.argv[4])

    #determine which files to look at


    snapi = 0
    gizmoDir = gizmodatadir.gizmoDir()
    movieDir = gizmodatadir.movieDir()
    fullDir = gizmoDir+"/"+run_id+"/"+output_dir

    fnames = os.listdir(fullDir)
    sort_nicely(fnames)
    fnames = np.array(fnames)
    #fnames.sort()
    snapshotfilebools = np.array([x.startswith("snapshot") for x in fnames])
    snapshotfiles = fnames[snapshotfilebools]

    snapf = 0
    ctime = os.path.getmtime(fullDir+"/snapshot_000.hdf5")
    for fname in snapshotfiles[1:]:
        new_snapf = int(fname[9:len(fname)-5])
        new_ctime = os.path.getmtime(fullDir+"/"+fname)
        if ( new_ctime>ctime ) :
            ctime = new_ctime
            snapf = new_snapf
    
    
# max run
    if ( maxsnapf>-1 and snapf>maxsnapf ):
        print("Forcing snapf down from {} to {}".format(snapf,maxsnapf))
        snapf = maxsnapf

    snapf = maxsnapf

    os.system("rm ../pics/sphplot"+run_id+output_dir+"???.png")


    print("nfiles:",snapf-snapi+1)
    
    infiles = [fullDir+"/snapshot_"+("%03d" % snapx)+".hdf5" for snapx in range(snapi,snapf+1)]
    outfiles = ["../pics/sphplot"+run_id+output_dir+"%03d.png"%snapx for snapx in range(snapi,snapf+1)]

    #Parallel(n_jobs=nprocs)(delayed(sph_frame.makesph_trhoz_frame)(infiles[i],outfiles[i],cmap=cmap,flat=flatPlot,ring=flatPlot,plot=toplot,L=L,scale=rad) for i in range(snapi,snapf+1))
    Parallel(n_jobs=nprocs)(delayed(sph_frame.makesph_trhoz_frame)(infiles[i],outfiles[i],cmap=cmap,flat=flatPlot,ring=flatPlot,plot=toplot,L=L,scale=rad,views=toview) for i in range(snapi,snapf+1))


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
    cmd = "ffmpeg -y -r 24 -i ../pics/sphplot"+run_id+output_dir+"%03d.png -c:v mpeg4 -q:v 1 "+movieDir+"/"+smooth_str+"sum_"+outp_plot+"giz_"+run_id+"_"+output_dir+".mp4"


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
