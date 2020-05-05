print("Importing")

# Import libraries to do our magic
import numpy as np
import sys
import os
from joblib import Parallel, delayed
from multiprocessing import Pool

#from sys import path
#path.append("../src/")
import sph_frame

from sys import path
path.append("../")
import gizmo_tools

global Fixed_Binary
Fixed_Binary = True
print("Running")

run_id = sys.argv[1]
output_dir = sys.argv[2]
snapx = int(sys.argv[3])
nprocs = int(sys.argv[4])
angle_in = float(sys.argv[5])
outdir = sys.argv[6]

# fullDir = "/export/1/djw/gizmos/"+run_id+"/"+output_dir

gizmoDir = gizmo_tools.getGizmoDir(run_id)
movieDir = gizmo_tools.getMovieDir()
fullDir = gizmoDir+"/"+run_id+"/"+output_dir

#snapf = gizmo_tools.lastConsecutiveSnapshot(run_id,output_dir,False)
angle = angle_in/360.*2.*np.pi

# max run
#     if ( maxsnapf>-1 and snapf>maxsnapf ):

#infiles = [fullDir+"/snapshot_"+("%03d" % snapx)+".hdf5" for snapx in isnaps]
#nfiles = len(infiles)

infile = fullDir+"/snapshot_%03d.hdf5" % snapx


# rots = np.linspace(0.,np.pi/2.,24) # representative
#rots = np.linspace(0.,np.pi*2.,360.) # smooth spin
# rots = np.linspace(0.,np.pi*2.,64) # smooth-ish
#rots = np.linspace(0.,np.pi*2.,2) # test

#rots+=10./360.*np.pi*2.
rots=[[0.,angle],[0.3927,angle],[0.7854,angle],[1.1781,angle]]

#titlesuffixes = [r" $\phi=%2d^\circ$"%(theta*180./np.pi) for theta in rots]

#outfiles = ["/sfs/fs5/home-sh/supas356/djw_gizmo_tools/pics_test/sphrotplot"+run_id+output_dir+"%03d.png"%(snapx) for snapx in isnaps]
outfiles = "/sfs/fs5/home-sh/supas356/djw_gizmo_tools/"+outdir+"/sphrotplot"+run_id+output_dir+"%2f_%03d.png"%(angle_in,snapx)

#os.system("rm ../pics_test/sphrotplot"+run_id+output_dir+"???.png")

#Parallel(n_jobs=nprocs)(delayed(sph_frame.makesph_binary_view)(infiles[snap],outfiles[snap],cmap='plasma',plot=['view'],L=800,views=['side'],rot=rots,scale=5.,visibleAxes=False) for snap in isnaps)
#sph_frame.makesph_binary_view(infile,outfiles,cmap='plasma',plot=['view'],L=800,views=['side'],rot=rots,scale=5.,visibleAxes=False)
#sph_frame.makesph_binary_view(infile,outfiles,cmap='plasma',plot=['dens'],L=800,views=['side'],rot=rots,scale=1.,visibleAxes=False)
sph_frame.makesph_trhoz_frame(infile,outfiles,cmap='plasma',flat=False,ring=False,plot=['vels'],L=800,views=['side'],rot=[0.,90.],scale=1.,visibleAxes=False)

#def frame_bin(snap):
#  return sph_frame.makesph_binary_view(infiles[snap],outfiles[snap],cmap='plasma',plot=['view'],L=800,views=['side'],rot=rots,scale=5.,visibleAxes=False)
#with Pool(processes=nprocs) as pool:
#  pool.map(frame_bin,range(snapf+1-snapi))

#Parallel(n_jobs=nprocs)(delayed(sph_frame.makesph_binary_view)(infile,outfiles[0],cmap='plasma',flat=True,ring=False,plot=['view'],L=800,views=['side'],rot=[[0.,0.],[0.7854,0.],[1.5708,0.],[2.3561,0.]],scale=5.,visibleAxes=False))
#dies ist die urspruengliche variante: Parallel(n_jobs=nprocs)(delayed(sph_frame.makesph_trhoz_frame)(infile,outfiles[irot],cmap='plasma',flat=True,ring=False,plot=['view'],L=800,views=['side'],rot=[0.,phi],scale=100.,visibleAxes=False) for irot,phi in enumerate(rots))

#cmd = "ffmpeg -y -r 24 -i ../pics_spin/sphrotplot"+run_id+output_dir+"%03d_"%snapx+"%03d.png -c:v mpeg4 -q:v 1 ../../movies/rotateview_"+run_id+"_"+output_dir+"_%03d"%snapx+".mp4"

#os.system(cmd)
