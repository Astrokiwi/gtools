print("Importing")

# Import libraries to do our magic
import numpy as np
import sys
import os
from joblib import Parallel, delayed

#from sys import path
#path.append("../src/")
import sph_frame

from sys import path
path.append("../")
import gizmo_tools

print("Running")

run_id = sys.argv[1]
output_dir = sys.argv[2]

snapx = int(sys.argv[3])


if ( len(sys.argv)>4 ):
    nprocs = int(sys.argv[4])
else:
    nprocs = 80

# fullDir = "/export/1/djw/gizmos/"+run_id+"/"+output_dir

gizmoDir = gizmo_tools.getGizmoDir(run_id)
movieDir = gizmo_tools.getMovieDir()
fullDir = gizmoDir+"/"+run_id+"/"+output_dir


infile = fullDir+"/snapshot_%03d.hdf5" % snapx

rots = np.linspace(0.,np.pi/2.,24) # representative
# rots = np.linspace(0.,np.pi*2.,160) # smooth
#rots = np.linspace(0.,np.pi*2.,2) # test

titlesuffixes = [r" $\phi=%2d^\circ$"%(theta*180./np.pi) for theta in rots]
# print(titlesuffixes)
# sys.exit()

# rots+=80./360.*2.*np.pi # start at an angle

outfiles = ["../pics/sphrotplot"+run_id+output_dir+"%03d_%03d.png"%(snapx,irot) for irot in range(rots.size)]

os.system("rm ../pics/sphrotplot"+run_id+output_dir+"%03d"%snapx+"_???.png")



# for irot,phi in enumerate(rots):
#     sph_frame.makesph_trhoz_frame(infile,outfile=outfiles[irot],cmap='plasma',flat=True,ring=False,plot=['view'],L=600,cols=1,rot=[0.,phi],scale=40.)

# Parallel(n_jobs=nprocs)(delayed(sph_frame.makesph_trhoz_frame)(infile,outfiles[irot],cmap='plasma',flat=True,ring=False,plot=['view'],L=800,views=['face'],rot=[0.,phi],scale=3.,visibleAxes=False) for irot,phi in enumerate(rots))


# Parallel(n_jobs=nprocs)(delayed(sph_frame.makesph_trhoz_frame)(infile,outfiles[irot],cmap='plasma',flat=True,ring=False,plot=['view'],L=32,views=['face'],rot=[0.,phi],scale=6.,pixsize=64,titlesuffix=titlesuffixes[irot]) for irot,phi in enumerate(rots))
# pretty plot for talks etc
Parallel(n_jobs=nprocs)(delayed(sph_frame.makesph_trhoz_frame)(infile,outfiles[irot],cmap='plasma',flat=True,ring=False,plot=['view'],L=32,pixsize=16,views=['face'],rot=[0.,phi],scale=.4,visibleAxes=True) for irot,phi in enumerate(rots))
# Parallel(n_jobs=nprocs)(delayed(sph_frame.makesph_trhoz_frame)(infile,outfiles[irot],cmap='plasma',flat=True,ring=False,plot=['faceTemp'],L=800,views=['face'],rot=[0.,phi],scale=.2,visibleAxes=True) for irot,phi in enumerate(rots))
# Parallel(n_jobs=nprocs)(delayed(sph_frame.makesph_trhoz_frame)(infile,outfiles[irot],cmap='plasma',flat=True,ring=False,plot=['temp'],L=800,views=['face'],rot=[0.,phi],scale=1.6,visibleAxes=False) for irot,phi in enumerate(rots))
cmd = "ffmpeg -y -r 24 -i ../pics/sphrotplot"+run_id+output_dir+"%03d_"%snapx+"%03d.png -c:v mpeg4 -q:v 1 ../../movies/rotateview_"+run_id+"_"+output_dir+"_%03d"%snapx+".mp4"

# pixelated plot for paper
# Parallel(n_jobs=nprocs)(delayed(sph_frame.makesph_trhoz_frame)(infile,outfiles[irot],cmap='plasma',flat=True,ring=False,plot=['view'],L=32,views=['face'],rot=[0.,phi],scale=6.,pixsize=64,titlesuffix=titlesuffixes[irot]) for irot,phi in enumerate(rots))

# Parallel(n_jobs=nprocs)(delayed(sph_frame.makesph_trhoz_frame)(infile,outfiles[irot],cmap='plasma',flat=True,ring=False,plot=['facetemp'],L=32,views=['face'],rot=[0.,phi],scale=6.,pixsize=16) for irot,phi in enumerate(rots))
# cmd = "ffmpeg -y -r 24 -i ../pics/sphrotplot"+run_id+output_dir+"%03d_"%snapx+"%03d.png -c:v mpeg4 -q:v 1 ../../movies/facetemprot_"+run_id+"_"+output_dir+"_%03d"%snapx+".mp4"

# Parallel(n_jobs=nprocs)(delayed(sph_frame.makesph_trhoz_frame)(infile,outfiles[irot],cmap='plasma',flat=True,ring=False,plot=['tdust'],L=32,views=['face'],rot=[0.,phi],scale=4.,pixsize=16) for irot,phi in enumerate(rots))
# cmd = "ffmpeg -y -r 24 -i ../pics/sphrotplot"+run_id+output_dir+"%03d_"%snapx+"%03d.png -c:v mpeg4 -q:v 1 ../../movies/tdustrot_"+run_id+"_"+output_dir+"_%03d"%snapx+".mp4"

# Parallel(n_jobs=nprocs)(delayed(sph_frame.makesph_trhoz_frame)(infile,outfiles[irot],cmap='plasma',flat=True,ring=False,plot=['dens'],L=900,views=['side'],rot=[0.,phi],scale=5.,visibleAxes=False) for irot,phi in enumerate(rots))
# cmd = "ffmpeg -y -r 24 -i ../pics/sphrotplot"+run_id+output_dir+"%03d_"%snapx+"%03d.png -c:v mpeg4 -q:v 1 ../../movies/rotatesmooth_"+run_id+"_"+output_dir+"_%03d"%snapx+".mp4"

os.system(cmd)
