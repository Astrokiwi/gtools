print("Importing")

# Import libraries to do our magic
import numpy as np
import sys
import os
#from joblib import Parallel, delayed

#from sys import path
#path.append("../src/")
import sph_frame

print("Running")

run_id = sys.argv[1]
output_dir = sys.argv[2]

snapx = int(sys.argv[3])

fullDir = "/export/1/djw/gizmos/"+run_id+"/"+output_dir

infile = fullDir+"/snapshot_"+("%03d" % snapx)+".hdf5"

#rots = np.linspace(0.,np.pi/2.,20)
rots = np.linspace(0.,np.pi*2.,80)

outfiles = ["../pics/sphrotplot"+run_id+output_dir+"%03d_%03d.png"%(snapx,irot) for irot in range(rots.size)]

os.system("rm ../pics/sphplot"+run_id+output_dir+"%03d"%snapx+"_???.png")


for irot,phi in enumerate(rots):
    #sph_frame.makesph_trhoz_frame(infile,outfile=outfiles[irot],cmap='plasma',flat=True,ring=True,plot=['dens'],L=600,cols=1,rot=[0.,phi])
    sph_frame.makesph_trhoz_frame(infile,outfile=outfiles[irot],cmap='plasma',flat=True,ring=False,plot=['view'],L=600,cols=1,rot=[0.,phi],scale=40.)

cmd = "ffmpeg -y -r 24 -i ../pics/sphrotplot"+run_id+output_dir+"%03d_"%snapx+"%03d.png -c:v mpeg4 -q:v 1 /export/1/djw/movies/rotateview_"+run_id+"_"+output_dir+"_%03d"%snapx+".mp4"
os.system(cmd)
