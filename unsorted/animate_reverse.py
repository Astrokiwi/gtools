import matplotlib as mpl
mpl.use('Agg')
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as P
from multiprocessing import Pool
import os

def make_frame(iframe):
    print("frame=",iframe)
    xyz = np.loadtxt("blob_reverse_frames/frame%06d.dat"%iframe,usecols=(0,1,2))
    P.close('all')
    fig = P.figure()
#     ax = fig.add_subplot(111)
#     R=np.sqrt(xyz[:,0]**2+xyz[:,1]**2)
#     ax.scatter(R,xyz[:,2])
#     ax.set_xlim(0.,1.e-1)
    ax = fig.add_subplot(111,projection='3d')
    size = 1.
    ax.set_xlim(-size,size)
    ax.set_ylim(-size,size)
    ax.set_zlim(-size,size)
    ax.scatter(xyz[:,0],xyz[:,1],xyz[:,2])
    P.savefig("pics/frame%06d.png"%iframe)

os.system("rm pics/frame*.png")
nprocs=128
pool = Pool(processes=nprocs)
pool.map(make_frame,range(991))

cmd = "ffmpeg -y -r 24 -i pics/frame%06d.png -c:v mpeg4 -q:v 1 ../movies/trajectory.mp4"
os.system(cmd)
