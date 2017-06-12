# This produces a "heat map", giving the number of points in each (2D) bin
# This is written in Python, which is extremely useful in some circumstances where Fortran is a bit annoying

#####
# Parameters
# YOU WILL NEED TO CHANGE THESE
fileName = "/export/1/djw/gizmo_public/disc_nocool_midres_out/snapshot_000.hdf5" # the file to read in
dataSetKey = "/PartType0/InternalEnergy"
pic_name = "pics/test.pdf" # set the name of the file to whatever you want
######

print("Importing")

# Import libraries to do our magic
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import pylab as P
import h5py
P.clf()

# Read in the data - because Python is magical, this is all done in one line
print("Reading")
#x_in,y_in = np.loadtxt(fileName, usecols = columns, unpack=True)

f = h5py.File(fileName,"r")
d = f["/PartType0/Coordinates"]
x_in = d[:,0]
y_in = d[:,1]

d = f[dataSetKey]
z_in = d[:]
z_in = z_in*100. # to K, approx

z_in = np.log10(z_in)

P.scatter(x_in,y_in,c=z_in,s=1.,marker=",",edgecolor='None')

# save the figure to a png file
P.savefig(pic_name,dpi=100) # dpi sets the size of the image. higher dpi = bigger image

P.close()