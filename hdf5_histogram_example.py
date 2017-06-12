# This produces a "heat map", giving the number of points in each (2D) bin
# This is written in Python, which is extremely useful in some circumstances where Fortran is a bit annoying

#####
# Parameters
# YOU WILL NEED TO CHANGE THESE
fileName = "/export/1/djw/gizmo_public/fortest_out/snapshot_000.hdf5" # the file to read in
dataSetKey = "/PartType0/Coordinates"
columns = (0,1) # List of the columns to plot. 0=1st column, 1=2nd column, etc
nbins = [1e3,1e3] # number of "bins" in each direciton
pic_name = "pics/test.png" # set the name of the file to whatever you want
x_cut = 17. # max extent for x axis
y_cut = 2. # max extent for y axis
xlabel = "X" # labels for axes
ylabel = "Y"
title = "test" # title of graph
logColour = False # True or False - should I log the colour-scale?
######

print("Importing")

# Import libraries to do our magic
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math
import h5py


# Read in the data - because Python is magical, this is all done in one line
print("Reading")
#x_in,y_in = np.loadtxt(fileName, usecols = columns, unpack=True)

f = h5py.File(fileName,"r")
d = f[dataSetKey]
x_in = d[:,columns[0]]
y_in = d[:,columns[1]]

#cut data
#data_ok = (x_in<x_cut) & (y_in<y_cut)
#x_in = x_in[data_ok]
#y_in = y_in[data_ok]


# bin the data - because Python is magical, this is all done in one line
print("Binning")
heatmap, yedges, xedges = np.histogram2d(y=x_in, x=y_in,bins=nbins)
myextent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
mywidth = myextent[1]-myextent[0]
myheight = myextent[2]-myextent[3]
myaspect = abs(mywidth/myheight)

# Plot the data! Okay, this takes a few lines, but only because we're setting up the axes etc
print("Plotting")
plt.clf() # clear the frame
if logColour:
    plt.imshow(heatmap,extent=myextent,norm=mpl.colors.LogNorm(), interpolation='none',aspect=myaspect) # plot the points with a grey colour scheme
else:
    plt.imshow(heatmap,extent=myextent, interpolation='none',aspect=myaspect) # plot the points with a grey colour scheme
plt.colorbar() # add a colour bar

plt.xlabel(xlabel) # Set the labels to whatever you want them to be
plt.ylabel(ylabel)
plt.title(title)

# save the figure to a png file
plt.savefig(pic_name,dpi=300) # dpi sets the size of the image. higher dpi = bigger image
