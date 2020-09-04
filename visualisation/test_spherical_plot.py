import matplotlib.pyplot as plt
import numpy as np
from sph_plotter import sph_plotter

N=10000

# x = np.random.random(N)-.5
# y = np.random.random(N)-.5
# z = np.random.random(N)-.5

theta = np.random.random(N)*2.*np.pi
z = np.random.random(N)*2.-1.
x = np.sqrt(1-z**2)*np.cos(theta)
y = np.sqrt(1-z**2)*np.sin(theta)

x+=2.

h = np.full(N,0.1)
m = np.ones(N)

corner = [0,-np.pi]
w = np.pi*2.
L = 128

y_edges = np.linspace(corner[0],corner[0]+w,L+1)
x_edges = np.linspace(corner[1],corner[1]+w,L+1)

g = sph_plotter.sph_dense_spherical(x,y,z,m,h,L,corner,w)

# rad_3d = np.sqrt(x**2+y**2+z**2)
# # rad_2d = np.sqrt(x**2+y**2)
# theta = np.arctan2(y,x) # azimuth
# phi = np.arccos(z/rad_3d) # inclination
# 

x_edges = np.degrees(x_edges)
y_edges = np.degrees(y_edges)

fig,sp = plt.subplots()
# sp.scatter(np.degrees(theta),np.degrees(phi))
sp.pcolormesh(x_edges,y_edges,g)
sp.set_xlabel("azimuth (°)")
sp.set_ylabel("inclination from pole (°)")
sp.set_ylim([0.,180.])
fig.savefig("../../figures/test.png")
plt.close(fig)