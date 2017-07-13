import numpy as np
import matplotlib as mpl
mpl.use('Agg')

import pylab as P
import itertools as it


def linearluminance(cmap):
    """Make the luminance linear"""
    cmap = P.cm.get_cmap(cmap)
    colors = cmap(np.arange(cmap.N))
    
    # convert RGBA to perceived greyscale luminance
    # cf. http://alienryderflex.com/hsp.html
    RGB_weight = [0.299, 0.587, 0.114]
    luminance = np.sqrt(np.dot(colors[:, :3] ** 2, RGB_weight))
    maxlum = np.max(luminance)
    for i in range(3):
       colors[:,i]/=luminance[:]
    
    minnewlum = 0.3
    newluminance = np.arange(cmap.N)+cmap.N*(minnewlum)/(1-minnewlum)
    for i in range(3):
       colors[:,i]*=newluminance
    
    maxcol = np.max(colors[:,:3])
    colors[:,:3]/=maxcol
    
    return cmap.from_list(cmap.name + "_linlum", colors, cmap.N)

def flatluminance(cmap):
    """Make the luminance flat"""
    cmap = P.cm.get_cmap(cmap)
    colors = cmap(np.arange(cmap.N))
    
    # convert RGBA to perceived greyscale luminance
    # cf. http://alienryderflex.com/hsp.html
    RGB_weight = [0.299, 0.587, 0.114]
    luminance = np.sqrt(np.dot(colors[:, :3] ** 2, RGB_weight))
    maxlum = np.max(luminance)
    colors[:,0] = colors[:,0]/luminance[:]
    colors[:,1] = colors[:,1]/luminance[:]
    colors[:,2] = colors[:,2]/luminance[:]
    
    maxcol = np.max(colors[:,:3])
    colors[:,:3]/=maxcol
    
    return cmap.from_list(cmap.name + "_delum", colors, cmap.N)


bins=200

x_basis = np.linspace(-2.,2.5,num=bins)
y_basis = np.linspace(0.,12.,num=bins)

#x, y = np.array([[x,y] for x,y in it.product(x_basis,y_basis)])

x,y = np.meshgrid(x_basis,y_basis)

#z=(5.9-(-2.7446+1.5439*xy[:,0]-.37046*xy[:,0]**2+.21642*xy[:,0]**3-.34755*xy[:,0]**4+.10114*xy[:,0]**5))-xy[:,1]
z = (5.9-(-2.7446+1.5439*x-.37046*x**2+.21642*x**3-.34755*x**4+.10114*x**5))-y

#z = x

this_cmap = P.get_cmap("jet_r")
this_cmap = flatluminance(this_cmap)

P.pcolormesh(x_basis+6.,y_basis,z,cmap=this_cmap,vmin=0.,vmax=9.)
P.colorbar(label=r"$\log_{10} t_\mathrm{sput}$ (yr)")
P.xlabel(r"$\log_{10} T$ (K)")
P.ylabel(r"$\log_{10} n_\mathrm{H}$ (cm$^{-3}$)")
P.suptitle("Sputtering time-scales (a/(da/dt)) from analytic fit of Tielens et al 1994 for silicates\n Assuming grain radius = 60e-4 cm")

P.savefig("pics/sputter_plot.png",dpi=300)

P.close()