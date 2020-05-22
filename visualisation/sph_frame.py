# Import libraries to do our magic
import tblib.pickling_support
tblib.pickling_support.install()

import numpy as np
import math
import h5py
import matplotlib as mpl
mpl.use('Agg')

from matplotlib import colors

import pylab as P
#from joblib import Parallel, delayed
from scipy.ndimage.filters import gaussian_filter

from sph_plotter import sph_plotter

from sys import path, version, exit
path.append("../src/")
path.append("../")
import tab_interp

import sys
import traceback

import numpy as np

from difflib import SequenceMatcher

import itertools

import gizmo_tools
import pandas as pd

from scipy import interpolate

class ExceptionWrapper(object):

    def __init__(self, ee):
        self.ee = ee
        __,  __, self.tb = sys.exc_info()

    def re_raise(self):
        raise self.ee.with_traceback(self.tb)


# from skimage import exposure


#from enum import Enum

def similar(a, b):
    return SequenceMatcher(None, a, b).ratio()

def raiseprint(a):
    print("EXCEPTION RAISED")
    print(a)
    raise a

# --- Custom colormaps ---


colsteps = [.375,.75]
stepsize = .001
cdict1 = {'red':   ((0.0, .5, .5),
                   (colsteps[0], .8, .8),
                   (colsteps[0]+stepsize, 0.0, 0.0),
                   (colsteps[1], 0.0, 0.0),
                   (colsteps[1]+stepsize, 0.0, 0.0),
                   (1.0, 0.0, 0.0)),

         'green':   ((0.0, 0.0, 0.0),
                   (colsteps[0], 0.0, 0.0),
                   (colsteps[0]+stepsize, 0.0, 0.0),
                   (colsteps[1], 0.0, 0.0),
                   (colsteps[1]+stepsize, .4, .4),
                   (1.0, .7, .7)),

         'blue':   ((0.0, 0.0, 0.0),
                   (colsteps[0], 0.0, 0.0),
                   (colsteps[0]+stepsize, .5, .5),
                   (colsteps[1], .8, .8),
                   (colsteps[1]+stepsize, 0.0, 0.0),
                   (1.0, 0.0, 0.0)),
        }

P.register_cmap(name='RGBsteps', data=cdict1)

dust_opacity_function = None
def load_interpolate_opacity(opac_mu,
                        opac_file="../prams/simple_dust_opacity.txt"):
    global dust_opacity_function
    if dust_opacity_function is None:
        dust_opacity_table = np.loadtxt(opac_file)
        dust_opacity_function = interpolate.interp1d(dust_opacity_table[:,0],dust_opacity_table[:,1])
    opacity=dust_opacity_function(opac_mu)
#     print("Setting opacity to ",opacity," cm^2/g")
    opacity*=0.000208908219 # convert to pc**2/solar mass for consistency
    return opacity



lines = ["co1","co2","hcn1","hcn2","h2_1","h2_2","h2_3"]
line_wavelengths = [866.727,433.438,845.428,422.796,2.121,28.18,9.66]
line_opacities = {line:load_interpolate_opacity(mu) for line,mu in zip(lines,line_wavelengths)}

# count "IRdust" as a line
lines = ["IRdust"]+lines
line_opacities.update({"IRdust":0.1*0.000208908219})

# class SliceType(Enum):
#     densslice = 0
#     weightslice = 1
#     zdensslice = 2
#     zweightslice = 3
#     vec2dslice = 4
#     viewslice = 5

densslice = 0
weightslice = 1
zdensslice = 2
zweightslice = 3
vec2dslice_flat = 4
viewslice = 5
minslice = 6
zminslice = 7
vorinoislice = 8
zvorinoislice = 9
maxslice = 10
zmaxslice = 11
maxdotslice = 12
mindotslice = 13
sdevslice = 14
weightviewslice = 15
zvec2dslice = 16


molecular_mass = 4./(1.+3.*.76)
proton_mass_cgs = 1.6726e-24
gamma_minus_one = 5./3.-1.
boltzmann_cgs = 1.38066e-16


# chTab = None
# interpTabVec = None

debug_mode = False

np.seterr(all='ignore') # don't worry about bad logs etc - we want to propagate NaNs

def if_not_debug(f):
    def empty_function(*args,**kwargs):
        pass
    
    if debug_mode:
        return empty_function
    else:
        return f

@if_not_debug
def verboseprint(*args,**kwargs):
    print(*args,**kwargs)

# adapted from http://www.astrobetter.com/wiki/tiki-index.php?page=python_radial_profiles
def azimuthalNorm(image, center=None):
    """
    Normalize profile by max in each radial bin.

    image - The 2D image
    center - The [x,y] pixel coordinates used as the center. The default is 
             None, which then uses the center of the image (including 
             fracitonal pixels).
    
    """
    # Calculate the indices from the image
    y, x = np.indices(image.shape)

    if not center:
        center = np.array([(x.max()-x.min())/2.0, (x.max()-x.min())/2.0])

    r = np.hypot(x - center[0], y - center[1])
    
    rmax = np.max(r)
    nbins = int(np.ceil(rmax))
    
    ind_r = r.astype(int)
    
    bins_n = [np.max(image[ind_r==ii]) for ii in range(nbins)]

    for ii in range(nbins):
        image[ind_r==ii] = image[ind_r==ii]/bins_n[ii]
    
    return image


class GadgetData:
    pass

class PlotParameters( object ):
    pass

#wrapper for pcolormesh so it doesn't die if there's nothing finite to plot
def safe_pcolormesh(sp,*args,**qwargs):
# args[2] is the actual map
    if np.sum(np.isfinite(args[2]))==0:
        args = list(args)
        args[2][:,:] = qwargs["vmin"]
        args = tuple(args)
    return sp.pcolormesh(*args,**qwargs)


def makesph_plot(fig,sp,cbax,x_p,y_p,z_p,zslice,val_p,m_p,h_p,L,mask,corner,width,cblabel,clow,chigh,cmap_label,mode,
                 dolog=True,
                 planenorm=False,
                 circnorm=False,
                 plusminus=False,
                 diverging=False,
                 visibleAxes=True,
                 cbar2=None,
                 cmap2=None,
                 contour=False,
                 gaussian=None,
                 symLog=None,
                 cax_orientation='vertical',
                 centrecross=True):

    cmap_label2=cmap2
    cbax2=cbar2

    this_cmap = P.get_cmap(cmap_label)

    #print(fig,sp,cbax,x_p,y_p,z_p,zslice,val_p,m_p,h_p,L,mask,corner,width,cblabel,clow,chigh,cmap_label,mode,dolog,planenorm,circnorm,plusminus,diverging,visibleAxes,cbar2,cmap2)
    if ( cmap_label2 ):
        this_cmap2 = P.get_cmap(cmap_label2)
    
    if ( not plusminus and not diverging ):
        this_cmap.set_bad('black',1.)
        this_cmap.set_under('black',1.)

    n = m_p.size
    
#     if np.unique(m_p).size>1:
#         print("not unique ",np.unique(m_p))
# 
#     for a,label in zip((x_p,y_p,m_p,h_p),("x","y","m","h")):
#         if np.any(~np.isfinite(a)):
#             print("bad ",label,n,np.sum(~np.isfinite(h_p)))
# 

    
#     if val_p is None:
#         raiseprint("val_p is None!")
    
#     print("calculating nerrors")
#     nerrors_in = np.sum(~np.isfinite(val_p))
    
    outmap = None
    
    if ( mode==weightslice ):
        map = sph_plotter.sph_weight(x_p,y_p,m_p,h_p,val_p,L,corner,width,mask,n)
    elif (mode==densslice ):
        map = sph_plotter.sph_dense(x_p,y_p,m_p,h_p,L,corner,width,mask,n)
    elif (mode==zdensslice ):
        map = sph_plotter.sph_dense_slice(x_p,y_p,m_p,h_p,L,corner,width,z_p,zslice,mask,n)
    elif (mode==zweightslice ):
        map = sph_plotter.sph_weight_slice(x_p,y_p,m_p,h_p,val_p,L,corner,width,z_p,zslice,mask,n)
    elif (mode==vec2dslice_flat ):
        map1 = sph_plotter.sph_weight(x_p,y_p,m_p,h_p,val_p[0],L,corner,width,mask,n)
        map2 = sph_plotter.sph_weight(x_p,y_p,m_p,h_p,val_p[1],L,corner,width,mask,n)
        norm_map = np.sqrt(map1**2+map2**2)
        map1/=norm_map
        map2/=norm_map
    elif (mode==zvec2dslice ):
        map1 = sph_plotter.sph_weight_slice(x_p,y_p,m_p,h_p,val_p[0],L,corner,width,z_p,zslice,mask,n)
        map2 = sph_plotter.sph_weight_slice(x_p,y_p,m_p,h_p,val_p[1],L,corner,width,z_p,zslice,mask,n)
        norm_map = np.sqrt(map1**2+map2**2)
        map1/=norm_map
        map2/=norm_map
    elif (mode==viewslice ):
        zarg = np.argsort(z_p)
#         map = sph_plotter.sph_optical_depth_los(x_p,y_p,m_p,h_p,val_p[0],val_p[1],L,corner,width,z_p,zarg,mask,n)
        map = sph_plotter.sph_optical_depth_los_area(x_p,y_p,m_p,h_p,val_p[0],val_p[1],L,corner,width,z_p,zarg,mask,n)
        
        total_luminosity = np.sum(map) * (width*2*3.086e+18/L)**2 / 3.839e33  # check normalisation
        print("log10(L/Lsun) = {}, mean flux = {} erg/s/cm**2, fullwidth = {} pc".format(np.log10(total_luminosity),np.mean(map),width*2))
    elif (mode==weightviewslice ):
        zarg = np.argsort(z_p)
#         map = sph_plotter.sph_optical_depth_los(x_p,y_p,m_p,h_p,val_p[0],val_p[1],L,corner,width,z_p,zarg,mask,n)

#         map = sph_plotter.sph_optical_depth_los_weight(x_p,y_p,m_p,h_p,val_p[0],val_p[1],val_p[2],L,corner,width,z_p,zarg,mask,n)
#         map = sph_plotter.sph_optical_depth_los_weight(x_p,y_p,m_p,h_p,val_p[0],val_p[1],val_p[2],L,corner,width,z_p,zarg,mask,n)
#         threshold_flux = 1.e-5
        threshold_flux = 1.e-7
        print("HACKY: USING THRESHOLD FLUX OF ",threshold_flux)
        map = sph_plotter.sph_optical_depth_los_weight_thresh(x_p,y_p,m_p,h_p,val_p[0],val_p[1],val_p[2],L,corner,width,z_p,zarg,mask,threshold_flux,n)


    elif (mode==minslice ):
        map = sph_plotter.sph_min(x_p,y_p,h_p,val_p,L,corner,width,mask,n)
    elif (mode==zminslice ):
        map = sph_plotter.sph_minslice(x_p,y_p,h_p,val_p,L,corner,width,mask,n)
    elif mode==vorinoislice:
        map = sph_plotter.sph_vorinoi(x_p,y_p,m_p,h_p,val_p,L,corner,width,mask,n)
    elif mode==zvorinoislice:
        map = sph_plotter.sph_vorinoi_slice(x_p,y_p,m_p,h_p,val_p,L,corner,width,z_p,zslice,mask,n)
    elif mode==maxslice:
        map = sph_plotter.sph_max(x_p,y_p,m_p,h_p,val_p,L,corner,width,mask,n)
    elif mode==zmaxslice:
        map = sph_plotter.sph_max_slice(x_p,y_p,m_p,h_p,val_p,L,corner,width,z_p,zslice,mask,n)
    elif mode==maxdotslice or mode==mindotslice:
        map = sph_plotter.sph_dot(x_p,y_p,val_p,L,corner,width,mode-maxdotslice,mask,n)
    elif mode==sdevslice:
        map = sph_plotter.sph_sdev(x_p,y_p,m_p,h_p,val_p,L,corner,width,mask,n)
    
#     nerrors = map[0,0]
    
#     print(nerrors_in,nerrors)
    
    
    if ( mode==vec2dslice_flat or mode==zvec2dslice ):
        #map1/=1.e5
        #map2/=1.e5
        map1 = map1.T
        map2 = map2.T
        
        norm_map = norm_map.T
        
        step = width/L
        xmids = np.arange(corner[0]+step,corner[0]+width+step,step)
        ymids = np.arange(corner[1]+step,corner[1]+width+step,step)
        #sp.set_axis_bgcolor('black') # deprecated apparently?
        sp.set_facecolor('black')
        norm_map = np.log10(norm_map)
        #qv = sp.quiver(xmids,ymids,map1,map2,norm_map,headwidth=10.,pivot='mid',cmap=this_cmap,clim=[clow,chigh])
        xedges = np.arange(corner[0],corner[0]+width,width/L)
        yedges = np.arange(corner[1],corner[1]+width,width/L)

        qv = safe_pcolormesh(sp,xedges,yedges,norm_map,cmap=this_cmap,vmin=clow,vmax=chigh)
        sp.streamplot(xmids,ymids,map1,map2)
        if ( visibleAxes ):
            cb = fig.colorbar(qv,label=cblabel,cax=cbax,orientation=cax_orientation)
        else:
            cbax.set_axis_off()
#         sp.set_xbound(corner[0],corner[0]+width)
#         sp.set_ybound(corner[1],corner[1]+width)
#         sp.set_xlim(corner[0],corner[0]+width)
#         sp.set_ylim(corner[1],corner[1]+width)
#         sp.set_xmargin(0)
#         sp.set_ymargin(0)
        outmap = [norm_map,map1,map2]
    else:
#         xedges = np.arange(corner[0],corner[0]+width,width/L)
#         yedges = np.arange(corner[1],corner[1]+width,width/L)
        step = width/L
        xedges = np.arange(corner[0]+step/2.,corner[0]+width+step/2.,step)
        yedges = np.arange(corner[1]+step/2.,corner[1]+width+step/2.,step)
        if ( planenorm ):
            #iy0 = np.argmin(np.abs(yedges))
            #map[:,:]=(map[:,:].T/map[:,iy0]).T
            map[:,:]=(map[:,:].T/np.max(map[:,:],axis=1)).T
        if ( circnorm ):
            map = azimuthalNorm(map)

        if ( plusminus ):
        
            plusmap = map.copy()
            minusmap = map.copy()
            
            isplus = (map>0.)
            isminus = (map<0.)
            
            if ( np.sum(isplus)>0 ):
                plusmap[isminus]=0.
                mesh1=safe_pcolormesh(sp,xedges,yedges,plusmap.T,cmap=this_cmap,vmin=clow,vmax=chigh,norm = colors.LogNorm())
                if ( visibleAxes ):
                    cb = fig.colorbar(mesh1,label=cblabel,cax=cbax,orientation=cax_orientation)
                else:
                    cbax.set_axis_off()
            
            if ( np.sum(isminus)>0 ):
                minusmap[isplus]=0.
                minusmap = -minusmap
            
                mesh2=safe_pcolormesh(sp,xedges,yedges,minusmap.T,cmap=this_cmap2,vmin=clow,vmax=chigh,norm = colors.LogNorm())
                if ( visibleAxes ):
                    cb2 = fig.colorbar(mesh2,label=cblabel,cax=cbax2,orientation=cax_orientation)
                else:
                    cbax2.set_axis_off()
            outmap = [plusmap.T,minusmap.T]
        else:
            if not debug_mode:
                verboseprint(np.nanmin(map),np.nanmax(map))
                
            if gaussian is not None:
                if gaussian>0.:
                    gaussian_pix = int(np.floor(L*gaussian/width))
                    if gaussian_pix>=1:
                        map = gaussian_filter(map,gaussian)
            

            if ( dolog ):
                map = np.log10(map)
#             finiteIndices = np.isfinite(map)
#             if ( np.sum(finiteIndices)==0 ):
#                 print(map)
#                 raise Exception("No finite values in map")
#             print(np.min(map[finiteIndices]),np.max(map[finiteIndices]))
#             map = exposure.equalize_hist(map)

            if symLog is not None:
                mesh = safe_pcolormesh(sp,xedges,yedges,map.T,cmap=this_cmap,vmin=clow,vmax=chigh,norm=colors.SymLogNorm(linthresh=symLog,linscale=symLog,vmin=clow,vmax=chigh))
            else:
                mesh = safe_pcolormesh(sp,xedges,yedges,map.T,cmap=this_cmap,vmin=clow,vmax=chigh)
            if ( visibleAxes ):
                cb = fig.colorbar(mesh,label=cblabel,cax=cbax,orientation=cax_orientation)
            else:
                cbax.set_axis_off()

            if ( contour ):
                sp.contour(xedges,yedges,gaussian_filter(map.T,5),vmin=clow,vmax=chigh,colors='white',levels=8,linewidths=.5)

            outmap = map.T
    
    sp.set_xlim([corner[0],corner[0]+width])
    sp.set_ylim([corner[1],corner[1]+width])
    if visibleAxes and centrecross:
        sp.plot([0],[0],'+g',markersize=10.,markeredgewidth=1.)
    #mesh = sp.pcolormesh(xedges,yedges,map.T,cmap=this_cmap) 
    #sp.axis('equal')
#     print(outmap.shape)
    return outmap
    
#     return nerrors


def load_gadget(infile, plot_thing,
                        centredens=False):
                        
#                         ,
#                         opac_mu=None,
#                         opac_file="../prams/simple_dust_opacity.txt"):
#     global chTab
#     global interpTabVec

    data = GadgetData()

    f = h5py.File(infile,"r")
    
    header = f["/Header"]
    time = header.attrs.get("Time")
    time*= 0.9778e9 # to yr
    time/=1.e6 # to Myr

    BH_data=f["/BH_binary"]
    Binary_pos_1 = BH_data.attrs.get("Binary_pos_1") 
    Binary_pos_2 = BH_data.attrs.get("Binary_pos_2")
    if(isinstance(Binary_pos_1,np.ndarray) & isinstance( Binary_pos_2,np.ndarray)):
        data.binary_positions = [Binary_pos_1,Binary_pos_2]
    else:
        data.binary_positions = None

    data.xyz = np.array(f["/PartType0/Coordinates"]) # kpc

    data.m_p = np.array(f["/PartType0/Masses"]) # 10^10 msun
    data.m_p*=1e+10 # 10^10 solar masses to solar masses
    
    n = data.m_p.size

    data.h_p = np.array(f["/PartType0/SmoothingLength"]) # kpc
    
    need_to_load = list(plot_thing)
    if centredens:
        need_to_load.append("nH")
    if ( "view" in need_to_load or "dusttau" in need_to_load or "emit" in need_to_load):
        need_to_load.append("tdust")
        need_to_load.append("opac")
        need_to_load.append("dg")
    if ( "facetemp" in need_to_load ):
        need_to_load.append("tdust")
        need_to_load.append("opac")
    if ( "vlos" in need_to_load ):
        need_to_load.append("opac")
        need_to_load.append("dg")
    if ( "vlos" in need_to_load or "vmag" in need_to_load or "vthin" in need_to_load):
        need_to_load.append("vels")
    if ( "emit" in need_to_load ):
        need_to_load.append("tdust")
        need_to_load.append("dg")
    if ( "tdust" in need_to_load ):
        need_to_load.append("table")
    if ( "dust" in need_to_load ):
        need_to_load.append("dg")
    if any(x in need_to_load for x in ("IRdust","IRdustm")):
        need_to_load.append("dg")
    if ( "dg" in need_to_load ):
        need_to_load.append("table")
    for line in lines:
        if "dv"+line in need_to_load:
            need_to_load.append("v"+line)
        if "v"+line in need_to_load:
            need_to_load.append("view"+line)
            need_to_load.append("opac")
            need_to_load.append("dg")
        if "view"+line in need_to_load:
#             need_to_load.append("v"+line)
            need_to_load.append(line+"m")
            need_to_load.append("opac")
            need_to_load.append("dg")
        if "vels"+line in need_to_load:
            need_to_load.append(line+"m")
        if line+"m" in need_to_load:
            need_to_load.append(line)
        if line in need_to_load and not "table" in need_to_load:
            need_to_load.append("table")
    if ( "col" in need_to_load ):
        if not "/PartType0/AGNColDens" in f:
            need_to_load.append("table")
    if ( "table" in need_to_load ):
        need_to_load.append("temp")
        need_to_load.append("nH")
        need_to_load.append("tau")
    if "vel0" in need_to_load or "age" in need_to_load:
        need_to_load.append("id")

#     if "nopac" in need_to_load:ga
#         need_to_load.append("opac")
#         need_to_load.append("nH")
    
    if any(x in need_to_load for x in ["vmag","vel_2d","vel_r","vel_x","vel_y","vel_z","vel_a"]+["v"+line for line in lines]+["vels"+line for line in lines]):
        need_to_load.append("vels")



    if "id" in need_to_load:
        data.id_p = np.array(f["/PartType0/ParticleIDs"]).astype(int)
#         data.id_p-=1

    if "age" in need_to_load:
        infile_split = infile.split("/")
        run_id = infile_split[-3]
        run_name = infile_split[-2]
        run_snapfile = infile_split[-1]
        run_isnap = int(run_snapfile[-8:-5])
        age_file = "../data/age_{}_{}_{}.dat".format(run_id,run_name,run_isnap)
        data.age = time-np.loadtxt(age_file)

    if ( "pres" in need_to_load ):
        data.pres = np.array(f["/PartType0/Pressure"])
#         data.pres*=1.989e+33 # from internal units to dyne/cm*8*2

    if ( "arad" in need_to_load ):
        data.arads = np.array(f["/PartType0/RadiativeAcceleration"])
        data.arads*=3.24086617e-12 # to cm/s/s
        data.arad = np.sqrt(np.sum(data.arads**2,axis=1))

    if ( "accel" in need_to_load ):
        data.accels = np.array(f["/PartType0/Acceleration"])
        data.accels*=3.24086617e-12 # to cm/s/s
        data.accel = np.sqrt(np.sum(data.accels**2,axis=1))


    if ( "depth" in need_to_load ):
        data.depth = np.array(f["/PartType0/AGNOpticalDepth"]) # Msun/kpc**2

    if ( "list" in need_to_load ):
        data.list = np.arange(n)

    if ( "rand" in need_to_load ):
        data.rand = np.random.random(n)
        #data.rand = np.random.randint(n,size=n)

    if ( "nneigh" in need_to_load ):
        data.nneigh = np.array(f["/PartType0/TrueNumberOfNeighbours"])
    
    if ( "temp" in need_to_load ):
        data.u_p = np.array(f["/PartType0/InternalEnergy"]) # 1e10 erg/g
        data.u_p*=1.e10 # to erg/g

        data.TK_p = (gamma_minus_one/boltzmann_cgs*(molecular_mass*proton_mass_cgs)*data.u_p)

    if ( "vels" in need_to_load ):
        data.vels = np.array(f["/PartType0/Velocities"]) # in km/s

    # doesn't work well
    if ( "dt" in need_to_load ):
        data.dt_p = np.array(f["/PartType0/TimeStep"])
        data.dt_p*=0.9778e9 # to yr

    # doesn't work well
    if ( "heat" in need_to_load ):
        data.heat = np.array(f["/PartType0/AGNHeat"])
        data.heat*=1e10/3.08568e+16# to erg/s/g
        
#         heatplus = (data.heat>0.)
#         heatminus = (data.heat<0.)
#         
#         data.heat[heatplus] = np.maximum(0.,np.log10(data.heat[heatplus]))
#         data.heat[heatminus] = -np.maximum(0.,np.log10(-data.heat[heatminus]))

#     if ( "depth" in need_to_load):
#         #rho_p = np.array(f["/PartType0/Density"])
#         #agn_heat_p = np.array(f["/PartType0/AGNHeat"])
#         data.depth_p = np.array(f["/PartType0/AGNDepth"]) # unitless

    if ( "nH" in need_to_load ):
        data.rho_p = np.array(f["/PartType0/Density"])
        data.rho_p*=6.77e-22 # to g/cm**3 
        data.nH_p = data.rho_p/(molecular_mass*proton_mass_cgs)

    if ( "tau" in need_to_load ):
        data.tau = np.array(f["/PartType0/AGNDepth"])

    if ( "AGNI" in need_to_load ):
        data.AGNI = np.array(f["/PartType0/AGNIntensity"])
        data.AGNI*= (1.989e53/(0.9778e9*3.154e7)/3.086e21**2) # convert from internal units (energy/Gyr-ish/kpc**2) to erg/s/cm**2

    if ( "table" in need_to_load ):
#         if ( not chTab ):
        verboseprint("Load dust tables")
#             chTab = tab_interp.CoolHeatTab( ("../coolheat_tab_marta/shrunk_table_labels_090818tau.dat"),
#                                             ("../coolheat_tab_marta/shrunk_table_090818_m0.0001_hsmooth_tau.dat"),
#                                             ("../coolheat_tab_marta/shrunk_table_labels_090818taunodust.dat"),
#                                             ("../coolheat_tab_marta/shrunk_table_090818_m0.0001_hsmooth_taunodust.dat"),
#                                             ("../coolheat_tab_marta/shrunk_table_labels_090818taudense.dat"),
#                                             ("../coolheat_tab_marta/shrunk_table_090818_m0.0001_hsmooth_taudense.dat")
#                                             )

      #      tableDate="281118"
#             tableDate="181018"
       #     tableRes="0.001"
       #     tableDate="010519"
#             tableRes="0.0001"
#             tableDate="281118"
#             tableRes="0.0001"
    #        tableDate="060319"
    #        tableRes="0.1"
#             tableDate="281118"
#             tableRes="0.0001"
        tableDate="060319"
        tableRes="0.1"
        cloudy_table = gizmo_tools.cloudy_table(tableDate,tableRes,"../")
#         cloudy_table = gizmo_tools.cloudy_table(tableDate,tableRes,"/sfs/fs2/work-sh1/supas356/tables/")
#             chTab = tab_interp.CoolHeatTab( ("../coolheat_tab_marta/shrunk_table_labels_"+tableDate+"tau.dat"),
#                                             ("../coolheat_tab_marta/shrunk_table_"+tableDate+"_m"+tableRes+"_hsmooth_tau.dat"),
#                                             ("../coolheat_tab_marta/shrunk_table_labels_"+tableDate+"taunodust.dat"),
#                                             ("../coolheat_tab_marta/shrunk_table_"+tableDate+"_m"+tableRes+"_hsmooth_taunodust.dat"),
#                                             ("../coolheat_tab_marta/shrunk_table_labels_"+tableDate+"taudense.dat"),
#                                             ("../coolheat_tab_marta/shrunk_table_"+tableDate+"_m"+tableRes+"_hsmooth_taudense.dat")
#                                             )
#             interpTabVec = np.vectorize(chTab.interpTab)
        data.flux_p = np.array(f["/PartType0/AGNIntensity"]) # energy per surface area per time
        data.flux_p*=1.989e+53/(3.086e21)**2/(3.08568e+16)
        
        table_particles = pd.DataFrame()
        table_particles["nH"] = data.nH_p
        table_particles["temp"] = data.TK_p
        table_particles["AGNIntensity"] = data.flux_p
        table_particles["AGNDepth"] = data.tau

        verboseprint("Calculating dust/cooling/heating properties from table")
#         tabStructs = interpTabVec(data.nH_p.astype(np.float64),data.TK_p.astype(np.float64),data.flux_p.astype(np.float64),data.tau.astype(np.float64))
        cloudy_table.interp(table_particles)
    
    if ( "col" in need_to_load ):
        if "/PartType0/AGNColDens" in f:
            data.coldens = np.array(f["/PartType0/AGNColDens"]) # Msun/kpc**2
            data.coldens*=(1.989e+43/3.086e+21**2) # to g/cm**2
            data.coldens/=(molecular_mass*proton_mass_cgs) # N in cm**(-2)
        else:
            data.coldens = table_particles["column_out"]
#             data.coldens = np.array(list(map(lambda y: y.column_out, tabStructs)))

    if ( "tdust" in need_to_load ):
        data.dustTemp = table_particles["dustT"]
#         data.dustTemp = map(lambda y: y.dustT, tabStructs)
#         data.dustTemp = np.array(list(data.dustTemp))
#         data.dustTemp = 10.**data.dustTemp
#         print(np.max(data.dustTemp),np.min(data.dustTemp))
#         data.dustTemp = data.dustTemp**4 # for test
    for line in lines[1:]:
        if line in need_to_load:
#             data.__dict__[line] = np.array(list(map(lambda y: y.__getattr__("line_"+line), tabStructs)))
            data.__dict__[line] = table_particles["line_"+line]
        if line+"m" in need_to_load:
            data.__dict__[line+"m"] = table_particles["line_"+line]*data.m_p*1.9891e33/9.52140614e36/(4.*np.pi) # erg/s/g to erg/s, extra factor for pc**2 to ster cm**2, output is erg/s/cm**2/ster
#     if ( "co1" in need_to_load ):
#         data.co1 = np.array(list(map(lambda y: y.line_co1, tabStructs))) # erg/s/g
# #         data.co1/=data.rho_p # to erg/g
#     if ( "co2" in need_to_load ):
#         data.co2 = np.array(list(map(lambda y: y.line_co2, tabStructs))) # erg/s/g
# #         data.co2/=data.rho_p # to erg/g
#     if ( "hcn1" in need_to_load ):
#         data.hcn1 = np.array(list(map(lambda y: y.line_hcn1, tabStructs))) # erg/s/g
# #         data.hcn1/=data.rho_p # to erg/g
#     if ( "hcn2" in need_to_load ):
#         data.hcn2 = np.array(list(map(lambda y: y.line_hcn2, tabStructs))) # erg/s/g
# #         data.hcn2/=data.rho_p # to erg/g
#     if ( "co1m" in need_to_load ):
#         data.co1m=data.co1*data.m_p*1.9891e33/9.52140614e36/(4.*np.pi) # erg/s/g to erg/s, extra factor for pc**2 to ster cm**2, output is erg/s/cm**2/ster
#     if ( "co2m" in need_to_load ):
#         data.co2m=data.co2*data.m_p*1.9891e33/9.52140614e36/(4.*np.pi) # erg/s/g to erg/s
#     if ( "hcn1m" in need_to_load ):
#         data.hcn1m=data.hcn1*data.m_p*1.9891e33/9.52140614e36/(4.*np.pi) # erg/s/g to erg/s
#     if ( "hcn2m" in need_to_load ):
#         data.hcn2m=data.hcn2*data.m_p*1.9891e33/9.52140614e36/(4.*np.pi) # erg/s/g to erg/s

    if ( "dg" in need_to_load ):
#         data.dg = map(lambda y: y.dg, tabStructs)
#         data.dg = np.array(list(data.dg))
        data.dg = table_particles["dg"]

    if ( "dust" in need_to_load ):
        data.dust = data.dg * data.m_p
    
    if ( "emit" in need_to_load ):
#         data.emissivity = 1.e-20*data.m_p * data.dustTemp**4. # arbitrary units
        data.emissivity = 5.67e-5 * data.dustTemp**4. * data.dg/np.nanmax(data.dg) # erg/s/cm^2

    if ( "opac" in need_to_load ):
        data.opac = np.array(f["/PartType0/AGNOpacity"]) # internal units: kpc**2/1e10 solar mass
#         data.opac*= 1.e-4 # to pc**2/solar mass
        data.opac*= 0.478679108 # to cm**2/g
#         print("faking opacity")
#         opacity = 65.2 # cm^2/g, somewhat arbitrary
#         opacity*=0.000208908219 # convert to pc**2/solar mass for consistency
#         data.opac = np.full((n),opacity)

    if ( "view" in need_to_load or "vlos" in need_to_load or any(x in need_to_load for x in ["view"+line for line in lines])
            or "dusttau" in need_to_load ):
        if ( "view" in need_to_load ):
            data.brightness = 5.67e-5 * data.dustTemp**4. * data.dg/np.nanmax(data.dg) # erg/s/cm^2
        #opacity = 652. # cm^2/g, somewhat arbitrary
        #opacity = 1. # cm^2/g, somewhat arbitrary

        # we want infrared opacity, not the flux-weighted opacity
#         opacity = 1. # cm^2/g, somewhat arbitrary
        
#         print("SUPER faking opacity")
#         opacity*=1.e3 # for lulz

#         if opac_mu and opac_file:
#             opacity = load_interpolate_opacity(opac_mu,opac_file)
# #             opac_data = np.loadtxt(opac_file)
# #             opacity=interpolate.interp1d(opac_data[:,0],opac_data[:,1])(opac_mu)
# #             print("Setting opacity to ",opacity," cm^2/g")
        print("faking opacity")
        opacity = 65.2 # cm^2/g, somewhat arbitrary
        print("Broad opacity is now ",opacity," cm^2/g")

        opacity*=0.000208908219 # convert to pc**2/solar mass for consistency
        data.opac = np.full((n),opacity)
        data.opac*=data.dg/np.nanmax(data.dg) # take into account dust fraction

        if "dusttau" in need_to_load:
            data.dusttau = data.opac * data.m_p
        for line in lines[1:]:
            if not "view"+line in need_to_load:
                continue
            data.__dict__[line+"opac"]=np.full((n),line_opacities[line])
            data.__dict__[line+"brightness"]=data.__dict__[line+"m"]/data.__dict__[line+"opac"] # erg/s/cm^2 - multiply by opacity to get actual emission
#             data.__dict__[line+"brightness"]=data.__dict__[line+"m"] # erg/s, gets SPH smoothed to get erg/s/cm**2

        
        if ( "view" in need_to_load ):
            #sputtered = (data.dustTemp>2.5e3) # or something, super arbitrary
            sputtered = (data.dustTemp>1.e5) # or something, super arbitrary
            data.brightness[sputtered] = 0.
            data.opac[sputtered] = 0.
    if any(x in need_to_load for x in ("IRdustm","IRdust","IRdustopac","IRdustbrightness","viewIRdust")):
        data.dustTemp = table_particles["dustT"]
        data.IRdustbrightness = 5.67e-5 * data.dustTemp**4. * data.dg/np.nanmax(data.dg)
        data.IRdustopac = np.full((n),line_opacities["IRdust"])
        data.IRdustm = data.IRdustbrightness*data.IRdustopac
        data.IRdust = data.IRdustm/(data.m_p*1.9891e33/9.52140614e36/(4.*np.pi))


    #rho_p*=6.77e-22 # to g/cm**3
    #m_p*=1.989e+43 # 10^10 solar masses to g
    #agn_heat_p*=1e10/3.08568e+16# to erg/s/g
    #agn_heat_p*=rho_p# to erg/cm**3/s
    
#     if ( 'agn_heat_p' in globals() ):
#         depth_p = agn_heat_p*rho_p*h_p*(16./3./np.pi)
    if "rad0" in need_to_load:
#         data.id_p = np.array(f["/PartType0/ParticleIDs"]).astype(int)
#         data.id_p-=1
        data.rad0 = np.load("rad0.npy")
        data.rad0 = data.rad0[data.id_p-1]
    
    return time,data

def pack_dicts(flatPlot=True,planenorm=False,cmap="viridis"):
    if ( flatPlot ):
        quantslice = weightslice
        dslice = densslice
        rhounit = r"$\log_{10} \Sigma$ (M$_\odot$/pc$^2$)"
        drange = [-2.,8.]
        thisminslice = minslice
        vslice = vorinoislice
        mslice = maxslice
        vec2dslice = vec2dslice_flat
#         drange = [-2.,2.] # TEMPORARY
        #drange = [2.,6.]
        #drange = [-1.,4.]
    else:
        quantslice = zweightslice
        dslice = zdensslice
        thisminslice = zminslice
        rhounit = r"$\log_{10} \rho$ (M$_\odot$/pc$^3$)"
        drange = [-2.,7.]
        vslice = zvorinoislice   
        mslice = zmaxslice
        vec2dslice = zvec2dslice
    if ( planenorm ):
        drange_s = [-3.6,0.]
        drange = [-4.5,0.]
    else:
        drange_s = drange


    # Pack the dictionaries!
    # We could do this with a bunch of if/endif statements, but
    # data-driven is better
    plotLabel = dict()
    plotLabel["temp"] = r"$\log_{10} T_g$ (K)"
    plotLabel["col"] = r"$\log_{10} N$ (cm$^{-2}$)"
    plotLabel["nH"] = r"$\log_{10} n_{H}$ (cm$^{-3}$)"
    plotLabel["heat"] = r"$H$"
    plotLabel["dens"] = rhounit
    plotLabel["dusttau"] = r"$\log_{10}\tau_d$"
    plotLabel["depth"] = r"$\tau$"
    plotLabel["vels"] = r"$\log_{10}v$ (km/s)"
    plotLabel["tdust"] = r"$\log_{10} T_d$ (K)"
    plotLabel["dg"] = r"$f_d$"
    plotLabel["dust"] = rhounit
    plotLabel["view"] = r"$\log_{10} F$ (arbitrary units)"
    plotLabel["vlos"] = r"$v_{LOS}$"
    plotLabel["vthin"] = r"$v_{LOS}$"
    plotLabel["emit"] = r"$\log_{10} F$"
    plotLabel["dt"] = r"$\log_{10}$ dt"
    plotLabel["vel_2d"] = r"$v_{2D}$"
    plotLabel["vel_x"] = r"$v_{x}$"
    plotLabel["vel_y"] = r"$v_{y}$"
    plotLabel["vel_z"] = r"$v_{z}$"
    plotLabel["vel_r"] = r"$v_{r}$ (km/s)"
    plotLabel["vmag"] = r"$|v|_{max}$ (km/s)"
    plotLabel["vel_a"] = r"$v_{\theta}$"
    plotLabel["arad"] = r"$a_{rad}$ (log cm/s/s)"
    plotLabel["accel"] = r"$a$ (log cm/s/s)"
    plotLabel["AGNI"] = r"Unextinced $I_{AGN}$ (log erg/s/cm^2)"
#     plotLabel["tau"] = r"$\log_{10}\tau$"
    plotLabel["tau"] = r"$\tau$"
    plotLabel["list"] = r"$i$"
    plotLabel["rand"] = r"$q$"
    plotLabel["opac"] = r"$\kappa$"
    plotLabel["smooth"] = r"$h_p$"
    plotLabel["facetemp"] = r"$\log_{10}T_d$ (K)"
    plotLabel["rad0"] = r"$R_0$ (pc)"
    plotLabel["nneigh"] = r"$N_n$"
#     plotLabel["pres"] = r"$P$ (dyne/cm$^2$)"
    plotLabel["pres"] = r"$P$ (internal)"
    plotLabel["age"] = r"$age$ (Myr)"
    for line in lines:
        plotLabel[line]=line
        plotLabel[line+"m"]=line+" line intensity - erg/s/cm**2/ster"
        plotLabel["v"+line]="v"+line
        plotLabel["dv"+line]="dv"+line
        plotLabel["view"+line]="view"+line
        plotLabel["vels"+line]=r"$\log_{10}v_{"+line+r"}$"
#     plotLabel["co1"] = r"co1"
#     plotLabel["co2"] = r"co2"
#     plotLabel["hcn1"] = r"hcn1"
#     plotLabel["hcn2"] = r"hcn2"
#     plotLabel["co1m"] = r"CO '1' line intensity - erg/s/cm**2/ster"
#     plotLabel["co2m"] = r"CO '2' line intensity - erg/s/cm**2/ster"
#     plotLabel["hcn1m"] = r"HCN '1' line intensity - erg/s/cm**2/ster"
#     plotLabel["hcn2m"] = r"HCN '2' line intensity - erg/s/cm**2/ster"
#     plotLabel["vco1"] = r"vco1"
#     plotLabel["vco2"] = r"vco2"
#     plotLabel["vhcn1"] = r"vhcn1"
#     plotLabel["vhcn2"] = r"vhcn2"
#     plotLabel["dvco1"] = r"dvco1"
#     plotLabel["dvco2"] = r"dvco2"
#     plotLabel["dvhcn1"] = r"dvhcn1"
#     plotLabel["dvhcn2"] = r"dvhcn2"

    plotRanges = dict()
#     plotRanges["temp"] = [.999,4.]*2
    plotRanges["temp"] = [.999,6.]*2
    #plotRanges["temp"] = [3.,4.]*2
    plotRanges["col"] = [18.,26.]*2
#     plotRanges["nH"] = [0.,8.]*2 # standard
    plotRanges["nH"] = [5.,10.]*2
#     plotRanges["nH"] = [5.,7.]*2
#     plotRanges["nH"] = [-1.,7.]*2
#     plotRanges["nH"] = [-7.,7.]*2
#     plotRanges["heat"] = [1e-2,1.e10]*2
#     plotRanges["heat"] = [1e-2,1.e10]*2
    plotRanges["heat"] = [-1e2,1.e2]*2
    plotRanges["dens"] = drange+drange_s
#     plotRanges["dusttau"] = [0.,10.]*2
    plotRanges["dusttau"] = [-3.,3.]*2
    plotRanges["depth"] = [0.,1.e3]*2
    plotRanges["vels"] = [0.,3.]*2
    plotRanges["tdust"] = [0.,3.]*2
#     plotRanges["tdust"] = [1.,3.]*2
    #plotRanges["tdust"] = [1.3,2.2]*2
#     plotRanges["tdust"] = [4.,8.]*2
#     plotRanges["dg"] = [0.00635,.006363]*2
    plotRanges["dg"] = [0.0,.000234]*2
    plotRanges["dust"] = drange+drange_s
#     plotRanges["view"] = [0.,4.]*2
    plotRanges["view"] = [1.,7.]*2
    plotRanges["vlos"] = [-1.2e2,1.2e2]*2
#     plotRanges["vthin"] = [-1.2e2,1.2e2]*2
    plotRanges["vthin"] = [-300.,300.]*2
    plotRanges["emit"] = [-1.,6.]*2
    plotRanges["dt"] = [-3.,3.]*2
    #plotRanges["vel_2d"] = [-25.,25.]*2
    plotRanges["vel_2d"] = [-250.,250.]*2
#     plotRanges["vel_r"] = [-25.,25.]*2
#     plotRanges["vel_r"] = [-3.e5,3.e5]*2
#     plotRanges["vel_r"] = [1.,4.]*2
#     plotRanges["vel_r"] = [-1.e2,1.e2]*2
    plotRanges["vel_r"] = [-25.,25.]*2
    plotRanges["vmag"] = [0.,120.]*2
    plotRanges["vel_x"] = [-200.,200.]*2
    plotRanges["vel_y"] = [-200.,200.]*2
#     plotRanges["vel_z"] = [-25.,25.]*2
    plotRanges["vel_z"] = [-2000.,2000.]*2
#     plotRanges["vel_a"] = [-100.,100.]*2
    plotRanges["vel_a"] = [-50.,50.]*2
    plotRanges["arad"] = [-9.,0.]*2
    plotRanges["accel"] = [-3.,0.]*2
#     plotRanges["AGNI"] = [4.,6.]*2
    plotRanges["AGNI"] = [-1.,6.]*2
#     plotRanges["tau"] = [0.,5.]*2
    plotRanges["tau"] = [0.,7.]*2
#     plotRanges["tau"] = [0.,1.]*2
#     plotRanges["tau"] = [-2,3.]*2
#     plotRanges["tau"] = [0.,250.]*2
#     plotRanges["tau"] = [0.,50.]*2
    plotRanges["list"] = [0.,1.e6]*2
    plotRanges["rand"] = [0.,1.]*2
#     plotRanges["opac"] = [-7.,-1.6]*2
    plotRanges["opac"] = [-3.,2]*2
#     plotRanges["opac"] = [-2.,-1.]*2
    plotRanges["smooth"] = [-4.,2.]*2
    plotRanges["facetemp"] = [0.,5.]*2
    plotRanges["rad0"] = [0.,4.]*2
    plotRanges["nneigh"] = [0,50]*2
#     plotRanges["pres"] = [36.3,38.6]*2
    plotRanges["pres"] = [3.3,16]*2
    plotRanges["age"] = [0.,3.]*2
    for line in lines:
#         plotRanges[line]=[-6.,-2.]*2
        plotRanges[line]=[-8.,2.]*2
#         plotRanges[line+"m"]=[-7.,-3.]*2
        plotRanges[line+"m"]=[-12.,0.]*2
#         plotRanges["v"+line]=[-100.,100.]*2
        plotRanges["v"+line]=[-200.,200.]*2
        plotRanges["dv"+line]=[0.,3.,0.,30.]
#         plotRanges["view"+line]=[-1.2e2,1.2e2]*2
#         plotRanges["view"+line]=[-12.,0.]*2
        plotRanges["view"+line]=[-20.,0.]*2
        plotRanges["vels"+line]=[0.,3.]*2
#         plotRanges["view"+line]=[-1.2e1,1.2e1]*2
#     plotRanges["co1"] = [-6.,-2.]*2
#     plotRanges["co2"] = [-6.,-2.]*2
#     plotRanges["hcn1"] = [-6.,-2.]*2
#     plotRanges["hcn2"] = [-6.,-2.]*2
#     plotRanges["co1m"] = [-7.,-3.]*2
#     plotRanges["co2m"] = [-7.,-3.]*2
#     plotRanges["hcn1m"] = [-7.,-3.]*2
#     plotRanges["hcn2m"] = [-7.,-3.]*2
#     plotRanges["vco1"] = [-100.,100.]*2
#     plotRanges["vco2"] = [-100.,100.]*2
#     plotRanges["vhcn1"] = [-100.,100.]*2
#     plotRanges["vhcn2"] = [-100.,100.]*2
#     plotRanges["dvco1"] = [0.,5.,0.,30.]
# #     plotRanges["dvco1"] = [-2.,2.]*2
#     plotRanges["dvco2"] = plotRanges["dvco1"]
#     plotRanges["dvhcn1"] = plotRanges["dvco1"]
#     plotRanges["dvhcn2"] = plotRanges["dvco1"]

    plotSliceTypes = dict()
    plotSliceTypes["temp"] = quantslice
    plotSliceTypes["col"] = quantslice
    plotSliceTypes["nH"] = quantslice
    plotSliceTypes["heat"] = quantslice
    plotSliceTypes["dens"] = dslice
    plotSliceTypes["dusttau"] = densslice
    plotSliceTypes["depth"] = zweightslice
    plotSliceTypes["vels"] = vec2dslice
    plotSliceTypes["tdust"] = quantslice
    plotSliceTypes["dg"] = quantslice
    plotSliceTypes["dust"] = dslice
    plotSliceTypes["view"] = viewslice
    plotSliceTypes["vlos"] = viewslice
    plotSliceTypes["emit"] = dslice
    plotSliceTypes["dt"] = thisminslice
    plotSliceTypes["vel_2d"] = quantslice
    plotSliceTypes["vel_r"] = quantslice
    plotSliceTypes["vmag"] = mslice # max
    plotSliceTypes["vel_x"] = quantslice
    plotSliceTypes["vel_y"] = quantslice
    plotSliceTypes["vel_z"] = quantslice
    plotSliceTypes["vthin"] = quantslice
    plotSliceTypes["vel_a"] = quantslice
    plotSliceTypes["arad"] = quantslice
    plotSliceTypes["accel"] = quantslice
    plotSliceTypes["AGNI"] = quantslice
    plotSliceTypes["tau"] = zweightslice
    plotSliceTypes["list"] = quantslice
    plotSliceTypes["opac"] = quantslice
    plotSliceTypes["smooth"] = quantslice
    plotSliceTypes["facetemp"] = viewslice
    plotSliceTypes["rand"] = quantslice
    plotSliceTypes["rad0"] = quantslice
    plotSliceTypes["nneigh"] = quantslice
    plotSliceTypes["pres"] = quantslice
    plotSliceTypes["age"] = quantslice
    for line in lines:
        plotSliceTypes[line]=quantslice
        plotSliceTypes[line+"m"]=dslice
#         plotSliceTypes["v"+line]=quantslice
        plotSliceTypes["v"+line]=weightviewslice
        plotSliceTypes["dv"+line]=sdevslice
#         plotSliceTypes["view"+line]=weightviewslice
        plotSliceTypes["view"+line]=viewslice
        plotSliceTypes["vels"+line]=vec2dslice
#     plotSliceTypes["co1"] = quantslice
#     plotSliceTypes["co2"] = quantslice
#     plotSliceTypes["hcn1"] = quantslice
#     plotSliceTypes["hcn2"] = quantslice
#     plotSliceTypes["co1m"] = dslice
#     plotSliceTypes["co2m"] = dslice
#     plotSliceTypes["hcn1m"] = dslice
#     plotSliceTypes["hcn2m"] = dslice
#     plotSliceTypes["vco1"] = quantslice
#     plotSliceTypes["vco2"] = quantslice
#     plotSliceTypes["vhcn1"] = quantslice
#     plotSliceTypes["vhcn2"] = quantslice
#     plotSliceTypes["dvco1"] = sdevslice
#     plotSliceTypes["dvco2"] = sdevslice
#     plotSliceTypes["dvhcn1"] = sdevslice
#     plotSliceTypes["dvhcn2"] = sdevslice
    
    plotCustomMass = dict()
    plotCustomMass["dust"] = "dust"
    plotCustomMass["dusttau"] = "dusttau"
    for line in lines:
        plotCustomMass[line+"m"]=line+"m"
        plotCustomMass["v"+line]=line+"m"
        plotCustomMass["dv"+line]=line+"m"
        plotCustomMass["vels"+line]=line+"m"
#         plotCustomMass["view"+line]=line+"m"
    
    plotData = dict()
    plotData["temp"] = "TK_p"
    plotData["col"] = "coldens"
    plotData["nH"] = "nH_p"
    plotData["heat"] = "heat"
    #plotData["dens"] = None
#     plotData["dusttau"] = "dusttau"
    plotData["depth"] = "depth_p"
    plotData["vels"] = ["vel_x","vel_y","vel2d","vel_z"]
    plotData["tdust"] = "dustTemp"
    plotData["dg"] = "dg"
#    plotData["dust"] = dslice
    plotData["view"] = ["brightness","opac","brightness","opac"]
    plotData["vlos"] = ["vthin","opac","vthin","opac"]
    plotData["emit"] = "emissivity"
    plotData["opac"] = "opac"
    plotData["dt"] = "dt_p"
    plotData["vel_2d"] = "vel2d"
    plotData["vel_r"] = "velr"
    plotData["vmag"] = "vmag"
    plotData["vel_x"] = "vel_x"
    plotData["vel_y"] = "vel_y"
    plotData["vel_z"] = "vel_z"
    plotData["vel_a"] = "vel_a"
    plotData["vthin"] = "vthin"
    plotData["arad"] = "arad"
    plotData["accel"] = "accel"
    plotData["tau"] = "tau"
    plotData["AGNI"] = "AGNI"
    plotData["list"] = "list"
    plotData["rand"] = "rand"
    plotData["smooth"] = "h_p"
    plotData["rad0"] = "rad0"
    plotData["facetemp"] = ["dustTemp","opac","dustTemp","opac"]
    plotData["nneigh"] = "nneigh"
    plotData["pres"] = "pres"
    plotData["age"] = "age"
#     plotData["co1"] = "co1"
#     plotData["co2"] = "co2"
#     plotData["hcn1"] = "hcn1"
#     plotData["hcn2"] = "hcn2"
    for line in lines:
        plotData[line]=line
        plotData[line+"m"]=line+"m"
        plotData["v"+line]=["vthin",line+"brightness",line+"opac","vythin",line+"brightness",line+"opac"]
        plotData["dv"+line]="vthin"
#         plotData["view"+line]=["vthin",line+"m","opac","vel_x",line+"m","opac"]
        plotData["view"+line] = [line+"brightness",line+"opac",line+"brightness",line+"opac"]
        plotData["vels"+line]=["vel_x","vel_y","vel2d","vel_z"]
#     plotData["vco1"] = "vthin"
#     plotData["vco2"] = "vthin"
#     plotData["vhcn1"] = "vthin"
#     plotData["vhcn2"] = "vthin"
#     plotData["dvco1"] = "vthin"
#     plotData["co1m"] = "co1m"
#     plotData["co2m"] = "co2m"
#     plotData["hcn1m"] = "hcn1m"
#     plotData["hcn2m"] = "hcn2m"
    
    logSliceTypes = [   "temp","col","nH","dens","dust",
                        "view","emit","dt","arad","accel",
                        "AGNI","opac","smooth","tdust","pres","dusttau"]
    logSliceTypes+=lines
    logSliceTypes+=[line+"m" for line in lines]
    logSliceTypes+=["view"+line for line in lines]
    symLogSliceTypes={"heat":1.e-4}
    # ,
#                         "co1","co2","hcn1","hcn2","co1m","co2m","hcn1m","hcn2m"
# #                         ,"dvco1","dvco2","dvhcn1","dvhcn2"
#                         ]
#     extraBarTypes = ["heat"]
#     plusMinusTypes = ["heat"]
    extraBarTypes = []
    plusMinusTypes = []

    
    divergingTypes = ["vel_r","vel_x","vel_y","vel_z","vel_2d","vel_a","vthin"]
    divergingTypes+=["v"+line for line in lines]
#     divergingTypes+=["view"+line for line in lines]
#     ,"vco1","vco2","vhcn1","vhcn2"]
    
    customCmaps = dict()
#     customCmaps["heat"] = 'Reds'
    customCmaps["dt"] = cmap+"_r"

    customCmaps2 = dict()
#     customCmaps2["heat"] = 'Blues'
    
    
    return plotLabel, plotRanges, plotSliceTypes, plotCustomMass, plotData, logSliceTypes, extraBarTypes, plusMinusTypes, divergingTypes, symLogSliceTypes, customCmaps, customCmaps2

def load_process_gadget_data(infile,rot,plot_thing,plotData,centredens=False,ringPlot=False,flatPlot=False,maskbounds=None):
#,opac_mu=None):
    need_to_load = list(plot_thing)
    if maskbounds:
        need_to_load.append(maskbounds[0])

    time,data = load_gadget(infile,need_to_load,centredens=centredens)#,opac_mu=opac_mu)

    n = data.h_p.size
    
    if data.binary_positions is not None:
      x_bin = binary_positions[0][0]-binary_positions[1][0]
      y_bin = binary_positions[0][1]-binary_positions[1][1]
      rot[0] = rot[0] - math.atan2(y_bin,x_bin)

        
    # convert to pc
    x = data.xyz[:,0]*1.e3
    y = data.xyz[:,1]*1.e3
    z = data.xyz[:,2]*1.e3
    data.h_p*=1.e3
    
    # rotate
    if ( rot[0]!=0. or rot[1]!=0. ):
#        yr = y*np.cos(rot[1]) - z*np.sin(rot[1])
#        zr = y*np.sin(rot[1]) + z*np.cos(rot[1])
#        z = zr
#
#        y = x*np.sin(rot[0]) + yr*np.cos(rot[0])
#        x = x*np.cos(rot[0]) - yr*np.sin(rot[0])

         xr = x*np.cos(rot[0]) - y*np.sin(rot[0])
         yr = x*np.sin(rot[0]) + y*np.cos(rot[0])
         x = xr
         
         y = yr*np.cos(rot[1]) - z*np.sin(rot[1])
         z = yr*np.sin(rot[1]) + z*np.cos(rot[1])
    
    if ( any(x in plot_thing for x in ["vels","vmag","vel_2d","vel_x","vel_y","vel_z","vel_r","vel_a","vthin","vlos"] + ["v"+line for line in lines] + ["dv"+line for line in lines] + ["vels"+line for line in lines]) ): #+ ["view"+line for line in lines] 
        #vel_mag = np.sqrt(np.sum(data.vels[:,:]**2,1))
        data.vel_a = (-x*data.vels[:,1]+y*data.vels[:,0])/np.sqrt(x**2+y**2)
        data.velr = (x*data.vels[:,0]+y*data.vels[:,1]+z*data.vels[:,2])/np.sqrt(x**2+y**2+z**2)
        data.vel_x = data.vels[:,0]
        data.vel_y = data.vels[:,1]
        data.vel_z = data.vels[:,2]
        data.vmag = np.sqrt(np.sum(data.vels**2,axis=1))
        if rot[0]!=0. or rot[1]!=0.:
            if ( ringPlot ):
                raise NotImplementedError()
            else:
                vyr = data.vel_y*np.cos(rot[1]) - data.vel_z*np.sin(rot[1])
                vzr = data.vel_y*np.sin(rot[1]) + data.vel_z*np.cos(rot[1])
                data.vel_z = vzr

                data.vel_y = data.vel_x*np.sin(rot[0]) + vyr*np.cos(rot[0])
                data.vel_x = data.vel_x*np.cos(rot[0]) - vyr*np.sin(rot[0])

                data.vel2d = data.vel_x
        else:
            if ( ringPlot ):
                data.vel2d = (x*data.vel_x+y*data.vel_y)/np.sqrt(x**2+y**2)
            else:
                data.vel2d = data.vel_x
    
    if any(x in plot_thing for x in ["vthin","vlos"]+["v"+line for line in lines] + ["dv"+line for line in lines]): # + ["view"+line for line in lines]
        if ringPlot:
            raise Exception("can't do optically thin line of sight velocity while integrating around ring")

        # rotation already done above
        data.vthin = data.vel_z
        data.vythin = data.vel_y
        
#         if rot[0]!=0. or rot[1]!=0.:
#             vyr = data.vel_x*np.sin(rot[0]) + data.vel_y*np.cos(rot[0])
#             data.vythin = vyr*np.cos(rot[1]) - data.vel_z*np.sin(rot[1])
# 
#             data.vthin = vyr*np.sin(rot[1]) + data.vel_z*np.cos(rot[1])
#         else:
#             data.vthin = data.vel_z
#             data.vythin = data.vel_y
        


#         for vthinplot in "vco1","vco2","vhcn1","vhcn2":
#             if vthinplot in plot_thing:
#                 data[vthinplot]*=data.vthin

    # set x coordinate - cartesian or cylindrical ("ringplot")
    if ( ringPlot ):
        if ( not flatPlot ):
            raise Exception("ring==true requires flat==true")
        rad2d = np.sqrt(x**2+y**2)
    else:   
        rad2d = x

    deep_face = z
    deep_side = y
    
    if maskbounds:
        if maskbounds[0] in plotData:
            v=data.__dict__[plotData[maskbounds[0]]]
            mask = (v>maskbounds[1]) & (v<maskbounds[2])
        else:
            raise Exception("mask value {} not found in plotData".format(maskbounds[0]))
    else:
        # mask out non-gas - currently everything is gas
        mask = np.full(n,True,dtype=bool)
        
    # flat weighting - dummy value required for some functions because I'm not qwarging properly yet
    n_ones = np.ones(n)

    return time,data,x,y,z,rad2d,deep_face,deep_side,mask,n,n_ones


def makesph_trhoz_frame(*args,**kwargs):
    try:
        return makesph_trhoz_frame_wrapped(*args,**kwargs)
    except Exception as e:
        print("returning exception")
        wropped = ExceptionWrapper(e)
        print("wrapped")
        return wropped
#         print(repr(e))
#         print("".join(traceback.format_exception(*sys.exc_info())))
#         print(str(e))
#         print(sys.exc_info())
#         traceback.print_exc(file=sys.stdout)
#         print(traceback.format_tb(sys.exc_info()[2]))
#         traceback.print_stack()
#         traceback.print_tb(sys.exc_info()[2])
#         print("Raising back up")
#         raise(e)
#         raise Exception("".join(traceback.format_exception(*sys.exc_info())))



def makesph_trhoz_frame_wrapped(infile,outfile,
# def makesph_trhoz_frame(infile,outfile,
                        scale = .7,
                        cmap="viridis",
                        L=256,
                        ring=False,
                        flat=False,
                        planenorm=False,
                        visibleAxes=True,
                        subsample=1,
                        pixsize=None,
                        plot=['dens','temp'],
                        rot=[0.,0.],
                        views=['face','side'],
                        centredens=False,
                        centrecom=False,
                        vorinoi=False,
                        maxmode=False,
                        dotmode=None,
                        titlesuffix="",
                        maskbounds=None,
                        data_ranges=None,
                        gaussian=None,
                        return_maps=False
                        #,opac_mu=None
                        ):
    if version[0]=='2':
        raise Exception("Requires Python 3.X. Current version:"+version)

    #cmap_r = cmap+"_r"
    #cmap_d = "coolwarm"
    plot_thing = plot
    flatPlot = flat
    ringPlot = ring
    
    if return_maps:
        out_maps = []
    
    if dotmode!='max' and dotmode!='min':
        dotmode = None

    if ( L%subsample!=0 ):
        raise Exception("subsample might divide evenly into L")

    if ( not pixsize ):
        pixsize = subsample

    if (len(rot)!=2):
        raise Exception("rot needs to be [theta,phi]")
        
    if ( not all(view=='face' or view=='side' for view in views) ):
        raise Exception("Views must be an array of size one or two containing only 'face' or 'side'")
    if ( not 'side' in views ):
        ringPlot = False # IRRELEVANT
    
    
    cols = len(views)
    if ( cols!=1 and cols!=2 ):
        raise Exception("len(views)=1 or =2")
    
    if data_ranges is not None:
        if len(data_ranges)!=len(plot):
            raise Exception("Need to give data ranges for ALL plots or none")

    plotLabel, plotRanges, plotSliceTypes, plotCustomMass, plotData, logSliceTypes, extraBarTypes, plusMinusTypes, divergingTypes, symLogSliceTypes, customCmaps, customCmaps2 = pack_dicts(flatPlot,planenorm,cmap)

    nrows = len(plot_thing)

    if ( visibleAxes ):
        # figure properties
#         if ( "heat" in plot_thing ):
#             fig, ax = P.subplots(nrows,3*cols, gridspec_kw = {'width_ratios':([1, 1, 16,16,1, 1])[0:3*cols]})
#             if ( nrows<2 ):
#                 ax = np.resize(ax, (3,6)) # because 1D arrays have different syntax, we have to pretend it's 3D
#             cbax2left_index = 0
#             cbaxleft_index = 1
#             spleft_index = 2
#             spright_index = 3
#             cbaxright_index = 4
#             cbax2right_index = 5
#         
#         else:
            fig, ax = P.subplots(nrows,2*cols, gridspec_kw = {'width_ratios':([1, 16,16,1])[0:2*cols]})
            if ( nrows<2 ):
                ax = np.resize(ax, (3,4)) # because 1D arrays have different syntax, we have to pretend it's 3D
            cbaxleft_index = 0
            spleft_index = 1
            spright_index = 2
            cbaxright_index = 3
    else:
        fig, ax = P.subplots(nrows,cols)
        if ( nrows<2 ):
            ax = np.resize(ax, (3,4)) # because 1D arrays have different syntax, we have to pretend it's 3D
        cbaxleft_index = 0
        spleft_index = 0
        spright_index = 1
        cbaxright_index = 1
        
    if not isinstance(infile,list):
        verboseprint("Loading",infile)
        time,data,x,y,z,rad2d,deep_face,deep_side,mask,n,n_ones = load_process_gadget_data(infile,rot,plot_thing,plotData,centredens,ringPlot,flatPlot,maskbounds)#,opac_mu=opac_mu)
#         time,data = load_gadget(infile,need_to_load,centredens=centredens)
        verboseprint("Plotting",infile,", t=%.4f Myr"% time)

        if ( visibleAxes ):
            if time<1.e-3:
                fig.suptitle(r"$T="+("%.4f" % (time*1.e3))+"$ kyr"+titlesuffix)
            else:
                fig.suptitle(r"$T="+("%.4f" % time)+"$ Myr"+titlesuffix)
#         fig.suptitle(infile)

    if ( visibleAxes ):
        fw_inches = 5.*cols
    else:
        fw_inches = 4.*cols
    fig.set_figwidth(fw_inches)
    fig.set_figheight(4.*nrows)

    # nerrors=0

    # do all subplots, calculating the full SPH smoothing each time
    for irow in range(nrows):
        
        if isinstance(infile,list):
            thisfile = infile[irow]
            verboseprint("Loading",thisfile)
            time,data,x,y,z,rad2d,deep_face,deep_side,mask,n,n_ones = load_process_gadget_data(thisfile,need_to_load,centredens,rot,plot_thing,plotData,ringPlot,flatPlot,maskbounds)
            verboseprint("Plotting",thisfile,", t=%.4f Myr"% time)

        # Process the tables we just made
        # Believe it or not, this is clearer than writing out each
        # plot command by hand - there's less repetition
        if ( not plot_thing[irow] in plotLabel ):
            errstr = "{} is not a valid plot type\n".format(plot_thing[irow])
            matchRatios = list(map(similar,itertools.repeat(plot_thing[irow]),plotLabel))
            bestFit = np.argmax(matchRatios)
            errstr+= "Valid plot types:\n"
            labelstrs = list(plotLabel)
            for i in range(len(labelstrs)):
                if ( i!=bestFit ):
                    if ( matchRatios[i]>.5 ):
                        labelstrs[i]+=" - possible fit?"
            labelstrs[bestFit]+=" - best match"
            errstr += "\n".join(labelstrs)
            raise Exception(errstr)
        thisPlotLabel = plotLabel[plot_thing[irow]]
        if data_ranges is not None:
            thisPlotRanges = data_ranges[irow]*2
        else:
            thisPlotRanges = plotRanges[plot_thing[irow]]
#             thisPlotRanges = plotRanges[plot_thing[irow]]*2

#         if isinstance(opac_mu,list):
            

        if isinstance(gaussian,list):
            plot_gaussian = gaussian[irow]
        else:
            plot_gaussian = gaussian


        # physical coordinates of region to plot, in pc
        if isinstance(scale,list):
            width = scale[irow]*2.
        else:
            width = scale*2.
        if centredens or centrecom :
            if centredens and centrecom :
                raise Exception("Can't set both centredens and centrecom")
            if centredens:
                i_maxdens = np.argmax(data.nH_p)
                corners = [x[i_maxdens]-width/2,y[i_maxdens]-width/2.]
                corners_side = [rad2d[i_maxdens]-width/2.,z[i_maxdens]-width/2.]
            if centrecom:
                x_com = np.mean(x)
                y_com = np.mean(y)
                corners = [x_com-width/2.,y_com-width/2.]
                r2d_com = np.mean(rad2d)
                z_com = np.mean(z)
                corners_side = [r2d_com-width/2.,z_com-width/2.]
        else:
            # centre on 0,0
            corners = [-width/2.,-width/2.]
            if ( ringPlot ):
                if ( not flatPlot ):
                    raise Exception("ring==true requires flat==true")
                corners_side = [0.,-width/2.]
            else:   
                corners_side = corners

        if dotmode=='max':
            thisSliceType = maxdotslice
        elif dotmode=='min':
            thisSliceType = mindotslice
        elif maxmode:
            thisSliceType = mslice
        elif vorinoi:
            thisSliceType = vslice
        else:
            thisSliceType = plotSliceTypes[plot_thing[irow]]
        thisDoLog = (plot_thing[irow] in logSliceTypes)
        thisDiverging = (plot_thing[irow] in divergingTypes)
        if plot_thing[irow] in symLogSliceTypes:
            thisSymLog = symLogSliceTypes[plot_thing[irow]]
        else:
            thisSymLog = None
        
        if ( plot_thing[irow] in plotCustomMass ):
            thisMass = data.__dict__[plotCustomMass[plot_thing[irow]]]
        else:
            thisMass = data.m_p

        if ( plot_thing[irow] in plusMinusTypes ):
            plusminus = True
        else:
            plusminus = False
        
        if ( plot_thing[irow] in plotData ):
            plotCommand = plotData[plot_thing[irow]]
            if ( type(plotCommand) is str ):
                thisPlotQuantityFace = data.__dict__[plotCommand]
                thisPlotQuantitySide = thisPlotQuantityFace
            elif ( type(plotCommand) is list ):
                if thisSliceType==viewslice or thisSliceType==vec2dslice_flat or thisSliceType==zvec2dslice:
                    if ( len(plotCommand)!=4 ):
                        raise Exception("{} must have length 4 for view or vec2d slice".format(plotCommand))
                    thisPlotQuantityFace = [data.__dict__[plotCommand[0]],data.__dict__[plotCommand[1]]]
                    thisPlotQuantitySide = [data.__dict__[plotCommand[2]],data.__dict__[plotCommand[3]]]
                elif thisSliceType==weightviewslice:
                    if ( len(plotCommand)!=6 ):
                        raise Exception("{} must have length 6 for weighted view slice".format(plotCommand))
                    thisPlotQuantityFace = [data.__dict__[plotCommand[0]],data.__dict__[plotCommand[1]],data.__dict__[plotCommand[2]]]
                    thisPlotQuantitySide = [data.__dict__[plotCommand[3]],data.__dict__[plotCommand[4]],data.__dict__[plotCommand[5]]]
                else:
                    raise Exception("{} can't be an array except for view and vec2d slices".format(plotCommand))
            else:
                raise Exception("{} is not a valid plot format".format(plotCommand))
        else:
            thisPlotQuantityFace = None
            thisPlotQuantitySide = None
#             raiseprint("{} data type not found! - {}".format(plot_thing[irow],infile))
        
        cbar2_axes = [None,None]
        if ( plot_thing[irow] in extraBarTypes ):
            cbar2_axes[0]=ax[irow,cbax2left_index]
            if ( cols==2 ):
                cbar2_axes[1]=ax[irow,cbax2right_index]
        
        if ( plot_thing[irow] in customCmaps ):
            this_cmap = customCmaps[plot_thing[irow]]
        else:
            this_cmap = cmap

        if ( plot_thing[irow] in customCmaps2 ):
            this_cmap2 = customCmaps2[plot_thing[irow]]
        else:
            this_cmap2 = None

        row_axes = [ax[irow,spleft_index],ax[irow,cbaxleft_index]]
        if ( cols==2 ):
            row_axes += ax[irow,spright_index],ax[irow,cbaxright_index]

        # actually do the plot
        for icol,view in enumerate(views):
            if ( view=='face' ):
                outmap=makesph_plot(fig,row_axes[icol*2],row_axes[icol*2+1],x,y,deep_face,0.,thisPlotQuantityFace,thisMass,data.h_p,L,mask,corners,width,thisPlotLabel,thisPlotRanges[0],thisPlotRanges[1],this_cmap,thisSliceType,thisDoLog,
                        cmap2=this_cmap2,circnorm=planenorm,cbar2=cbar2_axes[icol],plusminus=plusminus,visibleAxes=visibleAxes,diverging=thisDiverging,gaussian=plot_gaussian,symLog=thisSymLog)
                if(data.binary_positions):
                  row_axes[icol*2].scatter([data.binary_positions[0][0]],[data.binary_positions[1][1]],marker='x')
                  row_axes[icol*2].scatter([data.binary_positions[1][0]],[data.binary_positions[1][1]],marker='+')
            elif ( view=='side' ):
                outmap=makesph_plot(fig,row_axes[icol*2],row_axes[icol*2+1],rad2d,z,deep_side,0.,thisPlotQuantitySide,thisMass,data.h_p,L,mask,corners_side,width,thisPlotLabel,thisPlotRanges[2],thisPlotRanges[3],this_cmap,thisSliceType,thisDoLog,
                        cmap2=this_cmap2,planenorm=planenorm,cbar2=cbar2_axes[icol],plusminus=plusminus,visibleAxes=visibleAxes,diverging=thisDiverging,gaussian=plot_gaussian,symLog=thisSymLog)
            if return_maps:
                out_maps.append(outmap)

#         if ( "heat" in plot_thing and plot_thing[irow]!='heat' ):
#             for visax in [ax[irow,cbax2left_index],ax[irow,cbax2right_index]]:
#                 visax.set_frame_on(False)
#                 visax.axes.get_yaxis().set_visible(False)
#                 visax.axes.get_xaxis().set_visible(False)


    if ( visibleAxes ):
        if cols==2:
            for iax in range(nrows):
                this_ax = ax[iax,1]
                this_ax.yaxis.tick_right()
                this_ax.yaxis.set_visible(False)
            for iax in range(nrows):
                this_ax = ax[iax,0]
                this_ax.yaxis.tick_left()
                this_ax.yaxis.set_label_position("left")

        else:
            for iax in range(nrows):
                this_ax = ax[iax,0]
                this_ax.yaxis.tick_left()
                this_ax.yaxis.set_label_position("left")
            for iax in range(nrows):
                this_ax = ax[iax,1]
                this_ax.yaxis.tick_right()
                this_ax.yaxis.set_label_position("right")
#                 this_ax.set_ylabel('pc')
# 
#                 this_ax.set_xlabel('pc')
#                 this_ax.yaxis.set_visible(False)

    else:
        P.axis('off')

#     print(infile,nerrors)

    if ( visibleAxes ):
        if (nrows==2):
            fig.subplots_adjust(left=0.07,hspace=.07,bottom=.05,top=.95)
        else:
#             fig.subplots_adjust(left=0.07,hspace=.07,bottom=.05,top=.9)
            fig.subplots_adjust(left=0.12,hspace=.12,bottom=.1,top=.9)
    else:
        fig.subplots_adjust(left=0.0,hspace=.0,top=1.,bottom=.0,right=1.,wspace=0.)
        
    pos = ax[0,1].get_position()
    lpixx = (L/(pos.x1-pos.x0))
    my_dpi = int(np.floor(lpixx/fw_inches))*pixsize
    
    fig.savefig(outfile,dpi=my_dpi)
    P.close()
    if return_maps:
        return out_maps
    else:
        return None
