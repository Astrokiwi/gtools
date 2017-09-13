# Import libraries to do our magic
import numpy as np
import h5py
import matplotlib as mpl
mpl.use('Agg')

from matplotlib import colors

import pylab as P
#from joblib import Parallel, delayed

from sph_plotter import sph_plotter

from sys import path
path.append("../src/")
import tab_interp

import numpy as np

from difflib import SequenceMatcher

import itertools

#from enum import Enum

def similar(a, b):
    return SequenceMatcher(None, a, b).ratio()

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
vec2dslice = 4
viewslice = 5
minslice = 6
zminslice = 7

molecular_mass = 4./(1.+3.*.76)
proton_mass_cgs = 1.6726e-24
gamma_minus_one = 5./3.-1.
boltzmann_cgs = 1.38066e-16


chTab = None
interpTabVec = None


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

def makesph_plot(fig,sp,cbax,x_p,y_p,z_p,zslice,val_p,m_p,h_p,L,mask,corner,width,cblabel,clow,chigh,cmap_label,mode,
                 dolog=True,
                 planenorm=False,
                 circnorm=False,
                 plusminus=False,
                 diverging=False,
                 visibleAxes=True,
                 cbar2=None,
                 cmap2=None):

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
    
    if ( mode==weightslice ):
        map = sph_plotter.sph_weight(x_p,y_p,m_p,h_p,val_p,L,corner,width,mask,n)
    elif (mode==densslice ):
        map = sph_plotter.sph_dense(x_p,y_p,m_p,h_p,L,corner,width,mask,n)
    elif (mode==zdensslice ):
        map = sph_plotter.sph_dense_slice(x_p,y_p,m_p,h_p,L,corner,width,z_p,zslice,mask,n)
    elif (mode==zweightslice ):
        map = sph_plotter.sph_weight_slice(x_p,y_p,m_p,h_p,val_p,L,corner,width,z_p,zslice,mask,n)
    elif (mode==vec2dslice ):
        map1 = sph_plotter.sph_weight(x_p,y_p,m_p,h_p,val_p[0],L,corner,width,mask,n)
        map2 = sph_plotter.sph_weight(x_p,y_p,m_p,h_p,val_p[1],L,corner,width,mask,n)
        norm_map = np.sqrt(map1**2+map2**2)
        map1/=norm_map
        map2/=norm_map
    elif (mode==viewslice ):
        zarg = np.argsort(z_p)
        map = sph_plotter.sph_optical_depth_los(x_p,y_p,m_p,h_p,val_p[0],val_p[1],L,corner,width,z_p,zarg,mask,n)
    elif (mode==minslice ):
        map = sph_plotter.sph_min(x_p,y_p,h_p,val_p,L,corner,width,mask,n)
    elif (mode==zminslice ):
        map = sph_plotter.sph_minslice(x_p,y_p,h_p,val_p,L,corner,width,mask,n)
        
    if ( mode==vec2dslice ):
        #map1/=1.e5
        #map2/=1.e5
        map1 = map1.T
        map2 = map2.T
        
        norm_map = norm_map.T
        
        step = width/L
        xmids = np.arange(corner[0]+step,corner[0]+width+step,step)
        ymids = np.arange(corner[1]+step,corner[1]+width+step,step)
        sp.set_axis_bgcolor('black')
        norm_map = np.log10(norm_map)
        #qv = sp.quiver(xmids,ymids,map1,map2,norm_map,headwidth=10.,pivot='mid',cmap=this_cmap,clim=[clow,chigh])
        xedges = np.arange(corner[0],corner[0]+width,width/L)
        yedges = np.arange(corner[1],corner[1]+width,width/L)

        qv = sp.pcolormesh(xedges,yedges,norm_map,cmap=this_cmap,vmin=clow,vmax=chigh)
        sp.streamplot(xmids,ymids,map1,map2)
        if ( visibleAxes ):
            cb = fig.colorbar(qv,label=cblabel,cax=cbax)
        else:
            cbax.set_axis_off()
#         sp.set_xbound(corner[0],corner[0]+width)
#         sp.set_ybound(corner[1],corner[1]+width)
#         sp.set_xlim(corner[0],corner[0]+width)
#         sp.set_ylim(corner[1],corner[1]+width)
#         sp.set_xmargin(0)
#         sp.set_ymargin(0)
    else:
        xedges = np.arange(corner[0],corner[0]+width,width/L)
        yedges = np.arange(corner[1],corner[1]+width,width/L)
#         step = width/L
#         xedges = np.arange(corner[0]+step,corner[0]+width+step,step)
#         yedges = np.arange(corner[1]+step,corner[1]+width+step,step)
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
                mesh1=sp.pcolormesh(xedges,yedges,plusmap.T,cmap=this_cmap,vmin=clow,vmax=chigh,norm = colors.LogNorm())
                if ( visibleAxes ):
                    cb = fig.colorbar(mesh1,label=cblabel,cax=cbax)
                else:
                    cbax.set_axis_off()
            
            if ( np.sum(isminus)>0 ):
                minusmap[isplus]=0.
                minusmap = -minusmap
            
                mesh2=sp.pcolormesh(xedges,yedges,minusmap.T,cmap=this_cmap2,vmin=clow,vmax=chigh,norm = colors.LogNorm())
                if ( visibleAxes ):
                    cb2 = fig.colorbar(mesh2,label=cblabel,cax=cbax2)
                else:
                    cbax2.set_axis_off()
        else:

            if ( dolog ):
                map = np.log10(map)
            finiteIndices = np.isfinite(map)
            if ( np.sum(finiteIndices)==0 ):
                print(map)
                raise Exception("No finite values in map")
            print(np.min(map[finiteIndices]),np.max(map[finiteIndices]))

            mesh = sp.pcolormesh(xedges,yedges,map.T,cmap=this_cmap,vmin=clow,vmax=chigh)
            if ( visibleAxes ):
                cb = fig.colorbar(mesh,label=cblabel,cax=cbax)
            else:
                cbax.set_axis_off()
    
    sp.set_xlim([corner[0],corner[0]+width])
    sp.set_ylim([corner[1],corner[1]+width])
    if ( visibleAxes ):
        sp.plot([0],[0],'+g',markersize=10.,markeredgewidth=1.)
    #mesh = sp.pcolormesh(xedges,yedges,map.T,cmap=this_cmap) 
    #sp.axis('equal')

def load_gadget(infile, plot_thing):
    global chTab
    global interpTabVec

    data = GadgetData()

    f = h5py.File(infile,"r")
    
    header = f["/Header"]
    time = header.attrs.get("Time")
    time*= 0.9778e9 # to yr
    time/=1.e6 # to Myr

    data.xyz = np.array(f["/PartType0/Coordinates"]) # kpc

    data.m_p = np.array(f["/PartType0/Masses"]) # 10^10 msun
    data.m_p*=1e+10 # 10^10 solar masses to solar masses
    
    n = data.m_p.size

    data.h_p = np.array(f["/PartType0/SmoothingLength"]) # kpc
    
    need_to_load = list(plot_thing)
    if ( "view" in need_to_load ):
        need_to_load.append("tdust")
    if ( "vlos" in need_to_load or "vmag" in need_to_load ):
        need_to_load.append("vels")
    if ( "emit" in need_to_load ):
        need_to_load.append("tdust")
    if ( "tdust" in need_to_load ):
        need_to_load.append("table")
    if ( "dust" in need_to_load ):
        need_to_load.append("dg")
    if ( "dg" in need_to_load ):
        need_to_load.append("table")
    if ( "table" in need_to_load ):
        need_to_load.append("col")
        need_to_load.append("temp")
        need_to_load.append("nH")

    if ( "vel_2d" in need_to_load or
         "vel_r" in need_to_load or
         "vel_x" in need_to_load or
         "vel_y" in need_to_load or
         "vel_z" in need_to_load or
         "vel_a" in need_to_load):
        need_to_load.append("vels")

    if ( "arad" in need_to_load ):
        data.arads = np.array(f["/PartType0/RadiativeAcceleration"])
        data.arads*=3.24086617e-12 # to cm/s/s
        data.arad = np.sqrt(np.sum(data.arads**2,axis=1))

    if ( "col" in need_to_load ):
        data.coldens = np.array(f["/PartType0/AGNColDens"]) # Msun/kpc**2
        data.coldens*=(1.989e+43/3.086e+21**2) # to g/cm**2
        data.coldens/=(molecular_mass*proton_mass_cgs) # N in cm**(-2) 
    
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

    if ( "table" in need_to_load ):
        if ( not chTab ):
            print("Load dust tables")
            chTab = tab_interp.CoolHeatTab(("../coolheat_tab_marta/shrunk_table_labels_130617.dat"),("../coolheat_tab_marta/shrunk_table_130617.dat"))
            interpTabVec = np.vectorize(chTab.interpTab)
        data.flux_p = np.array(f["/PartType0/AGNIntensity"]) # energy per surface area per time
        data.flux_p*=1.989e+53/(3.086e21)**2/(3.08568e+16)
        

        print("Calculating dust/cooling/heating properties from table")
        tabStructs = interpTabVec(data.nH_p.astype(np.float64),data.TK_p.astype(np.float64),data.flux_p.astype(np.float64),data.coldens.astype(np.float64))
        

    if ( "tdust" in need_to_load ):
        data.dustTemp = map(lambda y: y.dustT, tabStructs)
        data.dustTemp = np.array(list(data.dustTemp))

    if ( "dg" in need_to_load ):
        data.dg = map(lambda y: y.dg, tabStructs)
        data.dg = np.array(list(data.dg))

    if ( "dust" in need_to_load ):
        data.dust = data.dg * data.m_p
    
    if ( "emit" in need_to_load ):
        data.emissivity = 1.e-20*data.m_p * data.dustTemp**4. # arbitrary units

    if ( "view" in need_to_load or "vlos" in need_to_load ):
        if ( "view" in need_to_load ):
            data.brightness = 5.67e-5 * data.dustTemp**4. # erg/s/cm^2
        #opacity = 652. # cm^2/g, somewhat arbitrary
        #opacity = 1. # cm^2/g, somewhat arbitrary
        opacity = 65.2 # cm^2/g, somewhat arbitrary
        opacity*=0.000208908219 # convert to pc**2/solar mass for consistency
        data.opac = np.full((n),opacity)
        
        if ( "view" in need_to_load ):
            sputtered = (data.dustTemp>2.5e3) # or something, super arbitrary
            data.brightness[sputtered] = 0.
            data.opac[sputtered] = 0.

    #rho_p*=6.77e-22 # to g/cm**3
    #m_p*=1.989e+43 # 10^10 solar masses to g
    #agn_heat_p*=1e10/3.08568e+16# to erg/s/g
    #agn_heat_p*=rho_p# to erg/cm**3/s
    
#     if ( 'agn_heat_p' in globals() ):
#         depth_p = agn_heat_p*rho_p*h_p*(16./3./np.pi)

    
    return time,data

def makesph_trhoz_frame(infile,outfile,
                        scale = 20.,
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
                        views=['face','side']
                        ):
    #cmap_r = cmap+"_r"
    #cmap_d = "coolwarm"
    plot_thing = plot
    flatPlot = flat
    ringPlot = ring
    width = scale*2.

    if ( L%subsample!=0 ):
        raise Exception("subsample might divide evenly into L")

    if ( not pixsize ):
        pixsize = subsample

    if (len(rot)!=2):
        raise Exception("rot needs to be [theta,phi]")
        
    if ( not all(view=='face' or view=='side' for view in views) ):
        raise Exception("Views must be an array of size one or two containing only 'face' or 'side'")
    
    cols = len(views)
    if ( cols!=1 and cols!=2 ):
        raise Exception("len(views)=1 or =2")

    nrows = len(plot_thing)

    print("Loading",infile)
    time,data = load_gadget(infile,plot_thing)
    print("Plotting",infile,", t=%.4f Myr"% time)



    n = data.h_p.size
        
    # convert to pc
    x = data.xyz[:,0]*1.e3
    y = data.xyz[:,1]*1.e3
    z = data.xyz[:,2]*1.e3
    data.h_p*=1.e3
    
    # rotate
    if ( rot[0]!=0. or rot[1]!=0. ):
        x  = x*np.cos(rot[0]) - y*np.sin(rot[0])
        yr = x*np.sin(rot[0]) + y*np.cos(rot[0])
        
        y = yr*np.cos(rot[1]) - z*np.sin(rot[1])
        z = yr*np.sin(rot[1]) + z*np.cos(rot[1])
    
    
    if ( any(x in plot_thing for x in ["vels","vmag","vel_2d","vel_x","vel_y","vel_z","vel_r","vel_a"]) ):
        #vel_mag = np.sqrt(np.sum(data.vels[:,:]**2,1))
        data.vel2d = (x*data.vels[:,0]+y*data.vels[:,1])/np.sqrt(x**2+y**2)
        data.vel_a = (-x*data.vels[:,1]+y*data.vels[:,0])/np.sqrt(x**2+y**2)
        data.velr = (x*data.vels[:,0]+y*data.vels[:,1]+z*data.vels[:,2])/np.sqrt(x**2+y**2+z**2)
        data.vel_x = data.vels[:,0]
        data.vel_y = data.vels[:,1]
        data.vel_z = data.vels[:,2]

    # mask out non-gas - currently everything is gas
    mask = np.full(n,True,dtype=bool)
    
    # flat weighting - dummy value required for some functions because I'm not qwarging properly yet
    n_ones = np.ones(n)

    if ( visibleAxes ):
        # figure properties
        if ( "heat" in plot_thing ):
            fig, ax = P.subplots(nrows,3*cols, gridspec_kw = {'width_ratios':([1, 1, 16,16,1, 1])[0:3*cols]})
            if ( nrows<2 ):
                ax = np.resize(ax, (3,6)) # because 1D arrays have different syntax, we have to pretend it's 3D
            cbax2left_index = 0
            cbaxleft_index = 1
            spleft_index = 2
            spright_index = 3
            cbaxright_index = 4
            cbax2right_index = 5
        
        else:
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
        
        

    if ( visibleAxes ):
        fig.suptitle(r"$T="+("%.4f" % time)+"$ Myr")
    if ( visibleAxes ):
        fw_inches = 5.*cols
    else:
        fw_inches = 4.*cols
    fig.set_figwidth(fw_inches)
    fig.set_figheight(4.*nrows)

    # physical coordinates of region to plot, in pc
    corners = [-width/2.,-width/2.]

    if ( ringPlot ):
        if ( not flatPlot ):
            raise Exception("ring==true requires flat==true")
        corners_side = [0.,-width/2.]
        rad2d = np.sqrt(x**2+y**2)
    else:   
        corners_side = corners
        rad2d = x

    deep_face = z
    deep_side = y
    if ( flatPlot ):
        quantslice = weightslice
        dslice = densslice
        rhounit = r"$\log_{10} \Sigma$ (M$_\odot$/pc$^2$)"
        drange = [-2.,8.]
        thisminslice = minslice
        #drange = [2.,6.]
        #drange = [-1.,4.]
    else:
        quantslice = zweightslice
        dslice = zdensslice
        thisminslice = zminslice
        rhounit = r"$\log_{10} \rho$ (M$_\odot$/pc$^3$)"
        drange = [-2.,7.]
    
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
    plotLabel["depth"] = r"$\tau$"
    plotLabel["vels"] = r"$\log_{10}v$ (km/s)"
    plotLabel["tdust"] = r"$\log_{10} T_d$ (K)"
    plotLabel["dg"] = r"$f_d$"
    plotLabel["dust"] = rhounit
    plotLabel["view"] = r"$\log_{10} F$"
    plotLabel["vlos"] = r"$v_{LOS}$"
    plotLabel["emit"] = r"$\log_{10} F$"
    plotLabel["dt"] = r"$\log_{10}$ dt"
    plotLabel["vel_2d"] = r"$v_{2D}$"
    plotLabel["vel_x"] = r"$v_{x}$"
    plotLabel["vel_y"] = r"$v_{y}$"
    plotLabel["vel_z"] = r"$v_{z}$"
    plotLabel["vel_r"] = r"$v_{r}$"
    plotLabel["vel_a"] = r"$v_{\theta}$"
    plotLabel["arad"] = r"$a_{rad}$ (log cm/s/s)"

    plotRanges = dict()
    plotRanges["temp"] = [.999,6.]*2
    plotRanges["col"] = [15.,26.]*2
    plotRanges["nH"] = [0.,8.]*2
    plotRanges["heat"] = [1e-2,1.e10]*2
    plotRanges["dens"] = drange+drange_s
    plotRanges["depth"] = [0.,1.e3]*2
    plotRanges["vels"] = [0.,3.]*2
    plotRanges["tdust"] = [1.,3.]*2
    plotRanges["dg"] = [0.00635,.006363]*2
    plotRanges["dust"] = drange+drange_s
    plotRanges["view"] = [0.,4.]*2
    plotRanges["vlos"] = [-1.2e2,1.2e2]*2
    plotRanges["emit"] = [-12.,-3.]*2
    plotRanges["dt"] = [1.,3.]*2
    #plotRanges["vel_2d"] = [-25.,25.]*2
    plotRanges["vel_2d"] = [-250.,250.]*2
    plotRanges["vel_r"] = [-25.,25.]*2
    plotRanges["vel_x"] = [-200.,200.]*2
    plotRanges["vel_y"] = [-200.,200.]*2
    plotRanges["vel_z"] = [-25.,25.]*2
    plotRanges["vel_a"] = [-100.,100.]*2
    plotRanges["arad"] = [-9.,-3.]*2

    plotSliceTypes = dict()
    plotSliceTypes["temp"] = quantslice
    plotSliceTypes["col"] = quantslice
    plotSliceTypes["nH"] = quantslice
    plotSliceTypes["heat"] = quantslice
    plotSliceTypes["dens"] = dslice
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
    plotSliceTypes["vel_x"] = quantslice
    plotSliceTypes["vel_y"] = quantslice
    plotSliceTypes["vel_z"] = quantslice
    plotSliceTypes["vel_a"] = quantslice
    plotSliceTypes["arad"] = quantslice
    
    plotCustomMass = dict()
    plotCustomMass["dust"] = "dust"
    
    plotData = dict()
    plotData["temp"] = "TK_p"
    plotData["col"] = "coldens"
    plotData["nH"] = "nH_p"
    plotData["heat"] = "heat"
    #plotData["dens"] = None
    plotData["depth"] = "depth_p"
    plotData["vels"] = ["vel_x","vel_y","vel2d","vel_z"]
    plotData["tdust"] = "dustTemp"
    plotData["dg"] = "dg"
#    plotData["dust"] = dslice
    plotData["view"] = ["brightness","opac","brightness","opac"]
    plotData["vlos"] = ["vel_z","opac","vel_x","opac"]
    plotData["emit"] = "emmissivity"
    plotData["dt"] = "dt_p"
    plotData["vel_2d"] = "vel2d"
    plotData["vel_r"] = "velr"
    plotData["vel_x"] = "vel_x"
    plotData["vel_y"] = "vel_y"
    plotData["vel_z"] = "vel_z"
    plotData["vel_a"] = "vel_a"
    plotData["arad"] = "arad"
    
    logSliceTypes = ["temp","col","nH","dens","tdust","dust","view","emit","dt","arad"]
    extraBarTypes = ["heat"]
    plusMinusTypes = ["heat"]
    
    divergingTypes = ["vel_r","vel_x","vel_y","vel_z","vel_2d","vel_a"]
    
    customCmaps = dict()
    customCmaps["heat"] = 'Reds'
    customCmaps["dt"] = cmap+"_r"

    customCmaps2 = dict()
    customCmaps2["heat"] = 'Blues'

    
    
    # do all subplots, calculating the full SPH smoothing each time
    for irow in range(nrows):
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
        thisPlotRanges = plotRanges[plot_thing[irow]]
        thisSliceType = plotSliceTypes[plot_thing[irow]]
        thisDoLog = (plot_thing[irow] in logSliceTypes)
        thisDiverging = (plot_thing[irow] in divergingTypes)
        
        if ( plot_thing[irow] in plotCustomMass ):
            thisMass = data.__dict__[plot_thing[irow]]
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
                if ( len(plotCommand)!=4 ):
                    raise Exception("{} must have length 4".format(plotCommand))
                thisPlotQuantityFace = [data.__dict__[plotCommand[0]],data.__dict__[plotCommand[1]]]
                thisPlotQuantitySide = [data.__dict__[plotCommand[2]],data.__dict__[plotCommand[3]]]
            else:
                raise Exception("{} is not a valid plot format".format(plotCommand))
        else:
            thisPlotQuantityFace = None
            thisPlotQuantitySide = None
        
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
                makesph_plot(fig,row_axes[icol*2],row_axes[icol*2+1],x,y,deep_face,0.,thisPlotQuantityFace,thisMass,data.h_p,L,mask,corners,width,thisPlotLabel,thisPlotRanges[0],thisPlotRanges[1],this_cmap,thisSliceType,thisDoLog,
                        cmap2=this_cmap2,circnorm=planenorm,cbar2=cbar2_axes[icol],plusminus=plusminus,visibleAxes=visibleAxes,diverging=thisDiverging)
            elif ( view=='side' ):
                makesph_plot(fig,row_axes[icol*2],row_axes[icol*2+1],rad2d,z,deep_side,0.,thisPlotQuantitySide,data.m_p,data.h_p,L,mask,corners_side,width,thisPlotLabel,thisPlotRanges[2],thisPlotRanges[3],this_cmap,thisSliceType,thisDoLog,
                        cmap2=this_cmap2,planenorm=planenorm,cbar2=cbar2_axes[icol],plusminus=plusminus,visibleAxes=visibleAxes,diverging=thisDiverging)



        if ( "heat" in plot_thing and plot_thing[irow]!='heat' ):
            for visax in [ax[irow,cbax2left_index],ax[irow,cbax2right_index]]:
                visax.set_frame_on(False)
                visax.axes.get_yaxis().set_visible(False)
                visax.axes.get_xaxis().set_visible(False)

    if ( visibleAxes ):
        for iax in range(nrows):
            this_ax = ax[iax,1]
            this_ax.yaxis.tick_right()
            this_ax.yaxis.set_visible(False)

        for iax in range(nrows):
            this_ax = ax[iax,0]
            this_ax.yaxis.tick_left()
            this_ax.yaxis.set_label_position("left")
    else:
        P.axis('off')

    if ( visibleAxes ):
        if (nrows==2):
            fig.subplots_adjust(left=0.07,hspace=.07,bottom=.05,top=.95)
        else:
            fig.subplots_adjust(left=0.07,hspace=.07,bottom=.05,top=.9)
    else:
        fig.subplots_adjust(left=0.0,hspace=.0,top=1.,bottom=.0,right=1.,wspace=0.)
        
    pos = ax[0,1].get_position()
    lpixx = (L/(pos.x1-pos.x0))
    my_dpi = int(np.floor(lpixx/fw_inches))*pixsize
    
    fig.savefig(outfile,dpi=my_dpi)
    P.close()
