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




densslice = 0
weightslice = 1
zdensslice = 2
zweightslice = 3
vec2dslice = 4
viewslice = 5

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

def makesph_plot(fig,sp,cbax,x_p,y_p,z_p,zslice,val_p,m_p,h_p,L,mask,corner,width,cblabel,clow,chigh,cmap_label,mode,dolog,**kwargs):
    if "planenorm" in kwargs:
        planenorm = kwargs["planenorm"]
    else:
        planenorm = False

    if "circnorm" in kwargs:
        circnorm = kwargs["circnorm"]
    else:
        circnorm = False

    if "plusminus" in kwargs:
        plusminus = kwargs["plusminus"]
        cbax2 = kwargs["cbar2"]
    else:
        plusminus = False
    

    this_cmap = P.get_cmap(cmap_label)
    if ( plusminus ):
        cmap_label2 = kwargs["cmap2"]
        this_cmap2 = P.get_cmap(cmap_label2)
    
    if ( not plusminus ):
        this_cmap.set_bad('black',1.)
        this_cmap.set_under('black',1.)
#         if ( not ("gray" in cmap_label or "Grey" in cmap_label) ):
#             this_cmap.set_bad('black',1.)
#             this_cmap.set_under('black',1.)
#         else:
#             this_cmap.set_bad('yellow',1.)
#             this_cmap.set_under('yellow',1.)

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
        qv = sp.quiver(xmids,ymids,map1,map2,norm_map,headwidth=10.,pivot='mid',cmap=this_cmap,clim=[clow,chigh])
        cb = fig.colorbar(qv,label=cblabel,cax=cbax)
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
        if ( circnorm ):
            azimuthalNorm(map)

        if ( plusminus ):
        
            plusmap = map.copy()
            minusmap = map.copy()
            
            isplus = (map>0.)
            isminus = (map<0.)
            
            if ( np.sum(isplus)>0 ):
                plusmap[isminus]=0.
                mesh1=sp.pcolormesh(xedges,yedges,plusmap.T,cmap=this_cmap,vmin=clow,vmax=chigh,norm = colors.LogNorm())
                cb = fig.colorbar(mesh1,label=cblabel,cax=cbax)
            
            if ( np.sum(isminus)>0 ):
                minusmap[isplus]=0.
                minusmap = -minusmap
            
                mesh2=sp.pcolormesh(xedges,yedges,minusmap.T,cmap=this_cmap2,vmin=clow,vmax=chigh,norm = colors.LogNorm())
                cb2 = fig.colorbar(mesh2,label=cblabel,cax=cbax2)
        else:

            if ( dolog ):
                map = np.log10(map)
            print(np.min(map),np.max(map))

            mesh = sp.pcolormesh(xedges,yedges,map.T,cmap=this_cmap,vmin=clow,vmax=chigh)
            cb = fig.colorbar(mesh,label=cblabel,cax=cbax)
    
    #sp.set_xlim([corner[0],corner[0]+width])
    #sp.set_ylim([corner[1],corner[1]+width])
    sp.plot([0],[0],'+g',markersize=10.,markeredgewidth=1.)
    #mesh = sp.pcolormesh(xedges,yedges,map.T,cmap=this_cmap) 
    sp.axis('equal')

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
    if ( "vlos" in need_to_load ):
        need_to_load.append("vels")
    if ( "emit" in need_to_load ):
        need_to_load.append("tdust")
    if ( "tdust" in need_to_load ):
        need_to_load.append("col")
        need_to_load.append("temp")

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
    if ( "heat" in need_to_load ):
        data.heat = np.array(f["/PartType0/AGNHeat"])
        data.heat*=1e10/3.08568e+16# to erg/s/g
        
#         heatplus = (data.heat>0.)
#         heatminus = (data.heat<0.)
#         
#         data.heat[heatplus] = np.maximum(0.,np.log10(data.heat[heatplus]))
#         data.heat[heatminus] = -np.maximum(0.,np.log10(-data.heat[heatminus]))

    if ( "depth" in need_to_load):
        #rho_p = np.array(f["/PartType0/Density"])
        #agn_heat_p = np.array(f["/PartType0/AGNHeat"])
        data.depth_p = np.array(f["/PartType0/AGNDepth"]) # unitless

    if ( "tdust" in need_to_load ):
        if ( not chTab ):
            print("Load dust tables")
            chTab = tab_interp.CoolHeatTab(("../coolheat_tab_marta/shrunk_table_labels_170517.dat"),("../coolheat_tab_marta/shrunk_table_170517.dat"))
            interpTabVec = np.vectorize(chTab.interpTab)
#         
#         if ( not hasattr(data,'coldens') ):
#             data.coldens = np.array(f["/PartType0/AGNColDens"]) # Msun/kpc**2
#             data.coldens*=(1.989e+43/3.086e+21**2) # to g/cm**2
#             data.coldens/=(molecular_mass*proton_mass_cgs) # N in cm**(-2) 
#         
#         if ( not hasattr(data,'u_p') ):
#             data.u_p = np.array(f["/PartType0/InternalEnergy"]) # 1e10 erg/g
#             data.u_p*=1.e10 # to erg/g
#             data.TK_p = (gamma_minus_one/boltzmann_cgs*(molecular_mass*proton_mass_cgs)*data.u_p)
#         
        data.flux_p = np.array(f["/PartType0/AGNIntensity"]) # energy per surface area per time
        data.flux_p*=1.989e+53/(3.086e21)**2/(3.08568e+16)
        
        data.rho_p = np.array(f["/PartType0/Density"])
        data.rho_p*=6.77e-22 # to g/cm**3 
        data.nH_p = data.rho_p/(molecular_mass*proton_mass_cgs)
        
        print("Calculating dust temperatures from table")
        tabStructs = interpTabVec(data.nH_p.astype(np.float64),data.TK_p.astype(np.float64),data.flux_p.astype(np.float64),data.coldens.astype(np.float64))
        data.dustTemp = map(lambda y: y.dustT, tabStructs)
        data.dustTemp = np.array(data.dustTemp)        
    
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

def makesph_trhoz_frame(infile,outfile,**kwargs):
    if "cmap" in kwargs:
        cmap = kwargs["cmap"]
    else:
        cmap = "viridis"
        
    cmap_r = cmap+"_r"
    
    cmap_d = "coolwarm"

    if "L" in kwargs:
        L = kwargs["L"]
    else:
        L = 256

    if "ring" in kwargs:
        ringPlot = kwargs["ring"]
    else:
        ringPlot = False


    if "flat" in kwargs:
        flatPlot = kwargs["flat"]
    else:
        flatPlot = False

    if "subsample" in kwargs:
        subsample = kwargs["subsample"]
        if ( L%subsample!=0 ):
            raise Exception("subsample might divide evenly into L")
    else:
        subsample = 1

    if "pixsize" in kwargs:
        pixsize = kwargs["pixsize"]
    else:
        pixsize = subsample


    if "plot" in kwargs:
        plot_thing = kwargs["plot"]
    else:
        plot_thing = ['dens','temp','depth']

    if "scale" in kwargs:
        width = kwargs["scale"]*2.
    else:
        width = 40.

    if "planenorm" in kwargs:
        planenorm = kwargs["planenorm"]
    else:
        planenorm = False

    if "cols" in kwargs:
        cols = kwargs["cols"]
        if ( cols!=1 and cols!=2 ):
            raise Exception("cols=1 or cols=2")
    else:
        cols = 2

    if "rot" in kwargs:
        rot = kwargs["rot"]
        if (len(rot)!=2):
            raise Exception("rot needs to be [theta,phi]")
    else:
        rot = [0.,0.]

    nrows = len(plot_thing)


    print("Loading",infile)
    time,data = load_gadget(infile,plot_thing)
    print("Plotting",infile)



    n = data.h_p.size
#    L = 256
    
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
    
    
    if ( "vels" in plot_thing ):
        #vel_mag = np.sqrt(np.sum(data.vels[:,:]**2,1))
        vel2d = (x*data.vels[:,0]+y*data.vels[:,1])/np.sqrt(x**2+y**2)
    
#    r2d = np.sqrt(x**2+y**2)
#    r3d_kpc = np.sqrt(x**2+y**2+z**2)/1.e3
#     if ( 'agn_heat_p' in globals() ):
#         depth_p *= r3d_kpc
#         depth_p /= (1.e7/1683) # take out luminosity (from solar L to internal units)
#         depth_p = -np.log(depth_p) # to optical depth
    #agn_heat_p*=rho_p*h_p*r3d*(16./3./np.pi) # to L*exp(-optical depth)
    
    # mask out non-gas - currently everything is gas
    mask = np.full(n,True,dtype=bool)
    
    # flat weighting - dummy value required for some functions because I'm not qwarging properly yet
    n_ones = np.ones(n)

    #corner = np.array([-.01,-.01])

    # figure properties
    if ( "heat" in plot_thing ):
        fig, ax = P.subplots(nrows,3*cols, gridspec_kw = {'width_ratios':([1, 1, 16,16,1, 1])[0:3*cols]})
        ax = np.resize(ax, (3,6)) # because 1D arrays have different syntax, we have to pretend it's 3D
        cbax2left_index = 0
        cbaxleft_index = 1
        spleft_index = 2
        spright_index = 3
        cbaxright_index = 4
        cbax2right_index = 5
        
    else:
        fig, ax = P.subplots(nrows,2*cols, gridspec_kw = {'width_ratios':([1, 16,16,1])[0:2*cols]})
        ax = np.resize(ax, (3,4)) # because 1D arrays have different syntax, we have to pretend it's 3D
        cbaxleft_index = 0
        spleft_index = 1
        spright_index = 2
        cbaxright_index = 3
        
    fig.suptitle(r"$T="+("%.4f" % time)+"$ Myr")
    fw_inches = 5.*cols
    fig.set_figwidth(fw_inches)
    fig.set_figheight(4.*nrows)

    #corner = np.array([0.,-.01])
    # physical coordinates of region to plot, in pc
#     corners = [-10.,-10.]
#     width = 20.
    corners = [-width/2.,-width/2.]
    #map = sph_plotter.sph_dense(r2d,z,m_p,h_p,L,corner,.02,mask,n)
    #map = sph_plotter.sph_weight(r2d,z,m_p,h_p,TK_p,L,corner,width,mask,n)

    if ( ringPlot ):
        if ( not flatPlot ):
            raise Exception("ring==true requires flat==true")
        corners_side = [0.,-width/2.]
        rad2d = np.sqrt(x**2+y**2)
    else:   
        corners_side = corners
        rad2d = x
#        deep_side = y

    deep_face = z
    deep_side = y
    if ( flatPlot ):
#         deep_face = np.zeros(n)
#         deep_side = np.zeros(n)
        quantslice = weightslice
        dslice = densslice
        rhounit = r"$\log_{10} \Sigma$ (M$_\odot$/pc$^2$)"
        drange = [-2.,8.]
        #drange = [2.,6.]
        #drange = [-1.,2.]
    else:
        quantslice = zweightslice
        dslice = zdensslice
        rhounit = r"$\log_{10} \rho$ (M$_\odot$/pc$^3$)"
        drange = [-2.,7.]
    
    if ( planenorm ):
        drange_s = [-3.6,0.]
        drange = [-4.5,0.]
    else:
        drange_s = drange
    
    # do all subplots, calculating the full SPH smoothing each time
    
    for irow in range(nrows):
        if ( plot_thing[irow]=='temp' ):
            makesph_plot(fig,ax[irow,spleft_index],ax[irow,cbaxleft_index],x,y,deep_face,0.,data.TK_p,data.m_p,data.h_p,L,mask,corners,width,r"$\log_{10} T_g$ (K)",1.,6.,cmap,quantslice,True)
            if ( cols==2 ):
                makesph_plot(fig,ax[irow,spright_index],ax[irow,cbaxright_index],rad2d,z,deep_side,0.,data.TK_p,data.m_p,data.h_p,L,mask,corners_side,width,r"$\log_{10} T_g$ (K)",1.,6.,cmap,quantslice,True)
            #makesph_plot(fig,ax[0,2],ax[0,3],x,z,y,0.,TK_p,m_p,h_p,L,mask,corners,width,r"$\log_{10} T$ (K)",1.,6.,cmap,zweightslice,True)
            #makesph_plot(fig,ax[1,1],ax[1,0],r2d,z,TK_p,m_p,h_p,L,mask,[0.,-10.],width,r"$\log_{10} T$ (K)",0.,6.,"plasma",True,True)
        elif ( plot_thing[irow]=='col' ):
            makesph_plot(fig,ax[irow,spleft_index],ax[irow,cbaxleft_index],x,y,deep_face,0.,data.coldens,data.m_p,data.h_p,L,mask,corners,width,r"$\log_{10} N$ (cm$^{-2}$)",15.,26.,cmap,quantslice,True)
            if ( cols==2 ):
                makesph_plot(fig,ax[irow,spright_index],ax[irow,cbaxright_index],rad2d,z,deep_side,0.,data.coldens,data.m_p,data.h_p,L,mask,corners_side,width,r"$\log_{10} N$ (cm$^{-2}$)",15.,26.,cmap,quantslice,True)
        elif ( plot_thing[irow]=='heat' ):
            makesph_plot(fig,ax[irow,spleft_index],ax[irow,cbaxleft_index],x,y,deep_face,0.,data.heat,data.m_p,data.h_p,L,mask,corners,width,r"$H$",.1e-2,1.e10,'Reds',quantslice,False,cmap2='Blues',plusminus=True,cbar2=ax[irow,cbax2left_index])
            if ( cols==2 ):
                makesph_plot(fig,ax[irow,spright_index],ax[irow,cbaxright_index],rad2d,z,deep_side,0.,data.heat,data.m_p,data.h_p,L,mask,corners_side,width,r"$H$",.1e-2,1.e10,'Reds',quantslice,False,cmap2='Blues',plusminus=True,cbar2=ax[irow,cbax2right_index])
        elif ( plot_thing[irow]=='dens' ):
            makesph_plot(fig,ax[irow,spleft_index],ax[irow,cbaxleft_index],x,y,deep_face,0.,None,data.m_p,data.h_p,L,mask,corners,width,rhounit,drange[0],drange[1],cmap,dslice,True,circnorm=planenorm)
            if ( cols==2 ):
                makesph_plot(fig,ax[irow,spright_index],ax[irow,cbaxright_index],rad2d,z,deep_side,0.,None,data.m_p,data.h_p,L,mask,corners_side,width,rhounit,drange_s[0],drange_s[1],cmap,dslice,True,planenorm=planenorm)
            #    makesph_plot(fig,ax[1,2],ax[1,3],r2d,z,None,m_p,h_p,L,mask,[0.,-10.],width,r"$\log_{10} \Sigma$ (M$_\odot$/kpc$^2$)",-3.,3.,"plasma",False,True)
        elif ( plot_thing[irow]=='depth' ):
            #makesph_plot(fig,ax[2,1],ax[2,0],x,y,z,0.,depth_p,m_p,h_p,L,mask,corners,width,r"$\tau$",0.,9.,cmap,zweightslice,False)
            #makesph_plot(fig,ax[2,2],ax[2,3],x,z,y,0.,depth_p,m_p,h_p,L,mask,corners,width,r"$\tau$",0.,9.,cmap,zweightslice,False)
            #makesph_plot(fig,ax[2,1],ax[2,0],x,y,deep_face,0.,np.log10(depth_p),m_p,h_p,L,mask,corners,width,r"$\log_{10} \tau$",-2.,7.,cmap_r,zweightslice,False)
            #makesph_plot(fig,ax[2,2],ax[2,3],rad2d,z,deep_side,0.,np.log10(depth_p),m_p,h_p,L,mask,corners_side,width,r"$\log_{10} \tau$",-2.,7.,cmap_r,zweightslice,False)
            makesph_plot(fig,ax[irow,spleft_index],ax[irow,cbaxleft_index],x,y,deep_face,0.,data.depth_p,data.m_p,data.h_p,L,mask,corners,width,r"$\tau$",0.,1.e3,cmap_r,zweightslice,False)
            if ( cols==2 ):
                makesph_plot(fig,ax[irow,spright_index],ax[irow,cbaxright_index],rad2d,z,deep_side,0.,data.depth_p,data.m_p,data.h_p,L,mask,corners_side,width,r"$\tau$",0.,1.e3,cmap_r,zweightslice,False)
        elif ( plot_thing[irow]=='vels' ):
            makesph_plot(fig,ax[irow,spleft_index],ax[irow,cbaxleft_index],x,y,deep_face,0.,[data.vels[:,0],data.vels[:,1]],data.m_p,data.h_p,L//subsample,mask,corners,width,r"$\log_{10}v$ (km/s)",-1.,3.,cmap,vec2dslice,False)
            if ( cols==2 ):
                makesph_plot(fig,ax[irow,spright_index],ax[irow,cbaxright_index],rad2d,z,deep_side,0.,[vel2d,data.vels[:,2]],data.m_p,data.h_p,L//subsample,mask,corners_side,width,r"$\log_{10}v$ (km/s)",-1.,3.,cmap,vec2dslice,False)
        elif ( plot_thing[irow]=='tdust' ):
            makesph_plot(fig,ax[irow,spleft_index],ax[irow,cbaxleft_index],x,y,deep_face,0.,data.dustTemp,data.m_p,data.h_p,L,mask,corners,width,r"$\log_{10} T_d$ (K)",1.,3.,cmap,quantslice,True)
            if ( cols==2 ):
                makesph_plot(fig,ax[irow,spright_index],ax[irow,cbaxright_index],rad2d,z,deep_side,0.,data.dustTemp,data.m_p,data.h_p,L,mask,corners_side,width,r"$\log_{10} T_d$ (K)",1.,3.,cmap,quantslice,True)
        elif ( plot_thing[irow]=='view' ):
            makesph_plot(fig,ax[irow,spleft_index],ax[irow,cbaxleft_index],x,y,z,0.,[data.brightness,data.opac],data.m_p,data.h_p,L,mask,corners,width,r"$\log_{10} F$",0.,4.,cmap,viewslice,True)
            if ( cols==2 ):
                makesph_plot(fig,ax[irow,spright_index],ax[irow,cbaxright_index],y,z,x,0.,[data.brightness,data.opac],data.m_p,data.h_p,L,mask,corners_side,width,r"$\log_{10} F$",0.,4.,cmap,viewslice,True)
        elif ( plot_thing[irow]=='vlos' ):
            makesph_plot(fig,ax[irow,spleft_index],ax[irow,cbaxleft_index],x,y,z,0.,[data.vels[:,2],data.opac],data.m_p,data.h_p,L,mask,corners,width,r"$v_{LOS}$",-1.2e2,1.2e2,cmap,viewslice,False)
            if ( cols==2 ):
                makesph_plot(fig,ax[irow,spright_index],ax[irow,cbaxright_index],y,z,x,0.,[data.vels[:,0],data.opac],data.m_p,data.h_p,L,mask,corners_side,width,r"$v_{LOS}$",-1.2e2,1.2e2,cmap,viewslice,False)
        elif ( plot_thing[irow]=='emit' ):
            makesph_plot(fig,ax[irow,spleft_index],ax[irow,cbaxleft_index],x,y,deep_face,0.,None,data.emissivity,data.h_p,L,mask,corners,width,r"$\log_{10} F$",-12.,-3.,cmap,dslice,True,circnorm=planenorm)
            if ( cols==2 ):
                makesph_plot(fig,ax[irow,spright_index],ax[irow,cbaxright_index],rad2d,z,deep_side,0.,None,data.emissivity,data.h_p,L,mask,corners_side,width,r"$\log_{10} F$",-12.,-3.,cmap,dslice,True,circnorm=planenorm)
        else:
            raise Exception(plot_thing[irow]+" not a valid plot type")
        
        if ( "heat" in plot_thing and plot_thing[irow]!='heat' ):
            for visax in [ax[irow,cbax2left_index],ax[irow,cbax2right_index]]:
                visax.set_frame_on(False)
                #visax.get_xaxis().tick_bottom()
                visax.axes.get_yaxis().set_visible(False)
                visax.axes.get_xaxis().set_visible(False)
    
    #makesph_plot(fig,ax[0,2],ax[0,3],x,y,agn_heat_p,m_p,h_p,L,mask,[-10.,-10.],width,r"$\log_{10} du/dt$ (erg/g/s)",-3.,3.,"plasma",True)
    #makesph_plot(fig,ax[1,2],ax[1,3],r2d,z,agn_heat_p,m_p,h_p,L,mask,[0.,-10.],width,r"$\log_{10} du/dt$ (erg/g/s)",-3.,3.,"plasma",True)
#     makesph_plot(fig,ax[0,2],ax[0,3],x,y,depth_p,m_p,h_p,L,mask,[-10.,-10.],width,r"$\tau$",0.,5.,"viridis",True,False)
#     makesph_plot(fig,ax[1,2],ax[1,3],r2d,z,depth_p,m_p,h_p,L,mask,[0.,-10.],width,r"$\tau$ (erg/g/s)",0.,5.,"viridis",True,False)

    #for this_ax in (ax[0,1],ax[1,1],ax[2,1]):
    for iax in range(nrows):
        this_ax = ax[iax,1]
        this_ax.yaxis.tick_right()
        this_ax.yaxis.set_visible(False)

    
#    for this_ax in (ax[0,0],ax[1,0],ax[2,0]):
    for iax in range(nrows):
        this_ax = ax[iax,0]
        this_ax.yaxis.tick_left()
        this_ax.yaxis.set_label_position("left")


    fig.subplots_adjust(left=0.07,hspace=.07,bottom=.05,top=.95)
    pos = ax[0,1].get_position()
    lpixx = (L/(pos.x1-pos.x0))
    my_dpi = int(np.floor(lpixx/fw_inches))*pixsize
    
    fig.savefig(outfile,dpi=my_dpi)
    P.close()
