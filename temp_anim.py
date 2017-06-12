print("Importing")

# Import libraries to do our magic
import numpy as np
import h5py
import sys
import matplotlib as mpl
mpl.use('Agg')

import pylab as P

import os
from joblib import Parallel, delayed

print("Running")


run_id = sys.argv[1]
output_dir = sys.argv[2]

#snapi = int(sys.argv[3])
#snapf = int(sys.argv[4])

snapi = 0
fullDir = "/export/1/djw/gizmos/"+run_id+"/"+output_dir

fnames = os.listdir(fullDir)
fnames = np.array(fnames)
fnames.sort()
snapshotfilebools = np.array([x.startswith("snapshot") for x in fnames])
snapshotfiles = fnames[snapshotfilebools]
snapf = 0
ctime = os.path.getmtime(fullDir+"/snapshot_000.hdf5")
for fname in snapshotfiles[1:]:
    new_snapf = int(fname[9:12])
    new_ctime = os.path.getmtime(fullDir+"/"+fname)
    if ( new_ctime>ctime ) :
        ctime = new_ctime
        snapf = new_snapf

#snapf = 8

print(snapi,snapf)
# if ( len(sys.argv)>3 ):
#     snapf = int(sys.argv[4])
# else:
#     snapf = snapi

animfiles = ""

os.system("rm pics/plot???.png")

def makeplot(snapx):
    print(snapx)

    #f = h5py.File("/export/1/djw/gizmo_public/disc_sf_out/snapshot_"+("%03d" % snapx)+".hdf5","r")
    f = h5py.File(fullDir+"/snapshot_"+("%03d" % snapx)+".hdf5","r")
    
    header = f["/Header"]
    time = header.attrs.get("Time")
    time*= 0.9778e9 # to yr
    time/=1.e6 # to Myr

    xyz = np.array(f["/PartType0/Coordinates"])

    u_p = np.array(f["/PartType0/InternalEnergy"])
    rho_p = np.array(f["/PartType0/Density"])

    u_p*=1.e10 # to erg/g
    rho_p*=6.77e-22 # to g/cm**3

    molecular_mass = 4./(1.+3.*.76)
    proton_mass_cgs = 1.6726e-24
    gamma_minus_one = 5./3.-1.
    boltzmann_cgs = 1.38066e-16

    TK_p = np.log10(gamma_minus_one/boltzmann_cgs*(molecular_mass*proton_mass_cgs)*u_p)
    nH_p = np.log10(rho_p/(molecular_mass*proton_mass_cgs))

    P.close('all')

    #fig = P.figure()
    fig, ax = P.subplots(2,4, gridspec_kw = {'width_ratios':[1, 16,16,1]})
    fig.suptitle(r"$T="+("%.4f" % time)+"$ Myr")
    fw_inches = 10.
    fig.set_figwidth(fw_inches)
    fig.set_figheight(8.)
    

    lres = 300#384#256
    plotrad = 1.e-2
    #cent = 0.

    x = xyz[:,0]
    y = xyz[:,1]
    z = xyz[:,2]
    
    r2d = np.sqrt(x**2+y**2)
    #theta = np.atan2(y,x)
#    z = TK_p
    
    makehist(fig,ax[0,1],ax[0,0],x,y,TK_p,lres,plotrad,[0.,0.],r"$\log_{10} T$ (K)",0.,6.,"cool")
    makehist(fig,ax[0,2],ax[0,3],x,y,nH_p,lres,plotrad,[0.,0.],r"$\log_{10} n_\mathrm{H}$ (cm$^{-3}$)",2.,7.,"inferno")
    makehist(fig,ax[1,1],ax[1,0],r2d,z,TK_p,lres,plotrad,[plotrad,0.],r"$\log_{10} T$ (K)",0.,6.,"cool")
    makehist(fig,ax[1,2],ax[1,3],r2d,z,nH_p,lres,plotrad,[plotrad,0.],r"$\log_{10} n_\mathrm{H}$ (cm$^{-3}$)",2.,7.,"inferno")
    for this_ax in (ax[0,1],ax[1,1]):
        this_ax.yaxis.tick_right()
        this_ax.yaxis.set_visible(False)
#        this_ax.yaxis.set_label_position("right")

    
    for this_ax in (ax[0,0],ax[1,0]):
        this_ax.yaxis.tick_left()
        this_ax.yaxis.set_label_position("left")


    outfile = "pics/plot"+("%03d" % snapx)+".png"
    #animfiles += outfile +" " 
    fig.subplots_adjust(left=0.07,hspace=.07,bottom=.05,top=.95)
    #print(ax[0,1].get_position())
    #print(ax[1,2].get_position())
    # get 1:1 pixel size
    pos = ax[0,1].get_position()
    lpixx = (lres/(pos.x1-pos.x0))
    my_dpi = int(np.floor(lpixx/fw_inches))
    #for this_ax in (ax[0,1],ax[0,2],ax[1,1],ax[1,2]):
    #    pos = this_ax.get_position()
    #    lpixx = (lres/(pos.x1-pos.x0))
    #    lpixy = (lres/(pos.y1-pos.y0))
    #    print(pos)
    #    print(lpixx,lpixy,lpixy/.8)

    fig.savefig(outfile,dpi=my_dpi)
    return outfile

def makehist(fig,sp,cbax,x,y,z,lres,plotrad,cent,cblabel,zlow,zhigh,cmap_label):

    this_cmap = P.get_cmap(cmap_label)
    this_cmap.set_bad('black',1.)
    this_cmap.set_under('black',1.)


    lx = np.floor((x+plotrad-cent[0])/(plotrad*2.)*lres).astype(int)
    ly = np.floor((y+plotrad-cent[1])/(plotrad*2.)*lres).astype(int)
    plotgrid = np.zeros((lres,lres))
    plotgrid[:,:] = -1.

    ingrid = (lx>=0) & (lx<lres) & (ly>=0) & (ly<lres)
    step = plotrad/lres
    xedges = np.arange(-plotrad+cent[0],plotrad+cent[0],plotrad*2./lres)
    yedges = np.arange(-plotrad+cent[1],plotrad+cent[1],plotrad*2./lres)
    xedges*=1.e3 # convert from kpc to pc
    yedges*=1.e3 # convert from kpc to pc

    #plotgrid[:,:] = -1
    plotgrid[lx[ingrid],ly[ingrid]] = z[ingrid]
    #plotgrid = np.log10(plotgrid)
    mesh = sp.pcolormesh(xedges,yedges,plotgrid.T,cmap=this_cmap,vmin=zlow,vmax=zhigh) 
    #sp.set_xlabel("x (kpc)")
    #sp.set_ylabel("y (kpc)")
    sp.axis('equal')
    #x=sp.scatter(xyz[:,0],xyz[:,1],c=TK_p,marker=",",edgecolor='',vmin=1.,vmax=8.,s=1.)
    cb = fig.colorbar(mesh,label=cblabel,cax=cbax)
    #sp.set_xlim([-1.e-2,1.e-2])
    #sp.set_ylim([-1.e-2,1.e-2])
    #sp.set_cbrange([1.,8.]) 

#for snapx in range(snapi,snapf+1):
#    print(snapx)
#    animfiles+=makeplot(snapx)+" "

animfiles = Parallel(n_jobs=32)(delayed(makeplot)(snapx) for snapx in range(snapi,snapf+1))
#animfiles = " ".join(animfiles)
#print(animfiles)
#print("to gif!")
#cmd = "convert "+animfiles+" pics/plot.gif"
#system(cmd)
print("to mp4!")
cmd = "ffmpeg -y -r 24 -i pics/plot%03d.png -c:v mpeg4 -q:v 1 /export/1/djw/movies/rhotempgiz_"+run_id+"_"+output_dir+".mp4"
os.system(cmd)
