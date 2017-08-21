print("Importing")

# Import libraries to do our magic
import numpy as np
import h5py
import sys
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.colors as colors
import matplotlib.ticker as ticker

import matplotlib.pyplot as P

from sys import path
path.append("src/")
import tab_interp

# luminosity = 6.2849e42 # erg/s

def savephaseplots(run_id,output_dir,snap_str,includedVals):
    print("plotting:"+"".join([run_id,output_dir,snap_str,"".join(includedVals)]))
    f = h5py.File("/export/1/djw/gizmos/"+run_id+"/"+output_dir+"/snapshot_"+snap_str+".hdf5","r")

    header = f["/Header"]
    time = header.attrs.get("Time")
    time*= 0.9778e9 # to yr
    time/=1.e6 # to Myr

    requiredVals = list(includedVals)

    if ( "rho_p" not in requiredVals ):
        requiredVals.append("rho_p")
    if ( "dt_p" not in requiredVals ):
        requiredVals.append("dt_p")

    preqs = dict()
    preqs["TK_p"]=["u_p"]
    if ( "/PartType0/Pressure" not in f ):
        preqs["p_p"]=["nH_p","TK_p"]
    preqs["dustTemp"] = ["nH_p","TK_p","flux_p","depth_p"]
    preqs["dHeat"] = ["nH_p","TK_p","flux_p","depth_p"]
    preqs["dHeat"] = ["u_p","agn_heat_p"]
    preqs["cs_p"] = ["u_p"]
    preqs["mJ_p"] = ["cs_p"]
    preqs["prat"] = ["p_p","nH_p","TK_p","mJ_p","m_p"]
    preqs["hz_rat"] = ["h_p","z_p"]

    iVal = 0
    while ( iVal<len(requiredVals) ):
        reqVal = requiredVals[iVal]
        if reqVal in preqs:
            for preq in preqs[reqVal]:
                if ( preq not in requiredVals ):
                    requiredVals.append(preq)
        iVal+=1

    values = dict()

    xyz = np.array(f["/PartType0/Coordinates"])
    rad_p = np.sqrt(xyz[:,0]**2+xyz[:,1]**2+xyz[:,2]**2)

    if ( "rad2d_p" in requiredVals ):
        values["rad2d_p"] = np.sqrt(xyz[:,0]**2+xyz[:,1]**2)
        values["rad2d_p"]*=1.e3 # to pc

    vel_p = np.array(f["/PartType0/Velocities"])
    values["vrad"] = np.sum(vel_p*xyz,axis=1)/rad_p

    values["z_p"] = np.abs(xyz[:,2])
    values["z_p"]*=1.e3 # to pc

    if ( "u_p" in requiredVals ):
        values["u_p"] = np.array(f["/PartType0/InternalEnergy"])
        values["u_p"]*=1.e10 # to erg/g

    if ( "rho_p" in requiredVals ):
        values["rho_p"] = np.array(f["/PartType0/Density"])
        values["rho_p"]*=6.77e-22 # to g/cm**3

    if ( "agn_heat_p" in requiredVals ):
        values["agn_heat_p"] = np.array(f["/PartType0/AGNHeat"])
        values["agn_heat_p"]*=1e10/3.08568e+16# to erg/s/g

    if ( "dt_p" in requiredVals ):
        values["dt_p"] = np.array(f["/PartType0/TimeStep"])
        values["dt_p"]*=0.9778e9 # to yr

    #values["p_p"] = (5./3.-1.)*values["u_p"]*values["rho_p"]

    if ( "h_p" in requiredVals ):
        values["h_p"] = np.array(f["/PartType0/SmoothingLength"])
        values["h_p"] *= 1.e3 # to pc

    N_part = values["rho_p"].size

    # if "/PartType0/AGNColDens" in f:
    #     depth_p = np.array(f["/PartType0/AGNColDens"]) # surface density units
    #     depth_exists = True
    # else:
    #     depth_p = np.zeros(N_part)
    #     depth_exists = False

    if ( "flux_p" in requiredVals ):
        values["flux_p"] = np.array(f["/PartType0/AGNIntensity"]) # energy per surface area per time
        values["flux_p"]*=1.989e+53/(3.086e21)**2/(3.08568e+16)

    if ( "depth_p" in requiredVals ):
        values["depth_p"] = np.array(f["/PartType0/AGNColDens"]) # surface density units
        values["depth_p"]*=(1.989e+43/3.086e+21**2) # to g/cm**2

    if ( "radrad_p" in requiredVals ):
        radaccel_p = np.array(f["/PartType0/RadiativeAcceleration"])
        radaccel_p*=3.0857e21/3.08568e+16**2 # to cm/s/s

        values["radrad_p"] = (xyz[:,0]*radaccel_p[:,0]+xyz[:,1]*radaccel_p[:,1]+xyz[:,2]*radaccel_p[:,2])/rad_p

    molecular_mass = 4./(1.+3.*.76)
    proton_mass_cgs = 1.6726e-24
    gamma_minus_one = 5./3.-1.
    boltzmann_cgs = 1.38066e-16

    if ( "depth_p" in requiredVals ):
        values["depth_p"]/=(molecular_mass*proton_mass_cgs) # N in cm**(-2)
    # flux_p = luminosity/(4.*np.pi*rad_p**2)
    # flux_p/=(3.086e21)**2 # erg/s/kpc**2 to erg/s/cm**2

    if ( "nH_p" in requiredVals ):
        values["nH_p"] = values["rho_p"]/(molecular_mass*proton_mass_cgs)

    if ( "TK_p" in requiredVals ):
        values["TK_p"] = gamma_minus_one/boltzmann_cgs*(molecular_mass*proton_mass_cgs)*values["u_p"]


    if ( "p_p" in requiredVals ):
        if ( "/PartType0/Pressure" in f ):
            values["p_p"] = np.array(f["/PartType0/Pressure"])
            values["p_p"] *=(1.989e+43/3.086e+21/3.08568e+16**2) # to dyne/cm**2 = g cm s**-2 cm**-2 = g s**-2 cm**-1
        else:
            values["p_p"] = values["nH_p"]*values["TK_p"]*boltzmann_cgs

    if ( "dustTemp" in requiredVals or "dHeat" in requiredVals ):
        #print("Loading table (short form)")
        chTab = tab_interp.CoolHeatTab(("coolheat_tab_marta/shrunk_table_labels_130617.dat"),("coolheat_tab_marta/shrunk_table_130617.dat"))
        interpTabVec = np.vectorize(chTab.interpTab)
        tabStructs = interpTabVec(values["nH_p"].astype(np.float64),values["TK_p"].astype(np.float64),values["flux_p"].astype(np.float64),values["depth_p"].astype(np.float64))

    if ( "dustTemp" in requiredVals ):
        values["dustTemp"] = map(lambda y: y.dustT, tabStructs)
        values["dustTemp"] = np.array(values["dustTemp"])

    if ( "dHeat" in requiredVals ):
        values["dHeat"] = map(lambda y: y.dHeat, tabStructs)
        values["dHeat"] = np.array(values["dHeat"])

    #values["dt_heat"] = values["u_p"]/(10.**values["dHeat"])*values["rho_p"]/3.154e7

    if ( "dt_heat" in requiredVals ):
        values["dt_heat"] = values["u_p"]/values["agn_heat_p"]/3.154e7

    if ( "cs_p" in requiredVals ):
        values["cs_p"] = np.sqrt(values["u_p"])


    G = 6.67e-8 # in cgs
    G_Jeansterm = G**1.5 * 1.989e33

    if ( "mJ_p" in requiredVals ):
        values["mJ_p"] = 2./3.*(np.pi)**2.5 * values["cs_p"]**3/G_Jeansterm/np.sqrt(values["rho_p"])
    #mJ_p/=1.989e33 # g to msun

    if ( "m_p" in requiredVals ):
        values["m_p"] = np.array(f["/PartType0/Masses"])*1.e10
    # ratio between ideal gas pressure and actual pressure
    if ( "prat" in requiredVals ):
        values["prat"] = values["p_p"]/(values["nH_p"]*values["TK_p"]*boltzmann_cgs)*(values["mJ_p"]/values["m_p"])**(2./3.)

    # ratio between softening and z coordinate
    if ( "hz_rat" in requiredVals ):
        values["hz_rat"] = values["h_p"]/values["z_p"]


    #nH_p = 10.**((4.-np.log10(TK_p))**2)

    # values = [nH_p,flux_p,TK_p,depth_p]
    # titles = ["nH","flux","TK","col"]
    # labels = [r"$\log n_H$ (cm$^{-3}$)",r"$\log \Phi$ (erg/cm$^{-2}$/s)",r"$\log T$ (K)",r"$\log N_H$ (cm$^{-2}$)"]

    # oldvalues = [nH_p,flux_p,depth_p,TK_p,rad2d_p,z_p,dustTemp,radrad_p,vrad,rad_p,mJ_p,p_p,dt_p]
    #values = [nH_p,flux_p,depth_p,TK_p]
    #titles = ["nH","flux","col","Tg","r","z","Td","arad","vrad","rad","mJ","p","dt"]
    # labels = [r"$\log n_H$ (cm$^{-3}$)",r"$\log \Phi$ (erg/cm$^{-2}$/s)",r"$\log N_H$ (cm$^{-2}$)",r"$\log T_g$ (K)",r"$R$ (kpc)",r"$|z|$ (kpc)",r"$T_d$ (K)",r"$a_{rad}$ (cm/s/s)",r"v_{rad}$ km/s",r"$R$ (kpc)",r"$M_\mathrm{J}$ (M$_\odot$)",r"$P$ (dyne/cm$^{-2}$)",r"dt (yr)"]
    # dolog = [True,True,True,True,False,False,False,True,False,False,True,True,True]
    # ranges = [[0,12],None,None,[1,5],[0,.02],[0.,.02],[0,200],None,[-500,500],None,None,None,None]
    #ranges = [None,None,None,None,[0,.02],[0.,.02],None,None]
    #ranges = [None,None,None,None,None,None,None,None]
    #ranges = [None]*4

    labels = dict()
    dolog = dict()
    ranges = dict()

    labels['nH_p'] = r"$\log n_H$ (cm$^{-3}$)"
    dolog['nH_p'] = True
    ranges['nH_p'] = [0,12]

    labels['flux_p'] = r"$\log \Phi$ (erg/cm$^{-2}$/s)"
    dolog['flux_p'] = True
    ranges['flux_p'] = None

    labels['depth_p'] = r"$\log N_H$ (cm$^{-2}$)"
    dolog['depth_p'] = True
    ranges['depth_p'] = None

    labels['TK_p'] = r"$\log T_g$ (K)"
    dolog['TK_p'] = True
    ranges['TK_p'] = [1,5]

    labels['rad2d_p'] = r"$R$ (pc)"
    dolog['rad2d_p'] = False
    ranges['rad2d_p'] = [0,20.]

    labels['z_p'] = r"$|z|$ (pc)"
    dolog['z_p'] = False
    ranges['z_p'] = [0,20.]

    labels['dustTemp'] = r"$T_d$ (K)"
    dolog['dustTemp'] = False
    ranges['dustTemp'] = [0,200]

    labels['radrad_p'] = r"$a_{rad}$ (cm/s/s)"
    dolog['radrad_p'] = True
    ranges['radrad_p'] = None

    labels['dHeat'] = r"$H$"
    dolog['dHeat'] = False
    ranges['dHeat'] = None

    labels['dt_heat'] = r"dt$_H$"
    dolog['dt_heat'] = True
    ranges['dt_heat'] = None


    labels['vrad'] = r"v_{rad}$ km/s"
    dolog['vrad'] = False
    ranges['vrad'] = [-500.,500.]

    labels['rad_p'] = r"$R$ (pc)"
    dolog['rad_p'] = False
    ranges['rad_p'] = None

    labels['mJ_p'] = r"$M_\mathrm{J}$ (M$_\odot$)"
    dolog['mJ_p'] = True
    #ranges['mJ_p'] = None
    ranges['mJ_p'] = [-1.,5.]

    labels['p_p'] = r"$P$ (dyne/cm$^{-2}$)"
    dolog['p_p'] = True
    ranges['p_p'] = None

    labels['dt_p'] = r"dt (yr)"
    dolog['dt_p'] = True
    ranges['dt_p'] = None

    labels['cs_p'] = r"$c_s$ (cm/s)"
    dolog['cs_p'] = True
    ranges['cs_p'] = None

    labels['h_p'] = r"$h$ (pc)"
    dolog['h_p'] = True
    ranges['h_p'] = None

    #labels['prat'] = r"$\log_{10} p/p$"
    #dolog['prat'] = True
    #ranges['prat'] = [-.5,3.]
    labels['prat'] = r"$\log_{10} [(p_p/p_g)(M_p/M_J)^{2/3}]$"
    dolog['prat'] = True
    ranges['prat'] = [1.,3.]

    labels['hz_rat'] = r"$h/z$"
    dolog['hz_rat'] = False
    ranges['hz_rat'] = [0.,2.]

    #labels["nH_p"] = 

    #includedVals = [0,1,2,3,4,5,6,7]
    #includedVals = [0,3]
    #includedVals = [0,11]
    #includedVals = [0,1,2,3]


    # if ( not depth_exists ):
    #     nv = len(values)-1
    # else:
    #     nv = len(values)
    nv = len(values)

    #bigslice = (dustTemp>235.) & (dustTemp<245.) & (depth_p<1.e15)
    #bigslice = (TK_p>1.e6)
    #bigslice = (agn_heat_p<-1.e6)
    #bigslice = (dt_p<1.e-10)
    #bigslice = (TK_p>10.**3.5) & (TK_p<10.**4.) & (depth_p>10.**22.) & (depth_p<10.**23.2)
    #bigslice = (TK_p==10.)

    bigslice = (values["dt_p"]>0.)
    #bigslice = (np.log10(values["dt_p"])<1.8)

    noranges = False

    outfiles = []

    i = 0
    #for iv in ["dt_p"]:
    for i,iv in enumerate(includedVals):
        for jv in includedVals[i+1:]:
    #for iv in range(7):
    #    for jv in [7]:
            fname = "pics/"+run_id+output_dir+iv+jv+snap_str+".png"
            outfiles.append(fname)
            P.figure()
            P.xlabel(labels[iv])
            P.ylabel(labels[jv])
            if ( dolog[iv] ):
                vx = np.log10(values[iv])
            else:
                vx = values[iv]
            if ( dolog[jv] ):
                vy = np.log10(values[jv])
            else:
                vy = values[jv]
            notnans = ((np.isfinite(vx)) & (np.isfinite(vy))) & bigslice
            if ( not noranges and not (ranges[iv] is None) ):
                notnans = notnans & (vx>=ranges[iv][0]) & (vx<=ranges[iv][1])
            if ( not noranges and not (ranges[jv] is None) ):
                notnans = notnans & (vy>=ranges[jv][0]) & (vy<=ranges[jv][1])
            vx = vx[notnans]
            vy = vy[notnans]
            # CHECK FOR RANGES
            if ( not (ranges[iv] is None) and not (ranges[jv] is None) and not noranges ):
                H,xedges,yedges = np.histogram2d(vx,vy,bins=(150,150),range=[ranges[iv],ranges[jv]])
            else:
                H,xedges,yedges = np.histogram2d(vx,vy,bins=(150,150))
            with np.errstate(divide='ignore'):
                H = np.log10(H).T
            #H = H.T
            #P.pcolormesh(xedges,yedges,H,cmap='plasma',vmin=np.min(H[np.isfinite(H)]),vmax=np.max(H[np.isfinite(H)])) #,norm=colors.LogNorm()
            P.pcolormesh(xedges,yedges,H,cmap='viridis',vmin=0.,vmax=3.5) #,norm=colors.LogNorm()
            P.xlim(xedges[0],xedges[-1])
            P.ylim(yedges[0],yedges[-1])
            ax = P.gca()
            #ax.xaxis.set_minor_locator(ticker.MultipleLocator(1.))
            #ax.yaxis.set_minor_locator(ticker.MultipleLocator(1.))
            P.tick_params(which='both',direction="out")
            P.colorbar(label=r"$\log$ count")
            #P.colorbar(label=r"count")
            P.suptitle(r"$t="+("%.4f" % time)+"$ Myr")
            P.savefig(fname,dpi=150)
            P.close()
    return outfiles


if __name__ == '__main__':
    run_id = sys.argv[1]
    output_dir = sys.argv[2]
    snap_str = sys.argv[3]

    print("Running")


    #includedVals = ["dt_p","nH_p","TK_p","rad2d_p","z_p","vrad","dHeat","dt_heat"]
    #includedVals = ["nH_p","TK_p","rad2d_p","z_p","vrad"]
    #includedVals = ["mJ_p","nH_p","TK_p","p_p"]
    #includedVals = ["mJ_p","TK_p"]
    #includedVals = ["rad2d_p","hz_rat"]
    includedVals = ["mJ_p","nH_p","TK_p"]
    #includedVals = ["prat"]
    #includedVals = ["nH_p","p_p"]
    #includedVals = ["h_p","nH_p"]
    #includedVals = ["mJ_p","prat"]
    #includedVals = ["cs_p","h_p"]

    x = savephaseplots(run_id,output_dir,snap_str,includedVals)
    print(x)

