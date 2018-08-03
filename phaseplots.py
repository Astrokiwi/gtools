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

import gizmo_tools

# luminosity = 6.2849e42 # erg/s

G_kms_pc_msun = 0.0043022682

labels = dict()
dolog = dict()
ranges = dict()

labels['nH_p'] = r"$\log n_H$ (cm$^{-3}$)"
dolog['nH_p'] = True
ranges['nH_p'] = [-4,12]

labels['flux_p'] = r"$\log \Phi$ (erg/cm$^{-2}$/s)"
dolog['flux_p'] = True
#ranges['flux_p'] = [3.,8.]
ranges['flux_p'] = [0.,8.]

labels['depth_p'] = r"$\log N_H$ (cm$^{-2}$)"
dolog['depth_p'] = True
ranges['depth_p'] = [16,27]

labels['TK_p'] = r"$\log T_g$ (K)"
# dolog['TK_p'] = True
# ranges['TK_p'] = [1,6]
# labels['TK_p'] = r"$T_g$ (K)"
# dolog['TK_p'] = False
# ranges['TK_p'] = [0,1e4]
dolog['TK_p'] = True
ranges['TK_p'] = [1.,7.]

labels['rad2d_p'] = r"$R$ (pc)"
dolog['rad2d_p'] = False
#ranges['rad2d_p'] = [0,150.]
ranges['rad2d_p'] = None

labels['z_p'] = r"$|z|$ (pc)"
dolog['z_p'] = False
ranges['z_p'] = [0,20.]

labels['dustTemp'] = r"$T_d$ (K)"
dolog['dustTemp'] = False
ranges['dustTemp'] = [0,150]

labels['radrad_p'] = r"$a_{rad,r}$ (cm/s/s)"
dolog['radrad_p'] = False
ranges['radrad_p'] = None

labels['arad_p'] = r"$a_{all,r}$ (cm/s/s)"
dolog['arad_p'] = False
ranges['arad_p'] = None

labels['dHeat'] = r"$H$"
dolog['dHeat'] = False
ranges['dHeat'] = None

labels['dt_heat'] = r"dt$_H$"
dolog['dt_heat'] = True
ranges['dt_heat'] = None


labels['vrad'] = r"$v_{rad}$ (km/s)"
dolog['vrad'] = False
#ranges['vrad'] = [-300,300]
#ranges['vrad'] = [-300,1.e3]
# ranges['vrad'] = [-50,300.]
ranges['vrad'] = [-50,600.]


labels['vcirc'] = r"$v_{circ}$ km/s"
dolog['vcirc'] = False
ranges['vcirc'] = [-20.,150.]

labels['vel'] = r"$v$ km/s"
dolog['vel'] = False
ranges['vel'] = [0.,300.]

# labels['rad_p'] = r"$\log R$ (pc)"
# dolog['rad_p'] = True
# ranges['rad_p'] = [-1.,1.9]

labels['rad_p'] = r"$R$ (pc)"
dolog['rad_p'] = False
ranges['rad_p'] = [0.,5.]

labels['mJ_p'] = r"$M_\mathrm{J}$ (M$_\odot$)"
dolog['mJ_p'] = True
#ranges['mJ_p'] = None
ranges['mJ_p'] = [-2.,8.]

labels['p_p'] = r"$P$ (dyne/cm$^{-2}$)"
dolog['p_p'] = True
ranges['p_p'] = None

labels['dt_p'] = r"dt (yr)"
dolog['dt_p'] = True
# ranges['dt_p'] = [-1.5,3.]
ranges['dt_p'] = [0.,3.]

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
#ranges['prat'] = [1.,3.]
ranges['prat'] = None

labels['hz_rat'] = r"$h/z$"
dolog['hz_rat'] = False
ranges['hz_rat'] = [0.,2.]

labels['tsf'] = r"$t_{sf}$"
dolog['tsf'] = True
ranges['tsf'] = [3.,8.]

labels['agn_heat_p'] = r"$H$ (erg/s/g)"
dolog['agn_heat_p'] = False
#ranges['agn_heat_p'] = None
ranges['agn_heat_p'] = [-300.,300.]

labels['cool_p'] = r"$H$ (erg/s/g)"
dolog['cool_p'] = True
ranges['cool_p'] = None

labels['opac'] = r"$\kappa$ (pc$^2$/M$_\odot$)"
dolog['opac'] = True
ranges['opac'] = [-4.,-1.]

def v_esc(r,m_bh,m_hern,a_hern):
    return np.sqrt(2*G_kms_pc_msun*(m_bh/r + m_hern/(a_hern+r)))

def loadvalues(run_id,output_dir,snap_str,includedVals,rcut=None):
    gizmoDir = gizmo_tools.getGizmoDir(run_id)
    fullDir = gizmoDir+"/"+run_id+"/"+output_dir
    f = h5py.File(fullDir+"/snapshot_"+snap_str+".hdf5","r")

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
    preqs["dustTemp"] = ["nH_p","TK_p","flux_p","tau"]
    preqs["dHeat"] = ["nH_p","TK_p","flux_p","tau"]
    preqs["dHeat"] = ["u_p","agn_heat_p"]
    preqs["cs_p"] = ["u_p"]
    preqs["mJ_p"] = ["cs_p"]
    preqs["prat"] = ["p_p","nH_p","TK_p","mJ_p","m_p"]
    preqs["hz_rat"] = ["h_p","z_p"]
    preqs["vcirc"] = ["rad2d_p"]

    iVal = 0
    while ( iVal<len(requiredVals) ):
        reqVal = requiredVals[iVal]
        if reqVal in preqs:
            for preq in preqs[reqVal]:
                if ( preq not in requiredVals ):
                    requiredVals.append(preq)
        iVal+=1
        
    if not rcut is None and not "rad_p" in requiredVals:
        requiredVals.append("rad_p")

    values = dict()

    xyz = np.array(f["/PartType0/Coordinates"])
    rad_p = np.sqrt(xyz[:,0]**2+xyz[:,1]**2+xyz[:,2]**2)
    
    if ( "rad_p" in requiredVals ):
        values["rad_p"] = rad_p*1.e3

    if ( "rad2d_p" in requiredVals ):
        values["rad2d_p"] = np.sqrt(xyz[:,0]**2+xyz[:,1]**2)
        values["rad2d_p"]*=1.e3 # to pc

    vel_p = np.array(f["/PartType0/Velocities"])
    values["vrad"] = np.sum(vel_p*xyz,axis=1)/rad_p

    values["z_p"] = np.abs(xyz[:,2])
    values["z_p"]*=1.e3 # to pc
    
    if "vcirc" in requiredVals:
        values["vcirc"] = (vel_p[:,0]*xyz[:,1]-vel_p[:,1]*xyz[:,0])*1.e3/values["rad2d_p"]

    if "vel" in requiredVals:
        values["vel"] = np.sqrt(np.sum(vel_p**2,axis=1))

    if ( "opac" in requiredVals ):
        values["opac"] = np.array(f["/PartType0/AGNOpacity"]) # internal units: kpc**2/1e10 solar mass
        values["opac"]*= 1.e-4 # to pc**2/solar mass
        print(np.max(values["opac"]),np.min(values["opac"]))

    if ( "u_p" in requiredVals ):
        values["u_p"] = np.array(f["/PartType0/InternalEnergy"])
        values["u_p"]*=1.e10 # to erg/g

    if ( "tau" in requiredVals ):
        values["tau"] = np.array(f["/PartType0/AGNDepth"])

    if ( "rho_p" in requiredVals ):
        values["rho_p"] = np.array(f["/PartType0/Density"])
        values["rho_p"]*=6.77e-22 # to g/cm**3

    if ( "agn_heat_p" in requiredVals ):
        values["agn_heat_p"] = np.array(f["/PartType0/AGNHeat"])
        values["agn_heat_p"]*=1e10/3.08568e+16# to erg/s/g

    if ( "cool_p" in requiredVals ):
        values["cool_p"] = np.array(f["/PartType0/AGNHeat"])
        values["cool_p"]*=-1e10/3.08568e+16# to erg/s/g

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

    if ( "arad_p" in requiredVals ):
        accel_p = np.array(f["/PartType0/Acceleration"])
        accel_p*=3.0857e21/3.08568e+16**2 # to cm/s/s

        values["arad_p"] = (xyz[:,0]*accel_p[:,0]+xyz[:,1]*accel_p[:,1]+xyz[:,2]*accel_p[:,2])/rad_p

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
        print("Loading table and extracting values")
        chTab = tab_interp.CoolHeatTab( ("coolheat_tab_marta/shrunk_table_labels_291117tau.dat"),
                                            ("coolheat_tab_marta/shrunk_table_291117_m0.04_hsmooth_tau.dat"),
                                            ("coolheat_tab_marta/shrunk_table_labels_011217taunodust.dat"),
                                            ("coolheat_tab_marta/shrunk_table_011217_m0.04_hsmooth_taunodust.dat")
                                            )
        interpTabVec = np.vectorize(chTab.interpTab)
        tabStructs = interpTabVec(values["nH_p"].astype(np.float64),values["TK_p"].astype(np.float64),values["flux_p"].astype(np.float64),values["tau"].astype(np.float64))

    if ( "dustTemp" in requiredVals ):
        values["dustTemp"] = list(map(lambda y: y.dustT, tabStructs))
        values["dustTemp"] = np.array(values["dustTemp"])
        #print(list(values["dustTemp"]))

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
    yr = 31556926

    if ( "mJ_p" in requiredVals ):
        values["mJ_p"] = 2./3.*(np.pi)**2.5 * values["cs_p"]**3/G_Jeansterm/np.sqrt(values["rho_p"])
    #mJ_p/=1.989e33 # g to msun

#     if ( "m_p" in requiredVals ):
    values["m_p"] = np.array(f["/PartType0/Masses"])*1.e10
    # ratio between ideal gas pressure and actual pressure
    if ( "prat" in requiredVals ):
        values["prat"] = values["p_p"]/(values["nH_p"]*values["TK_p"]*boltzmann_cgs)*(values["mJ_p"]/values["m_p"])**(2./3.)

    # ratio between softening and z coordinate
    if ( "hz_rat" in requiredVals ):
        values["hz_rat"] = values["h_p"]/values["z_p"]

    if ( "tsf" in requiredVals ):
        values["tsf"] = 1./np.sqrt(values["rho_p"]*G)/yr

    return time,values

# assumes that values are already loaded
def plot_phaseplot(sp,values,iv,jv,rcut=None,noranges=False,cmap='plasma',bins=(150,150)):

    bigslice = (values["dt_p"]>0.)

    sp.set_xlabel(labels[iv])
    sp.set_ylabel(labels[jv])
    if ( dolog[iv] ):
        vx = np.log10(values[iv])
    else:
        vx = values[iv]
    if ( dolog[jv] ):
        vy = np.log10(values[jv])
    else:
        vy = values[jv]
    notnans = ((np.isfinite(vx)) & (np.isfinite(vy))) & bigslice
    
    print(np.nanmin(vx),np.nanmax(vx),np.nanmin(vy),np.nanmax(vy))
    
    if ( not noranges and not (ranges[iv] is None) ):
        notnans = notnans & (vx>=ranges[iv][0]) & (vx<=ranges[iv][1])
    if ( not noranges and not (ranges[jv] is None) ):
        notnans = notnans & (vy>=ranges[jv][0]) & (vy<=ranges[jv][1])
    if rcut is not None:
        notnans = notnans & (values["rad_p"]<rcut)
    vx = vx[notnans]
    vy = vy[notnans]
    # CHECK FOR RANGES
    if ( not (ranges[iv] is None) and not (ranges[jv] is None) and not noranges ):
        H,xedges,yedges = np.histogram2d(vx,vy,bins=bins,range=[ranges[iv],ranges[jv]])
    else:
        H,xedges,yedges = np.histogram2d(vx,vy,bins=bins)
    H*=values["m_p"][0]
    with np.errstate(divide='ignore'):
        H = np.log10(H).T

    mappablePlot = sp.pcolormesh(xedges,yedges,H,cmap=cmap,vmin=-1.4,vmax=2.1) #,norm=colors.LogNorm()
    if iv == "rad_p" and jv == "vel":
        # plot escape velocities
#                 print("PLOTTING ESCAPE VELOCITY")
        y = v_esc(10.**xedges,1.e6,1.e9,250.)
        sp.plot(xedges,y)
    sp.set_xlim(xedges[0],xedges[-1])
    sp.set_ylim(yedges[0],yedges[-1])
    sp.tick_params(which='both',direction="out")
    
    return mappablePlot


def savephaseplots(run_id,output_dir,snap_str,includedVals,rcut=None):
    print("plotting:"+"".join([run_id,output_dir,snap_str,"".join(includedVals)]))
    
    time,values = loadvalues(run_id,output_dir,snap_str,includedVals,rcut)
    
    
    nv = len(values)

#     noranges = False

    outfiles = []

    i = 0
    #for iv in ["dt_p"]:
    for i,iv in enumerate(includedVals):
        for jv in includedVals[i+1:]:
    #for iv in range(7):
    #    for jv in [7]:
            fname = "pics/"+run_id+output_dir+iv+jv+snap_str+".png"
            outfiles.append(fname)
            fig,sp = P.subplots(1,1)
            
            mappablePlot = plot_phaseplot(sp,values,iv,jv,rcut)
            
#             ax = P.gca()
            #ax.xaxis.set_minor_locator(ticker.MultipleLocator(1.))
            #ax.yaxis.set_minor_locator(ticker.MultipleLocator(1.))
#             P.colorbar(label=r"$\log$ count")
            fig.colorbar(mappablePlot,label=r"$\log M$ (M$_\odot$)")
            #P.colorbar(label=r"count")
            #P.colorbar(label=r"count (capped)")
            fig.suptitle(r"$t="+("%.4f" % time)+"$ Myr")
            fig.savefig(fname,dpi=150)
            
            P.close()
    return outfiles


if __name__ == '__main__':
    run_id = sys.argv[1]
    output_dir = sys.argv[2]
    snap_str = sys.argv[3]

    print("Running")


    #includedVals = ["rad_p","radrad_p"]
#     includedVals = ["nH_p","TK_p","vel"]
    #includedVals = ["TK_p","agn_heat_p"]
#     includedVals = ["rad_p","vel"]
    #includedVals = ["TK_p","vrad"]
#     includedVals = ["arad_p","radrad_p"]
#     includedVals = ["TK_p","dustTemp","vrad"]
#     includedVals = ["TK_p","dustTemp","nH_p"]
#     includedVals = ["rad_p","opac"]
#     includedVals = ["TK_p","dustTemp","nH_p"]
    includedVals = ["dt_p","TK_p","nH_p","vrad","rad_p","agn_heat_p","cool_p","h_p"]
    
    x = savephaseplots(run_id,output_dir,snap_str,includedVals,rcut=80.)
    print(x)

