print("Importing")

# Import libraries to do our magic
import numpy as np
import h5py
import sys
import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as P

import os
this_dir, this_filename = os.path.split(__file__)

from . import gizmo_tools

import time as timer
import pandas as pd

from scipy import stats

# luminosity = 6.2849e42 # erg/s

los_angle = np.radians(10.)
cos_los = np.cos(los_angle)
sin_los = np.sin(los_angle)

G_kms_pc_msun = 0.0043022682

labels = dict()
dolog = dict()
ranges = dict()

lines = ["co1","co2","hcn1","hcn2","h2_1","h2_2","h2_3"]
lineLabels = ["CO(3-2)","CO(6-5)","HCN(4-3)","HCN(8-7)","H2 (1-0) S(1)","H2 (0-0) S(0)","H2 0-0 S(3)"]

for iline,lineVal in enumerate(lines):
    labels[lineVal] = lineLabels[iline]+" (erg/g/s)"
    dolog[lineVal] = True
#     ranges[lineVal] = [-10.,0.] # massy
#     ranges[lineVal] = [-40.,-20.] # volumetric
    ranges[lineVal] = [-20.,2.] # massy, big range

labels['nH_p'] = r"$\log n_H$ (cm$^{-3}$)"
dolog['nH_p'] = True
# ranges['nH_p'] = [-4,8]
# ranges['nH_p'] = [-4,12]
ranges['nH_p'] = [2,12]
# ranges['nH_p'] = [6,9]

labels['rho_p'] = r"$\log \rho$ (g cm$^{-3}$)"
dolog['rho_p'] = True
ranges['rho_p'] = [-35.,-19.]

labels['age'] = r"age (Myr)"
dolog['age'] = False
ranges['age'] = [0.,1.5]


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
ranges['TK_p'] = [1.,8.]
# ranges['TK_p'] = [1.,5.]
# ranges['TK_p'] = [.9,8.1]
# ranges['TK_p'] = [.9,8.1]

labels['rad2d_p'] = r"$R$ (pc)"
dolog['rad2d_p'] = False
# ranges['rad2d_p'] = [0,.05]
ranges['rad2d_p'] = [0,5.]
# ranges['rad2d_p'] = None

labels['z_p'] = r"$|z|$ (pc)"
dolog['z_p'] = False
ranges['z_p'] = [0,20.]

labels['dustTemp'] = r"$\log T_d$ (K)"
dolog['dustTemp'] = True
# ranges['dustTemp'] = [1.,4.]#[0.,0.5] #[0,150]
ranges['dustTemp'] = [1.,2.5]#[0.,0.5] #[0,150]
# ranges['dustTemp'] = [1.3,2.3]#[0.,0.5] #[0,150]

labels['radrad_p'] = r"$a_{rad,r}$ (cm/s/s)"
dolog['radrad_p'] = True
ranges['radrad_p'] = [-10.,1.]

labels['arad_p'] = r"$a_{all,r}$ (cm/s/s)"
dolog['arad_p'] = 1.e7
ranges['arad_p'] = [-13.,13.]

labels['dHeat'] = r"$H$"
dolog['dHeat'] = 10.
ranges['dHeat'] = None

labels['dt_heat'] = r"dt$_H$"
dolog['dt_heat'] = True
ranges['dt_heat'] = None


labels['vrad'] = r"$v_{rad}$ (km/s)"
dolog['vrad'] = False
# #ranges['vrad'] = [-300,300]
# ranges['vrad'] = [-300,1.e3]
# ranges['vrad'] = [-100,300.]
ranges['vrad'] = [-200,200.]
# # ranges['vrad'] = [-50,300.]
# # ranges['vrad'] = [-50,600.]
# ranges['vrad'] = [-4.,4.]

labels['vradpoly'] = r"$v_{rad}$ (km/s)"
dolog['vradpoly'] = False
# #ranges['vrad'] = [-300,300]
# #ranges['vrad'] = [-300,1.e3]
# # ranges['vrad'] = [-50,300.]
# # ranges['vrad'] = [-50,600.]
ranges['vradpoly'] = [-4.,22.]

labels['vlos'] = r"$v_{los}$ (km/s)"
dolog['vlos'] = False
ranges['vlos'] = [-200.,200.]

labels['xlos'] = r"$x_{los}$ (pc)"
dolog['xlos'] = False
ranges['xlos'] = [-20.,20.]

labels['ylos'] = r"$y_{los}$ (pc)"
dolog['ylos'] = False
ranges['ylos'] = [-20.,20.]

# polyfac = 3. # take cube root - need to do this explicitly

# labels['vrad'] = r"$v_{rad}$ (km/s)"
# dolog['vrad'] = False
# ranges['vrad'] = [-50.,50.]


labels['vcirc'] = r"$v_{circ}$ km/s"
dolog['vcirc'] = False
ranges['vcirc'] = [-20.,150.]

# labels['vel'] = r"$v$ km/s"
labels['vel'] = r"$\log_{10} v$ km/s"
# dolog['vel'] = True
dolog['vel'] = False
# ranges['vel'] = [0.,1.e3]
ranges['vel'] = [-1.e5,1.e5]
# ranges['vel'] = None

# labels['vel'] = r"$\log v$ km/s"
# dolog['vel'] = True
# ranges['vel'] = [1.,4.]

# labels['rad_p'] = r"$\log R$ (pc)"
# dolog['rad_p'] = True
# ranges['rad_p'] = [-1.,1.9]

labels['rad_p'] = r"$R$ (pc)"
dolog['rad_p'] = True
# ranges['rad_p'] = [0.,.3]
# ranges['rad_p'] = [0.,10.]
# ranges['rad_p'] = [0.,20.]
ranges['rad_p'] = [-3.,3.]

labels['tau'] = r"$\tau$"
dolog['tau'] = False
ranges['tau'] = [0.,7.]

labels['mJ_p'] = r"$M_\mathrm{J}$ (M$_\odot$)"
dolog['mJ_p'] = True
#ranges['mJ_p'] = None
ranges['mJ_p'] = [-2.,8.]

labels['p_p'] = r"$P$ (dyne/cm$^{-2}$)"
dolog['p_p'] = True
ranges['p_p'] = [-15.,-2.]

labels['prad'] = r"$P_{rad}$ (dyne/cm$^{-2}$)"
dolog['prad'] = True
ranges['prad'] = [-15.,-2.]

labels['pratio'] = r"$P_{rad}/P_{thermal}$"
dolog['pratio'] = True
# ranges['pratio'] = [-5.,5.]
ranges['pratio'] = [-10.,3.]

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
ranges['agn_heat_p'] = [-1.e3,1.e3]

labels['cool_p'] = r"$H$ (erg/s/g)"
dolog['cool_p'] = True
ranges['cool_p'] = None

labels['opac'] = r"$\kappa$ (pc$^2$/M$_\odot$)"
dolog['opac'] = True
ranges['opac'] = [-5.,-1.]

labels['phi'] = r"$\phi$ ($^\deg$)"
dolog['phi'] = False
ranges['phi'] = [0.,90.]

def symlog(a,threshold=100.):
    assert(threshold>=1.)
#     log_threshold = np.log10(threshold)
    b = a.copy()
    lin_slice = (a>=-threshold) & (a<=threshold)
    pos_log_slice = (a>threshold)
    neg_log_slice = (a<-threshold)
    b[lin_slice] = a[lin_slice]*np.log10(threshold)/threshold
    b[pos_log_slice] = np.log10(a[pos_log_slice])
    b[neg_log_slice] = -np.log10(-a[neg_log_slice])
    return b

def v_esc(r,m_bh,m_hern,a_hern):
    return np.sqrt(2*G_kms_pc_msun*(m_bh/r + m_hern/(a_hern+r)))

def loadvalues(run_id,output_dir,snap_str,includedVals,rcut=None):
    gizmoDir = gizmo_tools.getGizmoDir(run_id)
    fullDir = gizmoDir+"/"+run_id+"/"+output_dir
    fullFile = fullDir+"/snapshot_"+snap_str+".hdf5"

    f = h5py.File(fullFile,"r")

    header = f["/Header"]
    time = header.attrs.get("Time")
    time*= 0.9778e9 # to yr
    time/=1.e6 # to Myr

    requiredVals = list(includedVals)

    if ( "rho_p" not in requiredVals ):
        requiredVals.append("rho_p")
    if ( "dt_p" not in requiredVals ):
        requiredVals.append("dt_p")

#     requiredVals.append("hcn2") # for great justice


    preqs = dict()
    preqs["TK_p"]=["u_p"]
    if ( "/PartType0/Pressure" not in f ):
        preqs["p_p"]=["nH_p","TK_p"]
    for x in ["dustTemp","dHeat"]+lines:
        preqs[x] = ["table"]
    preqs["table"] = ["nH_p","TK_p","flux_p","tau"]
    preqs["dHeat"]+= ["u_p","agn_heat_p"]
    preqs["cs_p"] = ["u_p"]
    preqs["mJ_p"] = ["cs_p"]
    preqs["prat"] = ["p_p","nH_p","TK_p","mJ_p","m_p"]
    preqs["hz_rat"] = ["h_p","z_p"]
    preqs["vcirc"] = ["rad2d_p"]
    preqs["radrad_p"] = ["radaccel_p"]
    preqs["pratio"] = ["prad","p_p"]
    preqs["prad"] = ["m_p","radaccel_p","h_p"]
    preqs["phi"] = ["rad2d_p"]
    preqs["xlos"] = ["ylos","h_p"]
    preqs["ylos"] = ["xlos","h_p"]
#     preqs["age"] = ["id_p"]

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

#     if ( "m_p" in requiredVals ):
    values["m_p"] = np.array(f["/PartType0/Masses"])*1.e10
#     
#     if "id_p" in requiredVals:
# 

    if "age" in requiredVals:
        run_isnap = int(snap_str)
        age_file = "data/age_{}_{}_{}.dat".format(run_id,output_dir,run_isnap)
        values["age"]= time-np.loadtxt(age_file)
    
    if "phi" in requiredVals:
        values["phi"] = np.arctan(np.abs(xyz[:,2]*1.e3/values["rad2d_p"]))*180./np.pi
    
    if "vcirc" in requiredVals:
        values["vcirc"] = (vel_p[:,0]*xyz[:,1]-vel_p[:,1]*xyz[:,0])*1.e3/values["rad2d_p"]

    if "vlos" in requiredVals:
        values["vlos"] = vel_p[:,0]*sin_los + vel_p[:,1]*cos_los
    
    if "xlos" in requiredVals:
        values["xlos"] = xyz[:,0]*1.e3 # to pc

    if "ylos" in requiredVals:
        values["ylos"] = (xyz[:,1]*sin_los + xyz[:,2]*cos_los)*1.e3 # to pc
    
    if "vradpoly" in requiredVals:
        values["vradpoly"]=np.cbrt(values["vrad"])

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

    if ( "radaccel_p" in requiredVals ):
        radaccel_p = np.array(f["/PartType0/RadiativeAcceleration"])
        radaccel_p*=3.0857e21/3.08568e+16**2 # to cm/s/s

    if ( "radrad_p" in requiredVals ):
        values["radrad_p"] = (xyz[:,0]*radaccel_p[:,0]+xyz[:,1]*radaccel_p[:,1]+xyz[:,2]*radaccel_p[:,2])/rad_p

    if ( "prad" in requiredVals ):
        values["prad"] = np.sqrt(np.sum(radaccel_p**2,axis=1))
        values["prad"]*=values["m_p"]/(np.pi*values["h_p"]**2) * 0.00208908219 # units

    if ( "arad_p" in requiredVals ):
        accel_p = np.array(f["/PartType0/Acceleration"])
        accel_p*=3.0857e21/3.08568e+16**2 # to cm/s/s
        
        accel_p*=9.95839577e14 # to cm/yr/yr

        values["arad_p"] = -(xyz[:,0]*accel_p[:,0]+xyz[:,1]*accel_p[:,1]+xyz[:,2]*accel_p[:,2])/rad_p

    molecular_mass = 4./(1.+3.*.76)
    proton_mass_cgs = 1.6726e-24
    gamma_minus_one = 5./3.-1.
    boltzmann_cgs = 1.38066e-16

    print(gamma_minus_one/boltzmann_cgs*(molecular_mass*proton_mass_cgs)*1.e10)

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

    if ( "pratio" in requiredVals ):
        values["pratio"] = values["prad"]/values["p_p"]



    if ( "table" in requiredVals ):
#         tableDate="281118"
#         tableRes="0.0001"

        tableDate="130720" 
        tableRes="0.001"

#         print("Loading table and extracting values")
#         start = timer.time()
#         chTab = tab_interp.CoolHeatTab( ("coolheat_tab/shrunk_table_labels_"+tableDate+"tau.dat"),
#                                         ("coolheat_tab/shrunk_table_"+tableDate+"_m"+tableRes+"_hsmooth_tau.dat"),
#                                         ("coolheat_tab/shrunk_table_labels_"+tableDate+"taunodust.dat"),
#                                         ("coolheat_tab/shrunk_table_"+tableDate+"_m"+tableRes+"_hsmooth_taunodust.dat"),
#                                         ("coolheat_tab/shrunk_table_labels_"+tableDate+"taudense.dat"),
#                                         ("coolheat_tab/shrunk_table_"+tableDate+"_m"+tableRes+"_hsmooth_taudense.dat")
#                                         )
#         end = timer.time()
#         print("Table loaded",end-start)
# 
#         print("Interpolating from table")
#         start = timer.time()
#         interpTabVec = np.vectorize(chTab.interpTab)
#         tabStructs = interpTabVec(values["nH_p"].astype(np.float64),values["TK_p"].astype(np.float64),values["flux_p"].astype(np.float64),values["tau"].astype(np.float64))
#         end = timer.time()
#         print("Interpolation complete",end-start)

        print("Loading table and extracting values")
        start = timer.time()
        cloudy_table = gizmo_tools.cloudy_table(tableDate,tableRes,"coolheat_tab/")
        end = timer.time()
        print("Table loaded",end-start)
        
        print("Interpolating from table")
        start = timer.time()
        table_particles = pd.DataFrame()
        table_particles["nH"] = values["nH_p"]
        table_particles["temp"] = values["TK_p"]
        table_particles["AGNIntensity"] = values["flux_p"]
        table_particles["AGNDepth"] = values["tau"]

        cloudy_table.interp(table_particles)
        end = timer.time()
        print("Interpolation complete",end-start)

    if ( "dustTemp" in requiredVals ):
#         values["dustTemp"] = np.array(list(map(lambda y: y.dustT, tabStructs)))
        values["dustTemp"] = table_particles["dustT"]
#         data.dustTemp = map(lambda y: y.dustT, tabStructs)
#         data.dustTemp = np.array(list(data.dustTemp))
#         data.dustTemp = 10.**data.dustTemp
#         print(np.max(data.dustTemp),np.min(data.dustTemp))
#         data.dustTemp = data.dustTemp**4 # for test


#         values["dustTemp"] = 10.**np.array(values["dustTemp"])
        #print(list(values["dustTemp"]))
    
    
#     if ( "co1" in requiredVals ):
#         values["co1"] = np.array(list(map(lambda y: y.line_co1, tabStructs)))#/values["rho_p"]
#     if ( "co2" in requiredVals ):
#         values["co2"] = np.array(list(map(lambda y: y.line_co2, tabStructs)))#/values["rho_p"]
#     if ( "hcn1" in requiredVals ):
# #         values["hcn1"] = np.array(list(map(lambda y: y.line_hcn1, tabStructs)))#/values["rho_p"]
#         values["hcn1"] = np.array(list(map(lambda y: y.__getattr__("line_hcn1"), tabStructs)))#/values["rho_p"]
#     if ( "hcn2" in requiredVals ):
#         values["hcn2"] = np.array(list(map(lambda y: y.line_hcn2, tabStructs)))#/values["rho_p"]
    
    for lineVal in lines:
        if lineVal in requiredVals:
            values[lineVal] = table_particles["line_"+lineVal]

#             values[lineVal] = np.array(list(map(lambda y: y.__getattr__("line_"+lineVal), tabStructs)))


    if ( "dHeat" in requiredVals ):
        values["dHeat"] = table_particles["dHeat"]
#         values["dHeat"] = map(lambda y: y.dHeat, tabStructs)
#         values["dHeat"] = np.array(values["dHeat"])

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

    # ratio between ideal gas pressure and actual pressure
    if ( "prat" in requiredVals ):
        values["prat"]=values["p_p"]/(values["nH_p"]*values["TK_p"]*boltzmann_cgs)*(values["mJ_p"]/values["m_p"])**(2./3.)

    # ratio between softening and z coordinate
    if ( "hz_rat" in requiredVals ):
        values["hz_rat"] = values["h_p"]/values["z_p"]

    if ( "tsf" in requiredVals ):
        values["tsf"] = 1./np.sqrt(values["rho_p"]*G)/yr


    # HACKY CUT
#     print("HACKY CUT")
#     # for plotting ylos
#     slice = (values["h_p"]>np.abs(values["xlos"]))
#     # for plotting xlos
# #     slice = (values["h_p"]>np.abs(values["ylos"]))
#     for key in values:
#         values[key] = values[key][slice]
#     # for plotting xlos

    return time,values

# assumes that values are already loaded
# def plot_phaseplot(sp,values,iv,jv,rcut=None,noranges=False,cmap='plasma',bins=(150,150),m_bh=1.e6):
def plot_phaseplot(sp,values,iv,jv,rcut=None,noranges=False,vrange=None,cmap='plasma',bins=(300,300),m_bh=1.e6,weight=None,slice=None,drawmed=False):

    bigslice = (values["dt_p"]>0.)
    
#     bigslice = (bigslice) & (values["vrad"]<0.) & (values["rad_p"]>2.) & (values["z_p"]>0.5); print("MAKING BIG SLICE")
#     bigslice = (bigslice) & (values["hcn2"]>10.**-4)
    val_i = values[iv]
    val_j = values[jv]

    sp.set_xlabel(labels[iv])
    sp.set_ylabel(labels[jv])
    if type(dolog[iv]) is float:
        vx = symlog(val_i,dolog[iv])
    elif dolog[iv]:
        vx = np.log10(val_i)
    else:
        vx = val_i
    if type(dolog[jv]) is float:
        vy = symlog(val_j,dolog[jv])
    elif dolog[jv]:
        vy = np.log10(val_j)
    else:
        vy = val_j
    notnans = ((np.isfinite(vx)) & (np.isfinite(vy))) & bigslice
    
    print(np.nanmin(vx),np.nanmax(vx),np.nanmin(vy),np.nanmax(vy))
    
    if ( not noranges and not (ranges[iv] is None) ):
        notnans = notnans & (vx>=ranges[iv][0]) & (vx<=ranges[iv][1])
    if ( not noranges and not (ranges[jv] is None) ):
        notnans = notnans & (vy>=ranges[jv][0]) & (vy<=ranges[jv][1])
    if rcut is not None:
        notnans = notnans & (values["rad_p"]<rcut)
    
    if slice is not None:
        notnans = notnans & slice

    print("binning")
    vx = vx[notnans]
    vy = vy[notnans]
    # CHECK FOR RANGES
    if ( not (ranges[iv] is None) and not (ranges[jv] is None) and not noranges ):
        rangepair = [ranges[iv],ranges[jv]]
    else:
        rangepair = None
    H,xedges,yedges = np.histogram2d(vx,vy,bins=bins,range=rangepair) # N per bin
    if weight is not None:
        vweight = (values[weight])[notnans]
        print(np.min(vweight),np.max(vweight))
        Hweight,xedges,yedges = np.histogram2d(vx,vy,bins=bins,range=rangepair,weights=vweight) # N per bin
        print(np.max(Hweight),np.max(H))
        print(np.min(Hweight),np.min(H))
        H=Hweight/H
    else:
        H*=values["m_p"][0] # mass per bin
#     if ( not (ranges[iv] is None) and not (ranges[jv] is None) and not noranges ):
#         H,xedges,yedges = np.histogram2d(vx,vy,bins=bins,range=[ranges[iv],ranges[jv]])
#     else:
#         H,xedges,yedges = np.histogram2d(vx,vy,bins=bins)
    with np.errstate(divide='ignore'):
        H = np.log10(H).T

    if vrange is not None:
        vmin = vrange[0]
        vmax = vrange[1]
    elif weight is not None:
        vmin = np.nanmin(H)
        vmax = np.nanmax(H)
        print(vmin,vmax)
    else:
        vmin = np.log10(values["m_p"][0]/2.)
        vmax = np.log10(np.sum(values["m_p"]))
#     if iv=="vradpoly":
#         xedges=xedges**polyfac
#     if jv=="vradpoly":
#         yedges=yedges**polyfac
    print("Plotting")
    mappablePlot = sp.pcolormesh(xedges,yedges,H,cmap=cmap,vmin=vmin,vmax=vmax) #,norm=colors.LogNorm()
    if iv == "rad_p" and jv == "vel":
        # plot escape velocities
#                 print("PLOTTING ESCAPE VELOCITY")
#         y = v_esc(10.**xedges,m_bh,1.e9,250.)
        y = v_esc(xedges,m_bh,1.e9,250.)
        sp.plot(xedges,y)
    output = mappablePlot
    if drawmed:
        y,bin_edges,bindices = stats.binned_statistic(vx,vy,statistic='median',bins=xedges)
        print(y)
        xmids = (xedges[:-1]+xedges[1:])/2
        sp.plot(xmids,y)
        output = mappablePlot,xmids,y
    
    
    sp.set_xlim(xedges[0],xedges[-1])
    sp.set_ylim(yedges[0],yedges[-1])
    sp.tick_params(which='both',direction="out")
    if iv=="vradpoly" or jv=="vradpoly":
        newtick_vals = np.array([-100.,-10.,0.,10.,100.,1000.,5000.])
        newtick_locs = np.cbrt(newtick_vals)
        newtick_vals = ["%d"%x for x in newtick_vals]
        if iv=="vradpoly":
            sp.set_xticks(newtick_locs)
            sp.set_xticklabels(newtick_vals,rotation=90)
        else:
            sp.set_yticks(newtick_locs)
            sp.set_yticklabels(newtick_vals)
    
    return output


def savephaseplots(run_id,output_dir,snap_str,includedVals,rcut=None,m_bh=1.e6,gridMode=True,uniqueMode=False,weights=None):
    print("plotting:"+"".join([run_id,output_dir,snap_str,"".join(includedVals)]))
    if weights is not None:
        if type(weights) is list:
            loadingVals = includedVals+weights
        else:
            loadingVals = includedVals+[weights]
    else:
        weights = [None]
        loadingVals = includedVals
    
    time,values = loadvalues(run_id,output_dir,snap_str,loadingVals,rcut)
    

    
    nv = len(values)

#     noranges = False

    outfiles = []

    i = 0
    #for iv in ["dt_p"]:
    varOneRange = includedVals if gridMode else [includedVals[0]]
    
#     if uniqueMode:
#         varOneRange = includedVals[0::2]

    pic_dir = gizmo_tools.getPicDir()


    for i,iv in enumerate(varOneRange):
        if uniqueMode and i%2!=0:
            continue
        varTwoRange = [includedVals[i+1]] if uniqueMode else includedVals[i+1:] 
        for jv in varTwoRange:
            for weight in weights:
                if weight is not None:
                    weight_str = weight
                else:
#                     weights = [None]
                    weight_str = ""
        #for iv in range(7):
        #    for jv in [7]:
                fname = os.path.join(pic_dir, run_id + output_dir + iv + jv + weight_str + snap_str + ".png")
                outfiles.append(fname)
                fig,sp = P.subplots(1,1)
            
                mappablePlot = plot_phaseplot(sp,values,iv,jv,rcut,m_bh=m_bh,weight=weight)
            
    #             ax = P.gca()
                #ax.xaxis.set_minor_locator(ticker.MultipleLocator(1.))
                #ax.yaxis.set_minor_locator(ticker.MultipleLocator(1.))
    #             P.colorbar(label=r"$\log$ count")
                if weight is None:
                    fig.colorbar(mappablePlot,label=r"$\log M$ (M$_\odot$)")
                else:
                    fig.colorbar(mappablePlot,label=labels[weight])
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
#     includedVals = ["dt_p","TK_p","nH_p","vrad","rad_p","agn_heat_p","cool_p","h_p"]
#     includedVals = ["nH_p","agn_heat_p"]
#     includedVals = ["nH_p","TK_p"]
    includedVals = ["nH_p","dustTemp"]
#     includedVals = ["rho_p","dustTemp"]
#     includedVals = ["tau_p","dustTemp"]
#     x = savephaseplots(run_id,output_dir,snap_str,includedVals,rcut=80.,weights=lines)
    x = savephaseplots(run_id,output_dir,snap_str,includedVals,rcut=80.)
    
#     x = savephaseplots(run_id,output_dir,snap_str,includedVals,rcut=80.)
    print("saving to ",x)

