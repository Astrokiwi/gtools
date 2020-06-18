import socket
import getpass
import os
import numpy as np
import re
import pandas as pd
import h5py
import pickle

import pynbody
import matplotlib.cm
import string

from scipy import interpolate

import matplotlib.patches as patches

# for interpolating GIZMO tables
import ctypes
from sys import path
path.append("src/")
path.append("../src/")
import tab_interp

import matplotlib.pyplot as plt

molecular_mass = 4./(1.+3.*.76)
proton_mass_cgs = 1.6726e-24
gamma_minus_one = 5./3.-1.
boltzmann_cgs = 1.38066e-16
msun_g = 1.989e33


unit_conversions = None
unit_names = None

coord_suffixes = "xyz"

#https://gist.github.com/thriveth/8560036 - colourblind colours from here
colors = ['#377eb8', '#ff7f00', '#4daf4a',
          '#f781bf', '#a65628', '#984ea3',
          '#999999', '#e41a1c', '#dede00']
ncolors = len(colors)
markers = ["x","+",'|']
cmap = matplotlib.cm.get_cmap('coolwarm')

dust_opacity_function = None

binary_headers_unit_name_conversions = {
        "Binary_pos_1":("BH_pos_1",1000.,"pc"), # pc 
        "Binary_pos_2":("BH_pos_2",1000.,"pc"), # pc
        "Binary_vel_1":("BH_vel_1",1000./0.9778e9,"pc yr**-1"), # pc/yr
        "Binary_vel_2":("BH_vel_2",1000./0.9778e9,"pc yr**-1"), # pc/yr
        "Binary_force_1":("BH_force_1",1e10 * 1000.0/(0.9778e9)**2,"Msol pc yr**-2"),   # msun * pc / yr**2
        "Binary_force_2":("BH_force_2",1e10 * 1000.0/(0.9778e9)**2,"Msol pc yr**-2"),   # msun * pc / yr**2
        "Binary_mass_1":("BH_mass_1",1.e10,"Msol"), # msun
        "Binary_mass_2":("BH_mass_2",1.e10,"Msol")  # msun

    }

def load_interpolate_opacity(opac_mu,
                        opac_file="prams/simple_dust_opacity.txt"):
    global dust_opacity_function
    if dust_opacity_function is None:
        dust_opacity_table = np.loadtxt(opac_file)
        dust_opacity_function = interpolate.interp1d(dust_opacity_table[:,0],dust_opacity_table[:,1])
    opacity=dust_opacity_function(opac_mu)
#     print("Setting opacity to ",opacity," cm^2/g")
    opacity*=0.000208908219 # convert to pc**2/solar mass for consistency
    return opacity

line_bases = ["co1","co2","hcn1","hcn2","h2_1","h2_2","h2_3","12mic","8mic","850mic"]
line_codes = ["line_"+line for line in line_bases]
lineWavelengths_list = ['866.727', '433.438', '845.428', '422.796', '2.121', '28.18', '9.66', '12', '8', '850']
lineNames_list = ['CO(3-2)',
 'CO(6-5)',
 'HCN(4-3)',
 'HCN(8-7)',
 'H$_2$ (1-0) S(1)',
 'H$_2$ (0-0) S(0)',
 'H$_2$ (0-0) S(3)',
 'Continuum',
 'Continuum',
 'Continuum']
line_wavelengths = [866.727,433.438,845.428,422.796,2.121,28.18,9.66,12,8,850]

lineNamesFull = {code: name+f" ${wavelength}$ $\\mu$m"
                for code,wavelength,name in zip(line_bases,lineWavelengths_list,lineNames_list)}

def derive_opacities(opac_path=""):
    global line_opacities
    line_opacities = {line:load_interpolate_opacity(mu,opac_file=opac_path+"prams/simple_dust_opacity.txt") for line,mu in zip(line_bases,line_wavelengths)}


def sort_nicely( l ):
    """ Sort the given list in the way that humans expect.
    """
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    l.sort( key=alphanum_key )

# def getGizmoDir():
#     sname = socket.gethostname()
# 
#     if ( sname=="trillian" ):
#         return "/export/1/djw/gizmos"
#     elif ( sname=="srv01921" ):
#         return "/srv/djw1g16/gizmos"
#     
#     raise Exception("Unknown server; add server and directory to gizmodatadir.py")

def box_connected_two_axes(small_ax,big_ax,box_corner,box_dims,color='blue'):
    # draw boxes
    for ax in [small_ax,big_ax]:
        ax.add_patch(patches.Rectangle( box_corner,box_dims[0],box_dims[1],fill=False,color=color))
    # link them up
    for i0 in range(2):
        for j0 in range(2):
            corner_coordinates = (box_corner[0]+box_dims[0]*i0,box_corner[1]+box_dims[1]*j0)
            a=big_ax.add_artist(patches.ConnectionPatch(
                        xyA=corner_coordinates
                        ,xyB=corner_coordinates
                        ,coordsA="data"
                        ,coordsB="data"
                        ,axesA=big_ax
                        ,axesB=small_ax
                        ,color=color
                        ))
            a.set_in_layout(False) #otherwise matplotlib tries to stretch the plot to "fit" the lines in, and the imshow maps become tiny
    

def getGizmoDir(irun):
    sname = socket.gethostname()
    uname = getpass.getuser()
    
    if uname=='djw1g16' or uname=='djw':
        if sname=="trillian":
            paths = ["/export/1/djw/gizmos","/export/2/djw/gizmos"]
            for path in paths:
                if os.path.isdir("{}/{}".format(path,irun)):
                    return path
    #         return "/export/1/djw/gizmos"
            raise Exception("Directory not found on {}".format(sname))
        elif sname=="srv01921":
            return "/srv/djw1g16/gizmos"
        elif "cluster.local" in sname:
            return "/mainfs/VEILS/djw1g16/gizmos"

    if uname=='lb1g19':
        if sname=="trillian":
            return "/export/2/lb1g19"

    if uname=='supas356':
        if sname in ("nesh-fe1","nesh-fe2","nesh-fe3","nesh-fe4"):
            return "/sfs/fs2/work-sh1/supas356"

        if sname in ("neshcl226","neshcl227"):
            return "/sfs/fs2/work-sh1/supas356"
          
    raise Exception(f"Unknown server/username; add server ({sname}), username ({uname}) and directory (?) to getGizmoDir gizmo_tools.py")

def getMovieDir():
    sname = socket.gethostname()
    uname = getpass.getuser()

    if uname=='djw1g16' or uname=='djw':
        if ( sname=="trillian" ):
            return "/export/1/djw/movies"
        elif ( sname=="srv01921" ):
            return "/srv/djw1g16/movies"
        elif "cluster.local" in sname:
            return "/mainfs/VEILS/djw1g16/movies"

    if uname=='lb1g19':
        if sname=="trillian":
            return "/export/2/lb1g19"

    if uname=="supas356":
        if sname in ("nesh-fe1","nesh-fe2","nesh-fe3","nesh-fe4"):
            return "/sfs/fs2/work-sh1/supas356/movies"
        if sname in ("neshcl226","neshcl227"):
          return "/sfs/fs2/work-sh1/supas356/movies"

    raise Exception(f"Unknown server/username; add server ({sname}), username ({uname}) and directory (?) to getMovieDir gizmo_tools.py")

def lastConsecutiveSnapshot(run_id,output_dir,dumpsOrdered=True,gizmoDir=None):
    """When models are rerun, the snapshot file with the largest number (e.g. snapshot_100.dat)
    may be from a previous model. So we want the last snapshot that was made *after* the previous
    snapshot. However, this only works if the snapshots haven't been copied, because this kills
    the "last modified" data. So we can set dumpsOrdered=False and force us to just take the last
    snapshot numerically rather than by time"""

    if gizmoDir is None:
        gizmoDir = getGizmoDir(run_id)
    movieDir = getMovieDir()
    fullDir = gizmoDir+"/"+run_id+"/"+output_dir

    fnames = os.listdir(fullDir)
    sort_nicely(fnames)
    fnames = np.array(fnames)
    #fnames.sort()
    snapshotfilebools = np.array([x.startswith("snapshot") for x in fnames])
    snapshotfiles = fnames[snapshotfilebools]
    
    if not dumpsOrdered:
        fname = snapshotfiles[-1]
        snapf = int(fname[9:len(fname)-5])
        return snapf
#     print(snapshotfiles)
    
    snapf = 0

    ctime = os.path.getmtime(fullDir+"/snapshot_000.hdf5")
    for fname in snapshotfiles[1:]:
        new_snapf = int(fname[9:len(fname)-5])
        new_ctime = os.path.getmtime(fullDir+"/"+fname)
#         print(ctime,new_ctime)
        if new_ctime>ctime and snapf==new_snapf-1 :
            ctime = new_ctime
            snapf = new_snapf
    
    return snapf

def load_unit_conversions():
    global unit_conversions
    if unit_conversions is not None:
        return
    unit_conversions = dict()
    unit_names = dict()
    with open("prams/unit_conversions.dat") as f:
        for line in f:
            stripped_line = line.strip()
            if len(stripped_line)<=2:
                continue
            if stripped_line[0]=="#":
                continue
            split_line = stripped_line.split()
            unit_conversions[split_line[0]] = float(split_line[1])
            unit_names[split_line[0]] = split_line[2]
#     print(unit_conversions)

def check_requirements(gizmo_dataframe,reqs):
    if type(reqs)==str:
        return reqs in gizmo_dataframe
    # otherwise, assume iterable
    for req in reqs:
        if not (req in gizmo_dataframe):
            return False
    return True

def load_calc_reqs(gizmo_dataframe,reqs):
    if type(reqs) == str:
        return load_calc_req(gizmo_dataframe,reqs)
    # otherwise assume iterable
    status = True
    for req in reqs:
        status&=load_calc_req(gizmo_dataframe,req)
    return status
        
def load_calc_req(gizmo_dataframe,req):
    if req in gizmo_dataframe:
        return True

    if req=="rad2d":
        gizmo_dataframe["rad2d"] = np.sqrt(gizmo_dataframe["Coordinates_x"]**2+gizmo_dataframe["Coordinates_y"]**2)
    elif req=="rad3d":
        gizmo_dataframe["rad3d"] = np.sqrt(gizmo_dataframe["Coordinates_x"]**2+gizmo_dataframe["Coordinates_y"]**2+gizmo_dataframe["Coordinates_z"]**2)
    elif req=="nH":
        gizmo_dataframe["nH"] = gizmo_dataframe["Density"]/(molecular_mass*proton_mass_cgs)
    elif req=="temp":
        gizmo_dataframe["temp"] = gamma_minus_one/boltzmann_cgs*(molecular_mass*proton_mass_cgs)*gizmo_dataframe["InternalEnergy"]
    elif req=="v_circ":
        load_calc_req(gizmo_dataframe,"rad2d")
        gizmo_dataframe["v_circ"] = -(gizmo_dataframe["Coordinates_x"]*gizmo_dataframe["Velocities_y"]-gizmo_dataframe["Coordinates_y"]*gizmo_dataframe["Velocities_x"])/gizmo_dataframe["rad2d"]
    elif req=="v_rad":
        load_calc_req(gizmo_dataframe,"rad3d")
        gizmo_dataframe["v_rad"] = (gizmo_dataframe["Coordinates_x"]*gizmo_dataframe["Velocities_x"]+
                                     gizmo_dataframe["Coordinates_y"]*gizmo_dataframe["Velocities_y"]+
                                     gizmo_dataframe["Coordinates_z"]*gizmo_dataframe["Velocities_z"]
                                     )/gizmo_dataframe["rad3d"]

    if req in gizmo_dataframe:
        return True
    return False

def calculate_phi(gizmo_dataframe):
    load_calc_reqs(gizmo_dataframe,["Coordinates_x","Coordinates_y","Coordinates_z","rad2d"])
#     if not check_requirements(gizmo_dataframe,["Coordinates_x","Coordinates_y","Coordinates_z"]):
#         raise Exception("Coordinates not loaded")
#     if not ("rad2d" in gizmo_dataframe):
#         gizmo_dataframe["rad2d"] = np.sqrt(gizmo_dataframe["Coordinates_x"]**2+gizmo_dataframe["Coordinates_y"]**2)
    gizmo_dataframe["phi"] = np.arctan(np.abs(gizmo_dataframe["Coordinates_z"]/gizmo_dataframe["rad2d"]))*180./np.pi

def calculate_vrad(gizmo_dataframe):
    load_calc_reqs(gizmo_dataframe,["Coordinates_x","Coordinates_y","Coordinates_z","Velocities_x","Velocities_y","Velocities_z","rad3d"])
#     if not check_requirements(gizmo_dataframe,["Coordinates_x","Coordinates_y","Coordinates_z","Velocities_x","Velocities_y","Velocities_z"]):
#         raise Exception("Coordinates not loaded")
#     if not ("rad3d" in gizmo_dataframe):
#         gizmo_dataframe["rad3d"] = np.sqrt(gizmo_dataframe["Coordinates_x"]**2+gizmo_dataframe["Coordinates_y"]**2+gizmo_dataframe["Coordinates_z"]**2)
    gizmo_dataframe["vrad"] =  ((gizmo_dataframe["Coordinates_x"]*gizmo_dataframe["Velocities_x"])+
                                (gizmo_dataframe["Coordinates_y"]*gizmo_dataframe["Velocities_y"])+
                                (gizmo_dataframe["Coordinates_z"]*gizmo_dataframe["Velocities_z"]))/gizmo_dataframe["rad3d"]

def load_value(f,gizmo_dataframe,value,internal_units=False):
    hdf5_path = "/PartType0/"+value
    indata = np.array(f[hdf5_path])
    if indata.ndim>3:
        raise Exception(">=3D array found in hdf5 datafile")
    if indata.ndim<1:
        raise Exception("<1D array found in hdf5 datafile??")
    if indata.ndim==1:
        gizmo_dataframe[value] = indata
        if not internal_units and value in unit_conversions:
            gizmo_dataframe[value]*=unit_conversions[value]
    else: # i.e. 2D array
        for i in range(indata.shape[1]):
            key = value+"_"+coord_suffixes[i]
            gizmo_dataframe[key] = indata[:,i]
            if not internal_units and value in unit_conversions:
                gizmo_dataframe[key]*=unit_conversions[value]


column_conversions = {
    "Coordinates_x":"x",
    "Coordinates_y":"y",
    "Coordinates_z":"z",
    "Velocities_x":"v\dx",
    "Velocities_y":"v\dy",
    "Velocities_z":"v\dz",
    "Masses":"particle mass",
    "InternalEnergy":"u",
    "SmoothingLength":"h",
}

def dump_ascii(filename,gizmo_dataframe):

    particles_adjusted = gizmo_dataframe\
                    .join(pd.DataFrame({"itype":np.zeros(len(gizmo_dataframe))}))\
                    .rename(columns=column_conversions)
    # TORUS reads:
    
    # unread line (title)
    # unread line
    # unread line (time labels)
    # # time, time_unit (overwritten in gadget/gizmo mode)
    # unread line
    # list of particle types. first 8 characters are junk, the rest are length 13 strings giving particle labels
    # # list of numbers of each particle
    # unread line ("units" title)
    # # list of units for each column, 16 characters per unit (`1.0000000E+00   ` is fine, gets overwritten for Gadget/GIZMO)
    # unread line
    # unread line
    # # labels of each column (must be same as number of units), 16 characters each
    output_string = """# splash style ascii dump
#
# time:             time unit ()
#   1.0000000E-03   1.0000000E+00
#
# npart:          gas  dark matter   boundary 1   boundary 2         star sink / black
#       {:13d}            0            0            0            0            0
# units:
""".format(len(particles_adjusted))\
        +"#  "+"1.0000000E+00   "*(len(particles_adjusted.keys())-1)+"\n"\
        +"#  "+" "*16*(len(particles_adjusted.keys()))+"\n"\
        +"#\n"\
        +"# "+"".join("{:16s}".format(k) for k in particles_adjusted.keys())+"\n"\
        +"  "+\
            particles_adjusted.to_csv(
                header=False,
                index=False,
                float_format='%15.7e')\
                  .replace('\n','\n  ')\
                  .replace(',',' ')

    with open(filename,'w') as f:
        f.write(output_string)
    

def load_gizmo_pandas(run_id,output_dir,snap_str,values,internal_units = False,gizmoDir=None):
    global unit_conversions
    
    if not internal_units:
        load_unit_conversions()
    
    if gizmoDir is None:
        gizmoDir = getGizmoDir(run_id)
    fullDir = gizmoDir+"/"+run_id+"/"+output_dir
    fullFile = fullDir+"/snapshot_"+snap_str+".hdf5"
    
    header = dict()

    f = h5py.File(fullFile,"r")
    
    file_header = f["/Header"]
    header["time"]=file_header.attrs.get("Time")
    if not internal_units:
        header["time"]*= 0.9778e9 # to yr
        header["time"]/=1.e6 # to Myr
    
    gizmo_dataframe = pd.DataFrame()
    
    for value in values:
        load_value(f,gizmo_dataframe,value,internal_units=internal_units)
            
    
    return header,gizmo_dataframe

def load_gizmo_nbody(run_id,output_dir,snap_str,load_vals=None,gizmoDir=None,load_binary_headers=False,only_header=False):
    if gizmoDir is None:
        gizmoDir = getGizmoDir(run_id)
    fullDir = gizmoDir+"/"+run_id+"/"+output_dir
    fullFile = fullDir+"/snapshot_"+snap_str+".hdf5"

    header = dict()

    f = h5py.File(fullFile,"r")
    file_header = f["/Header"]
    header["time"]=file_header.attrs.get("Time")
    header["time"]*= 0.9778e9 # to yr
    header["time"]/=1.e6 # to Myr
    
    if load_binary_headers:
        binary_header = f["/BH_binary"]
        for key,(label,unit_conversion,unit_string) in binary_headers_unit_name_conversions.items():
            value = binary_header.attrs.get(key)
            if value is not None:
                value = pynbody.array.SimArray(value*unit_conversion,unit_string)
            header[label] = value
    f.close()
    
    if only_header:
        return header
    
    snap = pynbody.load(fullFile)
    snap.set_units_system(velocity="km s**-1",mass="1e10 Msol",distance="kpc")

    if load_vals is not None:
        nbody_calc_vals(snap,load_vals)
    
    return header,snap

def nbody_calc_vals(snap,vals,**kwargs):
    """nbody_calc_vals
    
    wrapper for nbody_calc_val for when vals is an iterable
    """
    for val in vals:
        nbody_calc_val(snap,val,**kwargs)

def nbody_quickhist2d(snap,val1,val2,weights=None,run_name="quickhist",cmap='viridis',log=False,logx=False,logy=False,**kwargs):
    fig,sp = plt.subplots()
    scale_x_func = np.log10 if logx else lambda x:x
    scale_y_func = np.log10 if logy else lambda x:x
    if type(val1)==str:
        a1 = snap[val1]
        sp.set_xlabel(val1)
    else:
        a1 = val1
    if type(val2)==str:
        a2 = snap[val2]
        sp.set_ylabel(val2)
    else:
        a2 = val2
    
    v1 = scale_x_func(a1)
    v2 = scale_y_func(a2)
    
    nonanslice = (np.isfinite(v1)) & (np.isfinite(v2))
    
    suffix=""
    
    if weights is None:
        hist,xs,ys = np.histogram2d(v1[nonanslice],v2[nonanslice],bins=(100,100),**kwargs)
    else:
        if type(weights)==str:
            w = snap[weights]
            suffix="_"+weights
        else:
            w = weights
        whist,xs,ys = np.histogram2d(v1[nonanslice],v2[nonanslice],weights=w,bins=(100,100),**kwargs)
        nhist,xs,ys = np.histogram2d(v1[nonanslice],v2[nonanslice],bins=(100,100),**kwargs)
        hist=whist/nhist
    if log:
        hist = np.log10(hist)
    mappable = sp.pcolormesh(xs,ys,hist.T,cmap=cmap)
    fig.colorbar(mappable)
    fig.savefig(f"../figures/hist2d_{run_name}_{val1}_{val2}{suffix}.png")
    plt.close()


def nbody_calc_val(snap
                   ,val
                   ,m_smbh = 1.e6
                   ,smbh_soft = 0.01e-3
                   ,m_hernquist = 1.e9
                   ,a_hernquist = 250.e-3):

    G_msun_kpcyr2_kpc2 = 4.49972292e-24
#     G_msun_kms2_kpc2 = 4.49972292e-24
    G_msun_km2s2_kpc = 4.3022682e-6

    def accel_smbh_hernquist(r):
        """accel_smbh_hernquist

        Calculates acceleration due to a Plummer softend SMBH plus a Hernquist bulge.

        Mass units are in Msun, length units are in kpc, time units are in years
        """
        x = r/a_hernquist
        r2 = r**2
        m =   m_smbh*r2/(r2+smbh_soft**2) + m_hernquist * (x/(1.+x))**2
        accel = -(G_msun_kpcyr2_kpc2 * m/r2)
        return pynbody.array.SimArray(accel,"kpc yr**-2")

    def grav_accel(snap):
        if "grav_accel" not in snap:
            snap["grav_accel"]=accel_smbh_hernquist(snap["r"])

    def t_dyn(snap):
        grav_accel(snap)
        if "tdyn" not in snap:
            snap["tdyn"]=(snap["r"]/-snap["grav_accel"])**(1,2)
 
    def phi(snap):
        if "phi" not in snap:
            snap["phi"]=pynbody.array.SimArray(-G_msun_km2s2_kpc*(
                        m_hernquist/(a_hernquist+snap["r"]) +
                        m_smbh/np.sqrt(snap["r"]**2+smbh_soft**2)
                        )
                        ,"km**2 s**-2")

    def v_esc(snap):
        phi(snap)
        if "vesc" not in snap:
            snap["vesc"]= np.sqrt(-2*snap["phi"])
    
    if val=='tdyn':
        t_dyn(snap)
    elif val=='phi':
        phi(snap)
    elif val=='vesc':
        v_esc(snap)
    elif val=='grav_accel':
        grav_accel(snap)
    elif val=='nH':
        snap["nH"] = snap.gas['rho'].in_units("g cm**-3")/(molecular_mass*pynbody.array.SimArray(proton_mass_cgs,'g'))
    elif val=="temp":
        snap["temp"] = (gamma_minus_one/pynbody.array.SimArray(boltzmann_cgs,"erg K**-1")*(molecular_mass*pynbody.array.SimArray(proton_mass_cgs,'g'))*snap["u"]).in_units('K') 
    

def gridsize_from_n(n,aspect=1.):
    nx = 1
    ny = 1
    
    while nx*ny<n:
        if nx>aspect*ny:
            ny+=1
        else:
            nx+=1
    return nx,ny

def load_run_parameters(rundir,dir=""):
    with open(f"{dir}data/runprams_{rundir}.pkl",'rb') as f:
        run_parameters = pickle.load(f) 

    for run_name,run_prams in run_parameters.items():
        if run_prams['outflowVel']<50.:
            run_prams['marker'] = markers[0]
        elif run_prams['outflowVel']<200.:
            run_prams['marker'] = markers[1]
        else:
            run_prams['marker'] = markers[2]
        run_prams['color'] = cmap((run_prams['outflowThetaCentre']-20.)/50.)

        run_prams['size'] = 24.*(1.+2.*np.log10(run_prams["outflowRate"]/0.0378))
        run_prams['thickness'] = (run_prams["outflowThetaWidth"]-20.)/30.*1. + 1.
        
    return run_parameters

def run_parameters_names(run_parameters):
    pairs = sorted(run_parameters.items(), key=lambda x: (-x[1]['outflowRate'],x[1]['outflowVel'],x[1]['outflowThetaCentre']))
    ordered_keys = [x[0] for x in pairs]

    for run_name,run_prams in run_parameters.items():
        if run_prams['outflowVel']<50.:
            vel_str = "slow"
        elif run_prams['outflowVel']<200.:
            vel_str = "vesc"
        else:
            vel_str = "rapid"
        theta_str = "{:02d}".format(int(90-run_prams['outflowThetaCentre']))
        mass_str = "heavy" if run_prams["outflowRate"]>0.1 else "light"
        thick_str = "thick" if run_prams["outflowThetaWidth"]>40. else "thin"
        
        run_prams["name"]="_".join([mass_str,vel_str,theta_str,thick_str])

#     thin_index = 0
#     thick_index = 0
#     all_names = []
#     
#     for key in ordered_keys:
#         run_prams = run_parameters[key]
#         if run_prams["outflowRate"]>0.1:
#             run_prams["name"] = string.ascii_uppercase[thick_index]
#             thick_index+=1
#         else:
#             run_prams["name"] = string.ascii_lowercase[thin_index]
#             thin_index+=1
# #         all_names+=run_prams["name"]
# #     sorted_names,sorted_keys = zip(*sorted(zip(all_names,ordered_keys)))   
    return ordered_keys

def run_parameters_angles(run_parameters):
    """Extra size of outflow, "radiance" etc - uses phi=0° is disk plane coordinates, while outflowThetaCentre uses theta=0° is polar"""
    for run_prams in run_parameters.values():
        run_prams["outflowPhiTop"]=90.-run_prams["outflowThetaCentre"]-run_prams["outflowThetaWidth"]/2.
        run_prams["outflowPhiBottom"]=90.-run_prams["outflowThetaCentre"]+run_prams["outflowThetaWidth"]/2.
        run_prams["openingAngle"]=run_prams["outflowPhiTop"]*2
        run_prams["coveringFraction"]=np.sin(np.radians(run_prams["outflowPhiBottom"]))-np.sin(np.radians(run_prams["outflowPhiTop"]))
        run_prams["coveringSolidAngle"]=run_prams["coveringFraction"]*4*np.pi

def run_parameters_table(run_parameters):
    sorted_keys=run_parameters_names(run_parameters)
    outp = ""
#     for run_name,run_prams in run_parameters.items():
    for key in sorted_keys:
        run_prams = run_parameters[key]
        outp+=  "{} & ${:6.4f}$ & ${:3.0f}$ & ${:2.0f}^\circ$ & ${:2.0f}^\circ$\\\\\n".format(
#                 run_prams["name"],
                run_prams["name"].replace("_","\_"),
                run_prams["outflowRate"],
                run_prams["outflowVel"],
                90-run_prams["outflowThetaCentre"],
                run_prams["outflowThetaWidth"]
        )
    return outp,sorted_keys

class cloudy_table:
    """cloudy_table - a class to interface with the dust_temp_interp c++ libraries extracting from the GIZMO code.
    Load in a set of tables with specified date and resolution strings, and with optional path, by instantiating this class.
    Then pass particles dataframe to `interp` (after loading nH,temp,AGNIntensity,AGNDepth) and it returns 
    """
    struct_attributes = [
             'arad',
             'column_out',
             'dCool',
             'dHeat',
             'dg',
             'dustT',
             'line_co1',
             'line_co2',
             'line_h2_1',
             'line_h2_2',
             'line_h2_3',
             'line_hcn1',
             'line_hcn2',
             'line_12mic',
             'line_8mic',
             'line_850mic',
             'line_hcn2',
             'line_hcn2',
             'opac_abs',
             'opac_scat'
             ]

    def double_pointer_to_array(self,x,n):
        """Turns C style arrays of known size into numpy arrays.
        For use with struct of arrays returned by swig.
        """
        ptr = int(x)
        type_size = ctypes.c_double * n
        base_array = type_size.from_address(ptr)
        numpy_array = np.array(base_array)
        return numpy_array


    def __init__(self,tableDate="281118",tableRes="0.0001",prefix=""):
        self.load_table(tableDate,tableRes,prefix)

    def load_table(self,tableDate="281118",tableRes="0.0001",prefix="../coolheat_tab_marta"):
        self.chTab = tab_interp.CoolHeatTab( (prefix+"shrunk_table_labels_"+tableDate+"tau.dat"),
                                        (prefix+"shrunk_table_"+tableDate+"_m"+tableRes+"_hsmooth_tau.dat"),
                                        (prefix+"shrunk_table_labels_"+tableDate+"taunodust.dat"),
                                        (prefix+"shrunk_table_"+tableDate+"_m"+tableRes+"_hsmooth_taunodust.dat"),
                                        (prefix+"shrunk_table_labels_"+tableDate+"taudense.dat"),
                                        (prefix+"shrunk_table_"+tableDate+"_m"+tableRes+"_hsmooth_taudense.dat")
                                        )
 
#    def load_table(self,tableDate="281118",tableRes="0.0001",prefix=""):
#        self.chTab = tab_interp.CoolHeatTab( (prefix+"coolheat_tab_marta/shrunk_table_labels_"+tableDate+"tau.dat"),
#                                        (prefix+"coolheat_tab_marta/shrunk_table_"+tableDate+"_m"+tableRes+"_hsmooth_tau.dat"),
#                                        (prefix+"coolheat_tab_marta/shrunk_table_labels_"+tableDate+"taunodust.dat"),
#                                        (prefix+"coolheat_tab_marta/shrunk_table_"+tableDate+"_m"+tableRes+"_hsmooth_taunodust.dat"),
#                                        (prefix+"coolheat_tab_marta/shrunk_table_labels_"+tableDate+"taudense.dat"),
#                                        (prefix+"coolheat_tab_marta/shrunk_table_"+tableDate+"_m"+tableRes+"_hsmooth_taudense.dat")
#                                        )
        
    
    def interp(self,particles):
        tabStructs=self.chTab.interpTabArray( 
                       particles["nH"].astype(np.float64),
                       particles["temp"].astype(np.float64),
                       particles["AGNIntensity"].astype(np.float64),
                       particles["AGNDepth"].astype(np.float64)
                       )
        
        n = len(particles)
        for attr in self.struct_attributes:
            particles[attr] = self.double_pointer_to_array(tabStructs.__getattr__(attr),n)


if __name__=='__main__':
    header,snap = load_gizmo_nbody("binary_ecc0_HPM_MM_radoff","binary_ecc0","120",gizmoDir="/export/2/lb1g19/data/",load_binary_headers=True)
#     run_table = load_run_parameters("3032")
#     run_parameters_angles(run_table)
#     run_parameters_names(run_table)


#     ordered_keys = run_parameters_names(run_table)
#     print([(x,run_table[x]['name']) for x in ordered_keys])
#     print(run_parameters_table(load_run_parameters("3032"))[0])
