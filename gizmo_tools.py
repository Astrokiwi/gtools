import socket
import getpass
import os
import numpy as np
import re
import pandas as pd
import h5py

import pynbody

# for interpolating GIZMO tables
import ctypes
from sys import path
path.append("src/")
import tab_interp

molecular_mass = 4./(1.+3.*.76)
proton_mass_cgs = 1.6726e-24
gamma_minus_one = 5./3.-1.
boltzmann_cgs = 1.38066e-16
msun_g = 1.989e33

unit_conversions = None
unit_names = None

coord_suffixes = "xyz"

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

    if uname=='lb1g19':
        if sname=="trillian":
            return "/export/2/lb1g19"

    raise Exception("Unknown server/username; add server, username and directory to gizmo_tools.py")

def getMovieDir():
    sname = socket.gethostname()
    uname = getpass.getuser()

    if uname=='djw1g16' or uname=='djw':
        if ( sname=="trillian" ):
            return "/export/1/djw/movies"
        elif ( sname=="srv01921" ):
            return "/srv/djw1g16/movies"

    if uname=='lb1g19':
        if sname=="trillian":
            return "/export/2/lb1g19"
    
    raise Exception("Unknown server; add server and directory to gizmo_tools.py")

def lastConsecutiveSnapshot(run_id,output_dir,dumpsOrdered=True):
    """When models are rerun, the snapshot file with the largest number (e.g. snapshot_100.dat)
    may be from a previous model. So we want the last snapshot that was made *after* the previous
    snapshot. However, this only works if the snapshots haven't been copied, because this kills
    the "last modified" data. So we can set dumpsOrdered=False and force us to just take the last
    snapshot numerically rather than by time"""

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

def load_gizmo_pandas(run_id,output_dir,snap_str,values,internal_units = False):
    global unit_conversions
    
    if not internal_units:
        load_unit_conversions()
    
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

def load_gizmo_nbody(run_id,output_dir,snap_str):
    gizmoDir = getGizmoDir(run_id)
    fullDir = gizmoDir+"/"+run_id+"/"+output_dir
    fullFile = fullDir+"/snapshot_"+snap_str+".hdf5"

    header = dict()

    f = h5py.File(fullFile,"r")
    file_header = f["/Header"]
    header["time"]=file_header.attrs.get("Time")
    header["time"]*= 0.9778e9 # to yr
    header["time"]/=1.e6 # to Myr
    f.close()
    
    snap = pynbody.load(fullFile)
    snap.set_units_system(velocity="km s**-1",mass="1e10 Msol",distance="kpc")

    
    return header,snap
    

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

    def load_table(self,tableDate="281118",tableRes="0.0001",prefix=""):
        self.chTab = tab_interp.CoolHeatTab( (prefix+"coolheat_tab_marta/shrunk_table_labels_"+tableDate+"tau.dat"),
                                        (prefix+"coolheat_tab_marta/shrunk_table_"+tableDate+"_m"+tableRes+"_hsmooth_tau.dat"),
                                        (prefix+"coolheat_tab_marta/shrunk_table_labels_"+tableDate+"taunodust.dat"),
                                        (prefix+"coolheat_tab_marta/shrunk_table_"+tableDate+"_m"+tableRes+"_hsmooth_taunodust.dat"),
                                        (prefix+"coolheat_tab_marta/shrunk_table_labels_"+tableDate+"taudense.dat"),
                                        (prefix+"coolheat_tab_marta/shrunk_table_"+tableDate+"_m"+tableRes+"_hsmooth_taudense.dat")
                                        )
        
    
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


