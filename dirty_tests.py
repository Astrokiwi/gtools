import gizmo_tools
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate


m_smbh = 1.e6
smbh_soft = 0.01e-3
m_hernquist = 1.e9
a_hernquist = 250.e-3

G_msun_kpcyr2_kpc2 = 4.49972292e-24
G_msun_km2s2_kpc = 4.3022682e-6

def calc_pot(x):
   return -G_msun_km2s2_kpc*(
                        m_hernquist/(a_hernquist+x) +
                        m_smbh/np.sqrt(x**2+smbh_soft**2)) # in km**2/s**2

def accel_smbh_hernquist(r):
    """accel_smbh_hernquist

    Calculates acceleration due to a Plummer softend SMBH plus a Hernquist bulge.

    Mass units are in Msun, length units are in pc, time units are in years
    """
    x = r/a_hernquist
    r2 = r**2
    m =   m_smbh*r2/(r2+smbh_soft**2) + m_hernquist * (x/(1.+x))**2
    accel = -(G_msun_kpcyr2_kpc2 * m/r2)
    return accel

r = np.logspace(-5,5,1000)
pot = calc_pot(r)

r_max_interp = interpolate.interp1d(pot,r,bounds_error=False)

launch_energy = (100.**2)/2.+calc_pot(1.e-3)
r_max = r_max_interp(launch_energy)

print(launch_energy,"kmkm/ss",r_max,"kpc")

tdyn_max = np.sqrt(r_max/-accel_smbh_hernquist(r_max))

print(tdyn_max,"yr")

# snap_strs = ["010","050","100","150","200"]
snap_strs = ["010"]
# 
for snap_str in snap_strs:
    print(snap_str)
    header,snap = gizmo_tools.load_gizmo_nbody("3032","longrun_medflow_vesc_defaultaniso_polar",snap_str,load_vals=["temp"])                                          
    snap['v_norm'] = snap['v2']**(1,2)
    snap['cs'].convert_units("km s**-1")
    snap['mach'] = snap['v_norm']/snap['cs']
    gizmo_tools.nbody_quickhist2d(snap,'v_norm','cs',log=True,logx=False,logy=False,run_name=f"dirty{snap_str}")
    gizmo_tools.nbody_quickhist2d(snap,'r','cs',log=True,logx=True,logy=False,run_name=f"dirty{snap_str}")
    gizmo_tools.nbody_quickhist2d(snap,'r','mach',log=True,logx=True,logy=False,run_name=f"dirty{snap_str}")
    print(snap['cs'])


#     snap['r_max']=r_max_interp(snap['te'])
#     snap['mte']=-snap['te']
# 
# #     gizmo_tools.nbody_quickhist2d(snap,'r','vesc',log=True,logx=False,logy=False,run_name=f"dirty{snap_str}")
# #     gizmo_tools.nbody_quickhist2d(snap,'r','r_max',log=True,logx=True,logy=True,run_name=f"dirty{snap_str}")
# #     gizmo_tools.nbody_quickhist2d(snap,'r','mte',log=True,logx=True,logy=True,run_name=f"dirty{snap_str}")
#     gizmo_tools.nbody_quickhist2d(snap,'r','te',log=True,logx=False,logy=False,run_name=f"dirty{snap_str}")
#     gizmo_tools.nbody_quickhist2d(snap,'r','phi',log=True,logx=False,logy=False,run_name=f"dirty{snap_str}")
#     gizmo_tools.nbody_quickhist2d(snap,'r','ke',log=True,logx=False,logy=False,run_name=f"dirty{snap_str}")
# #     gizmo_tools.nbody_quickhist2d(snap,'r','te',log=True,logx=False,logy=True,run_name=f"dirty{snap_str}")
# #     gizmo_tools.nbody_quickhist2d(snap,'r_max','mte',log=True,logx=True,logy=True,run_name=f"dirty{snap_str}")
# 
# #     snap["tdyn_max"] = (snap["r_max"]/-snap["grav_accel"])**(1,2)
# #     gizmo_tools.nbody_quickhist2d(snap,'r','tdyn_max',log=True,logx=True,logy=True,run_name=f"dirty{snap_str}")
