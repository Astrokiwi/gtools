import sys
import numpy as np
import h5py
import gizmo_tools
# import matplotlib as mpl
# mpl.use('Agg')
import matplotlib.pyplot as P


run_id = sys.argv[1]
output_dir = sys.argv[2]
snap_id = sys.argv[3]

gizmoDir = gizmo_tools.getGizmoDir(run_id)
fullDir = gizmoDir+"/"+run_id+"/"+output_dir
with h5py.File(fullDir+"/snapshot_"+snap_id+".hdf5","r") as f:
#     v = np.array(f["/PartType0/TimeStep"])
#     v = np.array(f["/PartType0/RadiativeAcceleration"])
#     v = np.array(f["/PartType0/AGNHeat"])
    v2 = np.array(f["/PartType0/TimeStep"])
    v = np.array(f["/PartType0/InternalEnergy"])
    v = v[v2>1.e-9]
    
    vsort = np.sort(v)
    print(vsort[-20:])
    print(vsort[:20])
    print(np.min(v),np.mean(v),np.median(v),np.max(v))
    
#     np.savetxt("data/quickdump.dat",[v,v2])
#     heatTime = v/v2
#     heatTimeSorted = np.sort(heatTime)
#     print(heatTimeSorted[140179],heatTimeSorted[-140179])
#     print(np.sum(heatTimeSorted<0.),np.sum(heatTimeSorted>0.),np.sum(np.abs(heatTimeSorted)<1.e-9))

    
#     heat = np.array(f["/PartType0/AGNHeat"])
#     heat_sort = np.sort(heat)
#     print(heat_sort[-20:],np.max(heat),np.min(heat),np.sum(heat!=0.))
#     m_p = np.array(f["/PartType0/Masses"])
#     xyz = np.array(f["/PartType0/Coordinates"])
#     rn_p = np.sqrt(xyz[:,0]**2+xyz[:,1]**2+xyz[:,2]**2) # in kpc
#     v = np.array(f["/PartType0/Velocities"])
#     vn_p = np.sqrt(v[:,0]**2+v[:,1]**2+v[:,2]**2) # in km/s
#     k_en = np.sum(m_p*vn_p**2/2) # in solar masses km**2/s**2
#     
#     m_bh = 1.e5
#     c_bh = 1.e-5 # 0.01 pc
#     pot_en = np.sum(m_p*(1./c_bh*np.arctan(rn_p/c_bh)-np.pi/2/c_bh-rn_p/(rn_p**2+c_bh**2)))
#     pot_en*=m_bh # solar masses**2/kpc
#     pot_en*=4.3022682e-6 # *G and unit conversions to solar masses km**2/s**2
#     e_tot = k_en+pot_en
#     print("K={0:+.6f} U={1:+.6f} E={2:+.6f}".format(k_en,pot_en,e_tot))

#     h_p = np.array(f["/PartType0/SmoothingLength"]) # kpc
#     m_p = np.array(f["/PartType0/Masses"]) # 10^10 msun
#     rho_p = np.array(f["/PartType0/Density"]) # 10^10 msun/kpc**3
#     h_p*=1.e3
#     m_p*=1.e10
#     rho_p*=6.77e-22 # to g/cm**2
#     molecular_mass = 4./(1.+3.*.76)
#     proton_mass_cgs = 1.6726e-24
#     rho_p/=(molecular_mass*proton_mass_cgs) # to particles/cm**2
#     np.savetxt("data/mhr{}.dat".format(output_dir),np.array([m_p,h_p,rho_p]).T)



#     h_p = np.array(f["/PartType0/SmoothingLength"])
#     h_sort = np.sort(h_p)
#     print(h_p[-20:])

#     xyz = np.array(f["/PartType0/Coordinates"])
#     rad2d_p = np.sqrt(xyz[:,0]**2+xyz[:,1]**2)
#     id_p = np.array(f["/PartType0/ParticleIDs"]).astype(int)
#     smallrad = (rad2d_p<2.e-4)
#     print(rad2d_p[smallrad])
#     print(id_p[smallrad])

#     xyz = np.array(f["/PartType0/Coordinates"])
#     rad3d_p = np.sqrt(xyz[:,0]**2+xyz[:,1]**2+xyz[:,2]**2) # in kpc
#     v = np.array(f["/PartType0/Velocities"])
#     v3d_p = np.sqrt(v[:,0]**2+v[:,1]**2+v[:,2]**2) # in km/s
#     a = np.array(f["/PartType0/Acceleration"])
#     a3d_p = np.sqrt(a[:,0]**2+a[:,1]**2+a[:,2]**2)
#     ra = np.array(f["/PartType0/RadiativeAcceleration"])
#     ra3d_p = np.sqrt(ra[:,0]**2+ra[:,1]**2+ra[:,2]**2) # in km/s
#     np.savetxt("data/rv{}.dat".format(output_dir),np.array([rad3d_p,v3d_p,a3d_p,ra3d_p]).T)
# 

#     ira = np.array(f["/PartType0/IRRadAccel"])
#     ira3d_p = np.sqrt(ira[:,0]**2+ira[:,1]**2+ira[:,2]**2) # in km/s
#     np.savetxt("data/rv{}.dat".format(output_dir),np.array([rad3d_p,v3d_p,a3d_p,ra3d_p,ira3d_p]).T)

#     outp = np.vstack((xyz.T,ra.T,ira.T))
#     np.savetxt("data/raira{}.dat".format(output_dir),outp.T)

#     #rad3d_p*=1.e3 # to pc
#     print(np.sum(rad3d_p>0.05))

#     v = np.array(f["/PartType0/TrueNumberOfNeighbours"])
#     v2 = np.array(f["/PartType0/Density"])
#     v3 = np.array(f["/PartType0/SmoothingLength"])
#     np.savetxt("data/vv2.dat",np.array([v,v2,v3]).T)
#     

#     v = np.array(f["/PartType0/InternalEnergy"])

#     vsort = np.sort(v)
#     print(vsort[-20:])
#     print(vsort[:20])
#     print(np.min(v),np.mean(v),np.median(v),np.max(v))
