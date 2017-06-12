print("Importing")

# Import libraries to do our magic
import numpy as np
import h5py
import sys

from sys import path
path.append("src/")
import tab_interp

print("Loading table (short form)")
chTab = tab_interp.CoolHeatTab(("coolheat_tab_marta/shrunk_table_labels_080517.dat"),("coolheat_tab_marta/shrunk_table_080517.dat"))
interpTabVec = np.vectorize(chTab.interpTab)

print("Running")


run_id = sys.argv[1]
output_dir = sys.argv[2]
snap_str = sys.argv[3]

luminosity = 6.2849e42 # erg/s, DEPENDS ON EACH RUN!


#f = h5py.File("/export/1/djw/gizmo_public/disc_sf_out/snapshot_"+snap_str+".hdf5","r")
#f = h5py.File("/export/1/djw/gizmos/1003/rearranged_test/snapshot_"+snap_str+".hdf5","r")
#f = h5py.File("/export/1/djw/gizmos/1003/control_test/snapshot_"+snap_str+".hdf5","r")
#f = h5py.File("/export/1/djw/gizmos/1003/tree_tests/snapshot_"+snap_str+".hdf5","r")
#f = h5py.File("/export/1/djw/gizmos/1003/isotropic_test/snapshot_"+snap_str+".hdf5","r")
f = h5py.File("/export/1/djw/gizmos/"+run_id+"/"+output_dir+"/snapshot_"+snap_str+".hdf5","r")

xyz = np.array(f["/PartType0/Coordinates"])
rad_p = np.sqrt(xyz[:,0]**2+xyz[:,1]**2+xyz[:,2]**2)
rad2d_p = np.sqrt(xyz[:,0]**2+xyz[:,1]**2)
z_p = np.abs(xyz[:,2])

N_part = rad_p.size

u_p = np.array(f["/PartType0/InternalEnergy"])
rho_p = np.array(f["/PartType0/Density"])
vel_p = np.array(f["/PartType0/Velocities"])
a_p = np.array(f["/PartType0/Acceleration"])
h_p = np.array(f["/PartType0/SmoothingLength"])
m_p = np.array(f["/PartType0/Masses"])
dt_p = np.array(f["/PartType0/TimeStep"])

id_p = np.array(f["/PartType0/ParticleIDs"])

if "/PartType0/AGNHeat" in f:
    agn_heat_p = np.array(f["/PartType0/AGNHeat"])
else:
    agn_heat_p = np.zeros(N_part)

if "/PartType0/RadiativeAcceleration" in f:
    radaccel_p = np.array(f["/PartType0/RadiativeAcceleration"])
else:
    radaccel_p = np.zeros((N_part,3))

# if "/PartType0/AGNDepth" in f:
#     depth_p = np.array(f["/PartType0/AGNDepth"]) # unitless
# else:
#     depth_p = np.zeros(N_part)

if "/PartType0/AGNColDens" in f:
    depth_p = np.array(f["/PartType0/AGNColDens"]) # surface density units
else:
    depth_p = np.zeros(N_part)

#print(depth_p)

agn_heat_p*=1e10/3.08568e+16# to erg/s/g

depth_p*=(1.989e+43/3.086e+21**2) # to g/cm**2

rho_p*=6.77e-22 # to g/cm**3
u_p*=1.e10 # to erg/g
rad_p*=1. # already in kpc
vel_p*=1. # in km/s
h_p*=1. # already in kpc
m_p*=1.989e+43 # 10^10 solar masses to g
radaccel_p*=3.0857e21/3.08568e+16**2 # to cm/s/s
a_p*=3.0857e21/3.08568e+16**2 # to cm/s/s
#dt_p*=3.08568e+16 # to seconds

G = 6.67259e-8 # in cgs

#Gintern = 7.e-45 # Gintern*(mass in g)/(distance in kpc)**2 = acceleration in cm/s/s

cs_p = np.sqrt(u_p) # cm/s
#a_p*=# 

vrad_p = (xyz[:,0]*vel_p[:,0]+xyz[:,1]*vel_p[:,1]+xyz[:,2]*vel_p[:,2])/rad_p
vcirc_p = (xyz[:,1]*vel_p[:,0]-xyz[:,0]*vel_p[:,1])/rad2d_p

omega_p = vcirc_p/(rad_p * 3.08567758e16) # in Hz

# in km/s/Gyr, roughly
arad_p = (xyz[:,0]*a_p[:,0]+xyz[:,1]*a_p[:,1]+xyz[:,2]*a_p[:,2])/rad_p
radrad_p = (xyz[:,0]*radaccel_p[:,0]+xyz[:,1]*radaccel_p[:,1]+xyz[:,2]*radaccel_p[:,2])/rad_p
acirc_p = (xyz[:,1]*a_p[:,0]-xyz[:,0]*a_p[:,1])/rad2d_p

lambda_cool = 1.e-24 # for testing

molecular_mass = 4./(1.+3.*.76)
proton_mass_cgs = 1.6726e-24
gamma_minus_one = 5./3.-1.
boltzmann_cgs = 1.38066e-16

depth_p/=(molecular_mass*proton_mass_cgs) # N in cm**(-2)

mJ_p = np.pi**(2.5)*cs_p**3*G**(-1.5)*rho_p**(-.5)/6.

nH_p = rho_p/(molecular_mass*proton_mass_cgs)

TK_p = gamma_minus_one/boltzmann_cgs*(molecular_mass*proton_mass_cgs)*u_p

cs_p = np.sqrt(10./9.*u_p) # in cm/s

rcool_p = lambda_cool*nH_p/(molecular_mass*proton_mass_cgs)

tcool_p = u_p/rcool_p

tcool_p/=3.154e+7 # convert from seconds to years

if "/PartType0/AGNIntensity" in f:
    flux_p = np.array(f["/PartType0/AGNIntensity"]) # energy per surface area per time
    flux_p*=1.989e+53/(3.086e21)**2/(3.08568e+16)
else:
    flux_p = luminosity/(4.*np.pi*rad_p**2)
    flux_p/=(3.086e21)**2 # erg/s/kpc**2 to erg/s/cm**2


tabStructs = interpTabVec(nH_p.astype(np.float64),TK_p.astype(np.float64),flux_p.astype(np.float64),depth_p.astype(np.float64))
dustTemp = map(lambda y: y.dustT, tabStructs)
dustTemp = np.array(dustTemp)

# dustTemp = np.zeros(N_part)
# #for ii in range(N_part):
# for ii in [0]:
#     x=chTab.interpTab(nH_p.astype(np.float64)[ii],TK_p.astype(np.float64)[ii],flux_p.astype(np.float64)[ii],depth_p.astype(np.float64)[ii])
#     #dustTemp[ii] = chTab.interpTab(nH_p.astype(np.float64)[ii],TK_p.astype(np.float64)[ii],flux_p.astype(np.float64)[ii],depth_p.astype(np.float64)[ii]).dustT
#     dustTemp[ii] = x.dustT
#     if ( dustTemp[ii]>1.e9 or dustTemp[ii]<1.e-10):
#         print("bad temp",ii,dustTemp[ii])
#         print(nH_p.astype(np.float64)[ii],TK_p.astype(np.float64)[ii],flux_p.astype(np.float64)[ii],depth_p.astype(np.float64)[ii])
#         #x=chTab.interpTab(nH_p.astype(np.float64)[ii],TK_p.astype(np.float64)[ii],flux_p.astype(np.float64)[ii],depth_p.astype(np.float64)[ii])
#         print(x.dustT,x.dHeat,x.dCool)
#         sys.exit()
# 
#print(np.log10([np.max(dustTemp),np.min(dustTemp)]))

#dout = [rad_p,z_p,rho_p,nH_p,TK_p,u_p,rcool_p,tcool_p,vrad_p,vcirc_p,arad_p,acirc_p,xyz[:,0],xyz[:,1],vel_p[:,0],vel_p[:,1],a_p[:,0],a_p[:,1]]

print("binning")

Nbins = 50
rmax = .02

#tobin = [m_p,omega_p]
#tobin = np.array(tobin)

bins = np.linspace(0,rmax,num=Nbins)
bindices = np.digitize(rad2d_p,bins)
#binmeans = np.array([np.nanmean(dout[:,bindices==bindex],axis=1) for bindex in range(bins.size)])
#binmeans = np.zeros((tobin.shape[0],bins.size))
#for bindex in range(bins.size):
    #binmeans[bindex] = np.mean(rad_p[bindices==bindex])
#    binmeans[:,bindex] = np.nanmean(tobin[:,bindices==bindex],axis=1)

pmass = m_p[0] # assume all gas!

binmasses = np.zeros((bins.size))
binvdisp =  np.zeros((bins.size))
for bindex in range(bins.size):
    if ( np.sum(bindices==bindex)>0 ):
        binmasses[bindex] = np.nansum(bindices==bindex)
        binvdisp[bindex] = np.nanstd(vrad_p[bindices==bindex])

binarea = np.zeros(Nbins)

binarea = np.pi*(bins[1:]**2-bins[:-1]**2)*1e6*(3.086e18)**2 # in cm**2

binarea = np.append(binarea,binarea[-1]*1.e5) # large area -> small surface density goes in denominator -> large Q out
binvdisp = np.append(binvdisp,1.e20)
binvdisp*=1.e5 # from km/s to cm/s

binsurf = binmasses/binarea

binsurf = np.insert(binsurf, 0, 1.e-20)

binstep = bins[1]-bins[0]
binmids = bins+binstep/2.

surfs = np.interp(rad2d_p,binmids,binsurf[1:])
vdisp_p = np.interp(rad2d_p,binmids,binvdisp[1:])

#vdisp_p = binvdisp[bindices]

vdisp_p = np.sqrt(vdisp_p**2+cs_p**2)

Q_approx_p = vdisp_p*omega_p/(pmass*surfs*np.pi*G)

print("packing")

dout = [rad_p,rad2d_p,xyz[:,0],xyz[:,1],xyz[:,2],vel_p[:,0],vel_p[:,1],vel_p[:,2],nH_p,TK_p,
        agn_heat_p,depth_p,dustTemp,flux_p,dt_p,h_p,u_p,m_p,mJ_p,cs_p,
        radaccel_p[:,0],radaccel_p[:,1],radaccel_p[:,2],radrad_p,a_p[:,0],a_p[:,1],a_p[:,2],arad_p,cs_p,omega_p,
        surfs,vdisp_p,Q_approx_p]
#dout = [rad_p,rad2d_p,xyz[:,0],xyz[:,1],xyz[:,2],rho_p,TK_p,agn_heat_p,depth_p]

#thinslice = (np.abs(xyz[:,2])<1.) & (rad2d_p<.02)
#dout = np.array(dout)[:,thinslice]

dout = np.array(dout)

#print(dout[:,0])
#print(dout[0,:])

#binmeans = np.mean(dout[



#dout = [u_p,m_p]
dout = dout.T
#np.savetxt("data/rhT.dat",dout)

#dout.sort(0)

print("Dumping")

np.savetxt("data/allp"+run_id+output_dir+"_"+snap_str,dout)

# binmeans = binmeans.T
# np.savetxt("data/prof"+run_id+output_dir+"_"+snap_str,binmeans)

#print(m_p[0])
