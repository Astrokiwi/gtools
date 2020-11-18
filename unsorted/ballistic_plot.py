#import pynbody
import gizmo_tools
import matplotlib.pyplot as plt
import numpy as np
import pynbody
import functools

m_smbh = 1.e6
smbh_soft = 0.01
m_hernquist = 1.e9
a_hernquist = 250.


def accel_smbh_hernquist(m_smbh,smbh_soft,m_hernquist,a_hernquist,r):
    """accel_smbh_hernquist

    Calculates acceleration due to a Plummer softend SMBH plus a Hernquist bulge.

    Mass units are in Msun, length units are in pc, time units are in years
    """
#     G_msun_km2s2_pc = 0.004302052
    G_msun_pcyr2_pc2 = 4.5e-15
    x = r/a_hernquist
    r2 = r**2
    m =   m_smbh*r2/(r2+smbh_soft**2) + m_hernquist * (x/(1.+x))**2
    accel = -(G_msun_pcyr2_pc2 * m/r2)
    return accel

accel = functools.partial(accel_smbh_hernquist,m_smbh,smbh_soft,m_hernquist,a_hernquist)

def t_dyn(r):
    return np.sqrt(r/-accel(r))

run_id = "3032"
run_name = "longrun_medflow_vesc_defaultaniso_polar_highedd"
snap_id = "205"

header,snap = gizmo_tools.load_gizmo_nbody(run_id,run_name,snap_id)
snap['t_dyn'] = pynbody.array.SimArray(t_dyn(snap['r']*.1e3),'yr')
time = header['time']

print(time,snap['t_dyn'].min(),snap['t_dyn'].max())

fig,sp = plt.subplots()
hist,xs, patches = sp.hist(snap['t_dyn'],bins=100,density=True)
sp.axvline(x=time*1.e6)
fig.savefig(f"../figures/t_dyn_{run_id}_{run_name}_{snap_id}.png")



hist,xs,ys = np.histogram2d(snap['r'],snap['vr'],bins=(100,100))
xcent = (xs[1:]+xs[:-1])/2.

def ballistic_vkms_kpc(r_kpc,v0=100.):
    r_pc = r_kpc*1000
    return np.sqrt(2*((v0**2)/2-21443. + 4302268/(250.+r_pc)+4302.3/np.sqrt(r_pc**2+0.01**2)))

fig,sp = plt.subplots()

sp.pcolormesh(xs,ys,np.log10(hist))
sp.set_xlabel(r'$r$ (kpc)')
sp.set_ylabel(r'$v_r$ (km/s)')
sp.plot(xcent,(xcent-3.e-2)*800.,ls='--')

# sp.plot(xcent,ballistic_vkms_kpc(10**xcent))
# print(xcent)
# print(ballistic_vkms_kpc(xcent))
# sp.plot(xcent,ballistic_vkms_kpc(10**xcent,v0=200.))
# sp.plot(xcent,ballistic_vkms_kpc(10**xcent,v0=300.))
# sp.plot(xcent,ballistic_vkms_kpc(10**xcent,v0=400.))
# sp.plot(xcent,ballistic_vkms_kpc(10**xcent,v0=1000.))

fig.savefig("../figures/ballistic.png")