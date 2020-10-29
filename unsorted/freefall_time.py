import numpy as np
import functools

m_smbh = 1.e6
smbh_soft = 0.01
m_hernquist = 1.e9
a_hernquist = 250.


def accel_smbh_hernquist(m_smbh,smbh_soft,m_hernquist,a_hernquist,r):
"""
accel_smbh_hernquist

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

def free_fall(r_in,r_cutoff=1.,dt_crit=0.01):
    v = 0
    r = r_in
    t=0
    while r>r_cutoff:
        a = accel(r)
        dt = dt_crit*np.sqrt(np.abs(r/a))
        r+=v*dt
        v+=a*dt
        t+=dt
    return t,v

def t_dyn(r):
    return np.sqrt(r/-accel(r))

r = np.arange(0.1,1.e4,1.)
outp = np.array([r,t_dyn(r)]).T

np.savetxt("data/tdyn.dat",outp)


# outp = []
# 
# # for r in np.arange(1.,100.,1.):
# for r in np.arange(1.,1.e5,1.e4):
#     t_ff,v_esc = free_fall(r)
#     outp+=[[r,t_ff/1.e6,v_esc*977813.106]]
# 
# np.savetxt("data/freefall.dat",outp)