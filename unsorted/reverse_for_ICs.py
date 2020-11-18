import numpy as np
import matplotlib.pyplot as P
from mpl_toolkits.mplot3d import Axes3D
import os
import sys
from multiprocessing import Pool

G_cgs = 6.67259e-8
solarmass_in_g = 1.989e33
pc_in_cm = 3.086e18
yr_in_s = 31556926
km_in_cm = 1.e5
G_internal = G_cgs/pc_in_cm**2*solarmass_in_g/km_in_cm

class grav_pot:
    def __init__(self,m_smbh,c_scale,m_hernquist,a_hernquist):
        self.m_smbh=m_smbh
        self.c_scale=c_scale
        self.m_hernquist=m_hernquist
        self.a_hernquist=a_hernquist
    
    def accel(self,xyz):
        r=np.sqrt(np.sum(xyz**2,axis=1))
        x=r/self.a_hernquist
        m=self.m_smbh*r**2/(r**2+self.c_scale**2) + self.m_hernquist*(x/(1.+x))**2
        return (-G_internal*m/r**3)[:,None] * xyz

def random_sphere(N,rad=1.):
    xyz = np.zeros((N,3))
    phi = 2.*np.pi*np.random.random(N)
    theta = np.arccos(2.*np.random.random(N)-1.)
    xyz[:,0] = rad*np.sin(theta)*np.sin(phi)
    xyz[:,1] = rad*np.sin(theta)*np.cos(phi)
    xyz[:,2] = rad*np.cos(theta)
    return xyz

pot = grav_pot(1.e6,1.e-2,1.e9,250.)

N = 1000
rmin = 0.005
inflow_rate=1. # clouds/year


# random initial positions, moving outwards at escape velocity, then add spin
xyz = random_sphere(N,rmin)
v_esc=np.sqrt(G_internal*pot.m_smbh/rmin*(pc_in_cm/km_in_cm))
# vel = random_sphere(N,v_esc)
r=np.sqrt(np.sum(xyz**2,axis=1))
vel=xyz*v_esc/r[:,None]

#distribute angular momentum randomly
mean_rot_vel = 100.
vcirc = (np.random.random([N])*3.-1.)*mean_rot_vel
R=np.sqrt(xyz[:,0]**2+xyz[:,1]**2)
vel[:,0]+=-vcirc*xyz[:,1]/R
vel[:,1]+=vcirc*xyz[:,0]/R
vcirc_avg = np.sum((-vel[:,0]*xyz[:,1] + vel[:,1]*xyz[:,0])/R)/N
#correct to make sure mean specific angular momentum is constant between runs
dvcirc=mean_rot_vel-vcirc_avg
vel[:,0]+=-dvcirc*xyz[:,1]/R
vel[:,1]+=dvcirc*xyz[:,0]/R

P.close('all')
fig = P.figure()
# ax = fig.add_subplot(111,projection='3d')
ax = fig.add_subplot(111)

os.system("rm pics/frame*.png")
# xyz in pc, vel in km/s, accel in km/s/s, time in years, I did not choose this well
dt = 1.e-2
t=0.
# tmax = N/inflow_rate
tmax = 1.e2
itick=0
iframe=0

minprocs=32
start_indices = np.arange(0,N,N//minprocs)
nprocs=start_indices.size
end_indices = np.zeros([nprocs],dtype=np.int)
end_indices[0:-1] = start_indices[1:]
end_indices[-1]=N

iprocs = np.arange(nprocs)

dump_steps=100

def step_particles(iproc):
    global xyz,vel
    start = start_indices[iproc]
    end = end_indices[iproc]
    for i in range(dump_steps):
        xyz[start:end,:]+=vel[start:end,:]*dt/977813.106
        vel[start:end,:]+=pot.accel(xyz[start:end,:])*dt/3.16887646e-8

pool=Pool(processes=nprocs)
while t<tmax:
#     if iframe%100==0: print(iframe,t,"/",tmax)
    t+=dt*dump_steps
    nslice=N-int(inflow_rate*t)
#         step_particles(0,nslice)
    pool.map(step_particles,iprocs)
#     if itick%100==0:
    print(itick,nslice,"/",N,iframe,t,"/",tmax)
    R=np.sqrt(xyz[:,0]**2+xyz[:,1]**2)
    ax.scatter(R,xyz[:,2])
#         ax.scatter(xyz[:,0],xyz[:,1],xyz[:,2])
    ax.set_xlim(0.,1.e-1)
    ax.set_ylim(-5.e-2,5.e-2)
#         ax.set_zlim(-5.e-2,5.e-2)
    P.savefig("pics/frame%06d.png"%iframe)
    iframe+=1
#     itick+=1
    itick+=dump_steps
#     ax.scatter(a[:,0],a[:,1],a[:,2])
    ax.clear()

cmd = "ffmpeg -y -r 24 -i pics/frame%06d.png -c:v mpeg4 -q:v 1 ../movies/trajectory.mp4"
os.system(cmd)

# ax.scatter(xyz[:,0],xyz[:,1],xyz[:,2])

# np.savetxt("test.dat",xyz)