
########!!!!!! UNFINISHED
# DOES NOTHING YET!

import numpy as np
import pylab as P

def vec_kernel(r,h):
    x = r/h
    k = np.zeros(x.size)
    k[x<.5] = 1.+6.*(x[x<.5]-1.)*x[x<.5]**2
    k[(x<1.)&(x>=.5)] = 2.*(1.-x[(x<1.)&(x>=.5)])**3
    # normalise
    k*=8./np.pi/h**3
    return k

def vec_kernel_inv(r,h):
    return np.where(r>0.,h/r,h/(1.e-9))/h**3

# def flatten_kernel(r,h):
#     dz = 1.e-2
#     rout = 0.
#     for zz in np.arange(-100.*h,100.*h,dz):
#         rout = rout + kernel(np.sqrt(r**2+zz**2),h)*dz
#     return dz
# print("Integrating")
# dr_flatkern = 1.e-3
# dz_flatkern = 1.e-3
# rflat = np.arange(0,10.,dr_flatkern)
# flat_kern = np.zeros(rflat.size)
# for zz in np.arange(-10.,10.,dz_flatkern):
#     flat_kern = flat_kern + vec_kernel(np.sqrt(rflat**2+zz**2),1.)*dz_flatkern
# 
# print("Integration done!")
# 
# def vec_flat_kernel(r,h):
#     return np.where(r/h<10.,flat_kern[np.floor(r/h/dr_flatkern).astype(int)]/h**2,0.)

# 
# def vec_flatten_kernel(r,h):
#     dz = 1.e-4
#     rout = np.zeros(r.size)
#     for zz in np.arange(-10.,10.,dz):
#         rout = rout + vec_kernel(np.sqrt(r**2+zz**2),h)*dz
#     return rout




slice_z = .0

r_p = np.random.random((N,3))
#h = np.random.random(N)*.3+.2
h = .2
r_p[:,2] = slice_z
#r_p[:,1] = .0
#r_p[:,0] = np.linspace(0.,1.,N)
#r_p[:,0] = r_p[:,0] * .95
#p_p = np.linspace(0.,1.,N)
#p_p = np.random.random(N)
#p_p[::2] = .5
#p_p[1::2] = 1.
#p_p[:] = 2.
#p_p[N/2:] = 1.

p_p = np.exp(-2.*((r_p[:,0]-.5)**2+(r_p[:,1]-.5)**2))*(np.sin(r_p[:,0]*np.pi*4.)+1)

#p_p = r_p[:,0]

# r_p[0,0] = .1
# r_p[0,1] = .1
# r_p[1,0] = .9
# r_p[1,1] = .1
# r_p[2,0] = .9
# r_p[2,1] = .9

#slice_L = 4096
slice_L = 128

#E_mat = np.zeros((3,3))
#2d
E_mat = np.zeros((2,2))

r_taylor = np.array([.5,.5,slice_z])

# t_rads = np.zeros(N)
# t_fits = np.zeros(N)


# for ip in range(N):
#     r0 = r_p[ip]
#     pp0 = p_p[ip]
#     h0 = h#[ip]
# 
#     dr3_p = r_p-r0
#     dr_p = np.linalg.norm(dr3_p,axis=1)
#     w_p = vec_kernel(dr_p,h)
# 
#     for ix in range(2):
#         for iy in range(2):
#             E_mat[ix,iy] = np.sum((r_p[:,ix]-r0[ix])*(r_p[:,iy]-r0[iy])*w_p[:])
#     
#     if ( np.sum(E_mat!=0.) ):
#         B_mat = np.linalg.inv(E_mat)
# 
#         gradvec = np.zeros(2)
# 
#         for ix in range(2):
#             for iy in range(2):
#                 gradvec[ix]+=np.sum(B_mat[ix,iy]*(p_p[:]-pp0)*(r_p[:,iy]-r0[iy])*w_p[:])
#     
#         #taylor fit
#         tf = np.sum((r_taylor-r0)[:2]*gradvec)+pp0
#         print(ip,tf,np.linalg.norm(r0-r_taylor))
#         t_rads[ip] = np.linalg.norm(r0-r_taylor)
#         t_fits[ip] = tf
#         #print(ip,r0[0:2],gradvec)
#     else:
#         print(ip,r0[0:2])
#         #print(ip,np.linalg.norm(r0-r_taylor))



# slice_line = np.zeros(slice_L)
# slice_xs = np.linspace(0.,1.,slice_L)
# 
# for ix in range(slice_L):
#     rg = np.array((slice_xs[ix],0.,0.))
#     dr_p = np.linalg.norm(r_p-rg,axis=1)
#     w_p = vec_kernel(dr_p,h)
#     w_t = np.sum(w_p)
#     slice_line[ix] = np.sum(w_p*p_p)/w_t
# 
# P.plot(slice_xs,slice_line)
# 
# print(np.sum(slice_line)/slice_line.size)
# 
# rflat = np.sqrt(r_p[:,2]**2+r_p[:,1]**2)
# lit = vec_flat_kernel(rflat,h)*p_p
# lit = np.sum(lit)/np.sum(vec_flat_kernel(rflat,h))
# print(lit)

dr3_p = r_p-r_taylor
dr_p = np.linalg.norm(dr3_p,axis=1)
w_p = vec_kernel(dr_p,h)

def pvecj(idex):
    if idex==0:
        return 1.
    else:
        return dr3_p[:,idex-1]

def pvecval(idex):
    if idex==0:
        return 1.
    else:
        #return r_taylor[idex-1]
        return 0.

# fit_mat = np.zeros((3,3))
# for ix in range(3):
#     for iy in range(3):
#         fit_mat[ix,iy] = np.sum(pvecj(ix)*pvecj(iy)*w_p)
# 
# inv_mat = np.linalg.inv(fit_mat)
# 
# phis = 0.
# for ix in range(3):
#     for iy in range(3):
#         phis+=inv_mat[ix,iy]*pvecj(ix)*pvecval(iy)*w_p
# 
# r_estimate = np.sum(phis*p_p)
# print(r_estimate)

###

slice_grid = np.zeros((slice_L,slice_L))
edges = np.linspace(0.,1.,slice_L)

for gx in range(slice_L):
#for gx in [slice_L/2]:
    print(gx)
    for gy in range(slice_L):
    #for gy in [slice_L/2]:
        r_taylor = np.array(((gx+.5)*1./slice_L,(gy+.5)*1./slice_L,slice_z))
        dr3_p = r_p-r_taylor
        dr_p = np.linalg.norm(dr3_p,axis=1)
        w_p = vec_kernel(dr_p,h)
        fit_mat = np.zeros((3,3))
        for ix in range(3):
            for iy in range(3):
                fit_mat[ix,iy] = np.sum(pvecj(ix)*pvecj(iy)*w_p)
                
        if ( np.linalg.matrix_rank(fit_mat)!=3 ):
            slice_grid[gx,gy] = -1
        else:
            inv_mat = np.linalg.inv(fit_mat)

            phis = 0.
            for ix in range(3):
                for iy in range(3):
                    phis+=inv_mat[ix,iy]*pvecj(ix)*pvecval(iy)*w_p

            r_estimate = np.sum(phis*p_p)
        
            slice_grid[gx,gy] = r_estimate
#         rg = np.array(((ix+.5)*1./slice_L,(iy+.5)*1./slice_L,slice_z))
#         dr_p = np.linalg.norm(r_p-rg,axis=1)
#         w_p = vec_kernel_inv(dr_p,h)
#         w_t = np.sum(w_p)
#         slice_grid[ix,iy] = np.sum(w_p*p_p)/w_t
        #slice_grid[ix,iy] = p_p[np.argmax(w_p)]

#print(slice_grid[:,slice_L/2])

#P.plot(edges,slice_grid[:,slice_L/2])

P.pcolormesh(edges,edges,slice_grid.T,vmin=-1.,vmax=2) #
P.colorbar()
P.scatter(r_p[:,0],r_p[:,1],c=p_p,vmin=-1.,vmax=2)

P.savefig("sliceplay.png",dpi=200)
P.close('all')
