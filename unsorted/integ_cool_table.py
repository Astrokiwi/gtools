print("Importing")

# Import libraries to do our magic
import numpy as np
from scipy import integrate

print("Running")

cooltab = np.loadtxt("data/cooltab_simple.dat")

temperatures = cooltab[:,0]

lambdas = -cooltab[:,1]
#currently in erg cm^6 /s. Convert to Gadget-ish units
year = 31556926 # seconds
tunit = .9778e9 * year
#lambdas*=(1.e9*year)/(1.e10) # units of 1.e10 erg, Gyr



molecular_mass = 4./(1.+3.*.76)
proton_mass_cgs = 1.6726e-24
munit = 1.9891e43 # g
eunit = 1.9891e53 # erg

inergs = temperatures/98.49184 # convert to Gadget/GIZMO internal units
#inergs = inergs*1.e10 # convert to cgs (erg/g)

lambdas = lambdas * tunit
lambdas = lambdas /proton_mass_cgs
lambdas = lambdas/1e10
invlambdas = 1./lambdas

lambint = integrate.cumtrapz(invlambdas,inergs,initial=0.)

dout = np.array([inergs,lambint])

# units of inergs are 1e10 erg/g
#units of lambint are 1.6726e-34 * cgs

np.savetxt("data/cooltab_int.dat",dout.T)

def du(dt, n, T0):
    u0 = T0/98.49184 # convert to internal energy
    print(u0)
    lamdbint0 = np.interp(u0,inergs,lambint,left=0.,right=-1)
    if ( lamdbint0==-1 ):
        return
    print(lamdbint0)
    #lamdbint0+= n**2*dt
    lamdbint0+= n*dt/(molecular_mass) # *proton_mass_cgs, take out for nicer units
    print(lamdbint0)
    u_out = np.interp(-lamdbint0,-lambint,inergs)
    print(u_out)
    T_out = u_out*98.49184 # convert to K
    print(T_out)
    
    
    #default method
    lam = np.interp(T0,temperatures,lambdas)
    du = lam*n*dt/(molecular_mass)
    #du/=1.e10 # convert to internal units
    u_out2 = u0+du
    T_out2 = u_out2*98.49184
    print(max(T_out2,10.))
    
    #full integration
    nintsteps = 10000
    littledt = dt/nintsteps
    u_full = u0
    u_out3 = u_full
    for iint in range(nintsteps):
        lam = np.interp(u_out3,inergs,lambdas)
        u_out3 += lam * n * littledt / (molecular_mass)
    T_out3 = u_out3*98.49184
    print(max(T_out3,10.))
    return T_out


du(1.e-7, 141242, 9.99768e+07)
#du(4.e-7, 1.e5, 1.e8)