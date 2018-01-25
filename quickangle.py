import sys
import gizmo_tools
import numpy as np
import h5py
from scipy import optimize
from multiprocessing import Pool
import itertools as it

#bins = np.arange(0.,91.,1.)
bins = np.arange(0.,91.,1.)

# def three_spaced_gaussian(x,a0,a1,s0,s1,x1):
#     return a0*np.exp(-(x[x>=0]/s0)**2)+a1*np.exp(-((x[x>=0]-x1)/s1)**2)+

# def three_spaced_gaussian(x,a0,a1,s0,s1,x1):
#     return a0*np.exp(-(x/s0)**2)+a1*np.exp(-((x-x1)/s1)**2)+a1*np.exp(-((x+x1)/s1)**2)
# 
# 
# def one_gaussian(x,a,s):
#     return a*np.exp(-(x/s)**2)
# 
# def fixed_peak_gaussian(x,s):
#     return fixed_peak_magnitude*np.exp(-(x/s)**2)
# 
# 
# def linear(x,m,c):
#     return m*x+c

def calc_hist(run_id,output_dir,snap_str,dumpall):
    gizmoDir = gizmo_tools.getGizmoDir()
    fullDir = gizmoDir+"/"+run_id+"/"+output_dir
    f = h5py.File(fullDir+"/snapshot_"+snap_str+".hdf5","r")

    header = f["/Header"]
    time = header.attrs.get("Time")
    time*= 0.9778e9 # to yr
    time/=1.e6 # to Myr

    xyz_p = np.array(f["/PartType0/Coordinates"])
    mass_p = np.array(f["/PartType0/Masses"])

    mass_p*=1.e10 # to Msun
    
    rad2d_p = np.sqrt(xyz_p[:,0]**2+xyz_p[:,1]**2)
    z_p = xyz_p[:,2]
    
    theta_p = np.abs(np.arctan2(z_p,rad2d_p)*180./np.pi)

        
    theta_histogram, theta_edges = np.histogram(theta_p,weights=mass_p,bins=bins)
    
    theta_centres = (theta_edges[1:]+theta_edges[:-1])/2.
    
#     fixed_peak_magnitude = np.max(theta_histogram)
#     fit_func = fixed_peak_gaussian
# 
#     popt, pcov = optimize.curve_fit(fit_func, theta_centres, theta_histogram)
#     print(popt)
#     print(pcov)
#   ,fit_func(theta_edges[:-1],*popt)

    grad = np.gradient(theta_histogram)
    gradgrad = np.gradient(grad)
    wind_theta = theta_centres[np.argmin(gradgrad)]
    
    if dumpall:
        return theta_edges,theta_histogram,wind_theta,grad,gradgrad
    else:
        return time,wind_theta
    
if __name__ == '__main__':
    run_id = sys.argv[1]
    output_dir = sys.argv[2]
    
    if len(sys.argv)>=4:
        print("Dumping one file")
        snap_str = sys.argv[3]
        theta_edges,theta_histogram,wind_theta,grad,gradgrad = calc_hist(run_id,output_dir,snap_str,True)
        output_data = np.array([theta_edges[:-1],theta_histogram,grad,gradgrad]).T
        np.savetxt("data/quickangle"+run_id+output_dir+snap_str+".dat",output_data)
    else:
        snapi = 0
        snapf = gizmo_tools.lastConsecutiveSnapshot(run_id,output_dir)
        snap_strs = ["%03d" % i for i in range(snapf+1)]
        with Pool(processes=80) as pool:
            output_data = pool.starmap(calc_hist,zip(it.repeat(run_id),it.repeat(output_dir),snap_strs,it.repeat(False)))
        np.savetxt("data/windangle_evolution"+run_id+output_dir+".dat",output_data)


