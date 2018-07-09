import sys
import gizmo_tools
import numpy as np
import h5py
from scipy import interpolate
from scipy import signal
from multiprocessing import Pool
import itertools as it
import sys

# bins = np.arange(0.,91.,1.)
bins = np.arange(0.,91.,.25)

# windcut = 10. # minimum outflow speed, km/s

G_kms_pc_msun = 0.0043022682

rcut = 80. # only count particles within 80 pc

def v_esc(r,m_bh,m_hern,a_hern):
    return np.sqrt(2*G_kms_pc_msun*(m_bh/r + m_hern/(a_hern+r)))


def smooth(x,y,width):
    sample_indices = np.unique(np.searchsorted(x,np.arange(x[0],x[-1],width)))
    x_samples = x[sample_indices]
    y_samples = y[sample_indices]
    print(x_samples)
    print(y_samples)
    f_interp = interpolate.interp1d(x_samples,y_samples,kind='cubic',fill_value="extrapolate")
    y_out = f_interp(x)
    return y_out


#     x_samples = x[::width]
#     y_samples = y[::width]
#     f_interp = interpolate.interp1d(x_samples,y_samples,kind='cubic',fill_value="extrapolate")
#     y_out = f_interp(x)
#     return y_out

#     lowx = x-width
#     highx = x+width
#     low_indices = np.searchsorted(x,lowx)
#     high_indices = np.searchsorted(x,highx)
#     smoothed = [np.mean(x[low_index:high_index]) for low_index,high_index in zip(low_indices,high_indices)]
#     
#     print(low_indices)
#     print(high_indices)
#     print(smoothed)
#     return smoothed

def load_gizmo_get_outflow(run_id,output_dir,snap_id):
    gizmoDir = gizmo_tools.getGizmoDir()
    fullDir = gizmoDir+"/"+run_id+"/"+output_dir
    snap_str = "%03d" % snap_id
    f = h5py.File(fullDir+"/snapshot_"+snap_str+".hdf5","r")

    header = f["/Header"]
    time = header.attrs.get("Time")
    time*= 0.9778e9 # to yr
    time/=1.e6 # to Myr

    xyz_p_in = np.array(f["/PartType0/Coordinates"])
    vel_p_in = np.array(f["/PartType0/Velocities"])
    mass_p_in = np.array(f["/PartType0/Masses"])
    id_p = np.array(f["/PartType0/ParticleIDs"])-1

#     id_p = np.zeros_like(id_p_in)
    xyz_p = np.zeros_like(xyz_p_in)
    vel_p = np.zeros_like(vel_p_in)
    mass_p = np.zeros_like(mass_p_in)
    
    xyz_p[id_p,:] = xyz_p_in
    vel_p[id_p,:] = vel_p_in
    mass_p[id_p] = mass_p_in
#     id_p[id_p_in] = id_p_in

    mass_p*=1.e10 # to Msun
    xyz_p*=1.e3 # to pc
    
    rad2d_p = np.sqrt(xyz_p[:,0]**2+xyz_p[:,1]**2)
    rad3d_p = np.sqrt(xyz_p[:,0]**2+xyz_p[:,1]**2+xyz_p[:,2]**2)
    z_p = xyz_p[:,2]

    vmag = np.sqrt(np.sum(vel_p**2,axis=1))
    
#     vrad = np.sum(xyz_p*vel_p,axis=1)/rad3d_p
#     outflowing = (vrad>windcut)
    outflowing = (v_esc(rad3d_p,1.e6,1.e9,250.)<vmag) & (rad3d_p<rcut)
    
    vtot = np.sum(vmag)
    
    return time,z_p,rad2d_p,outflowing,mass_p,vtot
    

def calc_hist(run_id,output_dir,snap_id,dumpall,bins,step):
    time,z_p,rad2d_p,outflowing,mass_p,vtot = load_gizmo_get_outflow(run_id,output_dir,snap_id)

    z_p = z_p[outflowing]
    rad2d_p = rad2d_p[outflowing]
    mass_p = mass_p[outflowing]
    
#    theta_p = np.abs(np.arctan2(z_p,rad2d_p)*180./np.pi)
    theta_p = np.abs(np.arctan2(z_p,rad2d_p)*180./np.pi)
    
    if not dumpall:
        m_out = np.sum(mass_p)
        
        if m_out>0.:
            wind_theta = np.median(theta_p)
        
            last_time,z_p_old,rad2d_p_old,last_outflowing,dummy,vtot_old = load_gizmo_get_outflow(run_id,output_dir,snap_id-step)
        
            newly_outflowing = outflowing & ~last_outflowing
            nolonger_outflowing = last_outflowing & ~outflowing
            particle_mass = mass_p[0]
        
            dt = (time-last_time)*1.e6
#             
#             out_data1 = np.array([z_p_old[newly_outflowing],rad2d_p_old[newly_outflowing]]).T
#             np.savetxt("data/newflow%03d.dat"%snap_id,out_data1)
#             out_data2 = np.array([z_p_old[nolonger_outflowing],rad2d_p_old[nolonger_outflowing]]).T
#             np.savetxt("data/oldflow%03d.dat"%snap_id,out_data2)
            
        
            outflow_rate = particle_mass*np.sum(newly_outflowing)/dt
            inflow_rate = particle_mass*np.sum(nolonger_outflowing)/dt
            dv = (vtot-vtot_old)/dt
            
        else:
            wind_theta = -1.
            outflow_rate = -1.
            inflow_rate = -1.
            dv = -1.

        
        
#         m_out = np.sum(outflowing)
        return snap_id,time,wind_theta,m_out,outflow_rate,inflow_rate,dv

    
#     indices = np.argsort(theta_p)
    theta_out = np.sort(theta_p)
    
#     bin_size = 1000
#     theta_edges = theta_out[::bin_size]
#     theta_histogram = bin_size/(theta_edges[1:]-theta_edges[:-1])
#     theta_cents = (theta_edges[1:]+theta_edges[:-1])/2.
#     smoothed = theta_histogram
#     grad = np.gradient(theta_histogram)/np.gradient(theta_cents)
#     gradgrad = np.gradient(grad)/np.gradient(theta_cents)
    
    
#     theta_cumsum = np.arange(theta_out.size)


    theta_histogram, theta_edges = np.histogram(theta_p,weights=mass_p,bins=bins)
    
    slice = (theta_histogram>0)
    theta_histogram = theta_histogram[slice]

#     log_hist = np.log10(theta_histogram)
#     smoothed = theta_histogram
    theta_cents = (theta_edges[1:]+theta_edges[:-1])/2.
    theta_cents = theta_cents[slice]
    smoothed = theta_histogram
    grad = np.gradient(smoothed)/np.gradient(theta_cents)
    gradgrad = np.gradient(grad)/np.gradient(theta_cents)
    return theta_cents,theta_histogram,smoothed,grad,gradgrad
    
#     knots = np.arange(theta_cents[1],theta_cents[-1],2.)
    
#     log_hist = np.log10(theta_histogram)
#     theta_fit = interpolate.splrep(theta_cents,log_hist,t=knots,k=5)
#     print(theta_fit)
#     print(theta_fit[1].size)
#     smoothed = interpolate.splev(theta_cents,theta_fit)
#     grad = interpolate.splev(theta_cents,theta_fit,der=1)
#     gradgrad = interpolate.splev(theta_cents,theta_fit,der=2)
    
    #smoothed = log_hist


#     fixed_peak_magnitude = np.max(theta_histogram)
#     fit_func = fixed_peak_gaussian
# 
#     popt, pcov = optimize.curve_fit(fit_func, theta_centres, theta_histogram)
#     print(popt)
#     print(pcov)
#   ,fit_func(theta_edges[:-1],*popt)

##    theta_fit = interpolate.splrep(theta_out,theta_cumsum,t=theta_out[1:-1:10000])
#     theta_fit = interpolate.splrep(theta_out,theta_cumsum,t=np.arange(1.e-1,10.,1.))
#     print(theta_fit)
#     print(theta_fit[1].size)
#     smoothed = interpolate.splev(theta_out,theta_fit)
#     grad = interpolate.splev(theta_out,theta_fit,der=1)
#     gradgrad = interpolate.splev(theta_out,theta_fit,der=2)

#     smoothed = signal.savgol_filter(theta_cumsum,polyorder=3,window_length=99)
#     smoothed = signal.savgol_filter(theta_out,polyorder=3,window_length=99)
#     smoothed = theta_cumsum
#     order = 6
#     grad = (theta_cumsum[order:]+theta_cumsum[:-order])/order
#     first_value = grad[0]
#     last_value = grad[-1]
#     grad = np.insert(grad,np.full((order//2),first_value),0)
#     grad = np.append(grad,np.full((order//2),last_value))
#     gradgrad = np.gradient(grad)
#     smoothed = smooth(theta_out,theta_cumsum,10000)
#     smoothed = smooth(theta_out,theta_cumsum,.01)

#     grad = np.gradient(smoothed)/np.gradient(theta_out)
# #     
# #     grad = signal.savgol_filter(grad,polyorder=3,window_length=99)
#     
#     smoothgrad = smooth(theta_out,grad,.01)
#     #gradgrad = smoothgrad
#     gradgrad = np.gradient(smoothgrad)/np.gradient(theta_out)

#     grad = signal.savgol_filter()

#     grad = 1./np.gradient(theta_out)
#     gradgrad = np.gradient(grad)
#     wind_theta = theta_centres[np.argmin(gradgrad)]
    
#     if dumpall:
# #         return theta_edges,theta_histogram,wind_theta,grad,gradgrad
# #         return theta_out,theta_cumsum,smoothed,grad,gradgrad
#         return theta_cents,theta_histogram,smoothed,grad,gradgrad
#     else:
#         return time,wind_theta

def dump_full_evolution(run_ids,output_dirs):
    for run_id,output_dir in zip(run_ids,output_dirs):
        print("Dumping full evolution",run_id,output_dir)
        
#             print("Loading and parsing info.txt")
#             gizmoDir = gizmo_tools.getGizmoDir()
#             fullDir = gizmoDir+"/"+run_id+"/"+output_dir
#             info_filename = fullDir+"/info.txt"
#             info_data_str = np.loadtxt(info_filename,delimiter=",",usecols=(0,1),dtype=str)
#             info_syncpoints = np.loadtxt(info_data_str[:,0],usecols=[1])
#             info_times = np.loadtxt(info_data_str[:,1],usecols=[1])
#             info_times*=0.9778e9 # to yr
#             info_times/=1.e6 # to Myr
#             print(info_syncpoints)
#             print(info_times)
        
        print("Reading and analysing angles")

        snapf = gizmo_tools.lastConsecutiveSnapshot(run_id,output_dir)
#             step = 10
        step = 1
        snapi = step
#         snapf = 999
#             snap_strs = ["%03d" % i for i in range(0,snapf,step)]
        snap_ids = range(snapi,snapf,step)
        bins = np.arange(0.,91.,1.)
        with Pool(processes=80) as pool:
#             with Pool(processes=8) as pool:
            output_data = pool.starmap(calc_hist,zip(it.repeat(run_id),it.repeat(output_dir),snap_ids,it.repeat(False),it.repeat(bins),it.repeat(step)))
#         print(output_data)
        output_data = np.array(output_data)
        
#             print(output_data[:,1])
        
#             output_time_ids = np.searchsorted(info_times,output_data[:,1])
#             output_syncpoints = info_syncpoints[output_time_ids]
#             output_data = np.vstack([output_data.T,output_syncpoints]).T
        
#             mdot = np.gradient(output_data[:,3])/np.gradient(output_data[:,1])
#             mdot/=1.e6 # Msun/yr
#             mdotdot = np.gradient(mdot)/np.gradient(output_data[:,1])#/mdot#/np.gradient(output_data[:,0])
#             output_data = np.vstack([output_data.T,mdot,mdotdot]).T
        np.savetxt("data/windangle_evolution"+run_id+output_dir+".dat",output_data)
    
#         cut = (output_data[:,0]>324) & (mdotdot>33.5)
#         t_jumps = output_data[cut,0]
#         np.savetxt("data/t_jumps.dat",np.array([range(t_jumps.size),t_jumps]).T)

 
if __name__ == '__main__':
    if len(sys.argv)>=4:
        print("Dumping one file")
        run_id = sys.argv[1]
        output_dir = sys.argv[2]
        snap_str = sys.argv[3]
        bins = np.arange(0.,91.,.25)
#         theta_edges,theta_histogram,wind_theta,grad,gradgrad = calc_hist(run_id,output_dir,snap_str,True,bins)
#         output_data = np.array([theta_edges[:-1],theta_histogram,grad,gradgrad]).T
        theta,cumsum,smoothed,grad,gradgrad = calc_hist(run_id,output_dir,snap_str,True,bins)
        output_data = np.array([theta,cumsum,smoothed,grad,gradgrad]).T
        np.savetxt("data/quickangle"+run_id+output_dir+snap_str+".dat",output_data)
    else:
        if len(sys.argv)>1:
            run_ids = [sys.argv[1]]
            output_dirs = [sys.argv[2]]
        else:
#             run_id = "2014"
#             output_dirs = ["q2edd05redo","q2edd10_aniso1redo","q2edd10_aniso3redo","q2edd10redo","q2edd20redo","q2redo"]
#             run_ids     = ["2014","2015","2015","2014","2014","2014"]
#             output_dirs = ["q2edd05redo","q2edd10_aniso1fixed","q2edd10_aniso3fixed","q2edd10redo","q2edd20redo","q2redo"]
#             run_ids     = ["1055","1055"]
#             output_dirs = ["q2_SN","q2_SN_slow"]
#             run_ids     = ["2015","2015"]
#             output_dirs = ["q2edd10_aniso1fixed","q2edd10_aniso3fixed"]
#             run_ids     = ["2015"]
#             output_dirs = ["q2edd10_aniso1fixed"]
#             run_ids     = ["2020"]*2
#             output_dirs = ["restest0m01_small","restest0m02_small"]
#             run_ids     = ["2020"]*3
#             output_dirs = ["restest0m005_tiny","restest0m01_tiny","restest0m02_tiny"]
#             output_dirs = ["restest0m02","restest0m04"]
#             run_ids = ["2018"]
#             output_dirs = ["treetest"]
#             run_ids = ["2019"]*5
#             output_dirs = ["restest0m"+x for x in ["02","04","06","08","1"]]
            run_ids = ["2020"]*5
#             output_dirs = ["aniso_"+x for x in ["0","1","2","3","4"]]
            output_dirs = ["aniso_4_biggerrad","aniso_4_bigrad"]
            print(run_ids)
            print(output_dirs)
            dump_full_evolution(run_ids,output_dirs)
    

