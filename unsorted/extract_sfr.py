import numpy as np
import gizmo_tools

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

munit = 1.e10 # solar masses
tunit = 0.9778e9 # years

# def moving_mean(x,N=5):
#     return np.convolve(x, np.ones((N,))/N, mode='same')

# class LineSkipper:
#     def __init__(self,f,nskip=1000):
#         self.f = f
#         self.nskip = nskip
#         
#     
#     def __iter__(self):
#         return self
#     
#     def __next__(self):
#         for iline in range(self.nskip):
#             linetxt = self.f.readline()
#         return linetxt
# #         for iline,linetxt in enumerate(self.f):
# #             if iline%nskip==0:
# #                 yield linetxt
# #         raise StopIteration()
#     
#     def next(self):
#         for iline in range(self.nskip):
#             linetxt = self.f.readline()
#         return linetxt
# #         try:
# #             while True:
# #                 yield self.__next__()
# #         except StopIteration:
# #             raise StopIteration()

def extract_sfr(run_id,output_dir,save=True,nskip=10):
#         print("Plotting: ",output_dir)
    gizmoDir = gizmo_tools.getGizmoDir(run_id)
    fullDir = gizmoDir+"/"+run_id+"/"+output_dir
    infilename = fullDir+"/info.txt"
    print("Loading and parsing")
    sf_data_str = np.loadtxt(infilename,delimiter=",",usecols=(1,4),dtype=str)
#         sf_data_str = np.loadtxt(LineSkipper(open(infilename)),delimiter=",",usecols=(1,4),dtype=str)
#         sf_data_str = np.loadtxt(LineSkipper(open(infilename)),delimiter=",",dtype=str)
#         exit()
#         nskip = 1000
#     nskip = 10
    timeskip = 0.01*1.e6/tunit
#     timeskip = 0.001*1.e6/tunit
    
    print("input size ",sf_data_str.shape)
    
    sf_data_str = sf_data_str[::nskip,:]

    times = np.loadtxt(sf_data_str[:,0],usecols=[1])
    dt = np.gradient(times)
    good_dt = dt>0.
    times = times[good_dt]

    final_time = times[-1]
    first_time = times[0]
    sample_times = np.arange(first_time,final_time,timeskip)
    sample_locs = np.searchsorted(times,sample_times)

#         print(sample_times)
#         
#         
#         print(sample_locs)

    sf_data_str = sf_data_str[good_dt,:][sample_locs,:]

    print("resampled size ",sf_data_str.shape)

    
    if sf_data_str.shape[0]>3:
    
        times = np.loadtxt(sf_data_str[:,0],usecols=[1])
        sf_cumulative = np.loadtxt(sf_data_str[:,1],usecols=[1])
    
        times *= tunit
        sf_cumulative *= munit
    
#         times = times[good_dt]
#         sf_cumulative = sf_cumulative[good_dt]
        # fix discontinuities from restarts losing count of cumulative sf
        dsf = np.diff(sf_cumulative)
        if np.any(dsf<0.):
            jump_args = np.argwhere(dsf<0.)
            print(jump_args)
            for jump_arg_array in jump_args:
                jump_arg = jump_arg_array[0]
                jump_size = sf_cumulative[jump_arg]
                sf_cumulative[jump_arg+1:]+=jump_size
#                 print(jump_arg)
#                 print(sf_cumulative[jump_arg-1:jump_arg+2])
#                 print(dsf[jump_arg-1:jump_arg+2])
        sfr = np.gradient(sf_cumulative)/np.gradient(times)
        
        times /= 1e6 # Myr for plot
    else:
        times = [0,0]
        sf_cumulative = [0,0]
        sfr = [0,0]
        print("No times")
    print("Dumping")
#         P.plot(times,sfr,label=runs[irun])
    
    if save:
        output_data = np.vstack((times,sf_cumulative,sfr)).T
        np.savetxt("data/sf"+run_id+output_dir+".dat",output_data)
    return times,sf_cumulative,sfr
    
#     P.xlabel("$t$ (Myr)")
#     P.ylabel("SFR (M$_\odot$/yr)")
#     P.legend(loc='best')
#     P.yscale('log')
#     P.savefig("../figures/sfr_q2.png",dpi=200)
#     P.close()
if __name__ == '__main__':
#     run_id = "2014"
#     runs = ["q2edd05redo","q2edd10_aniso1redo","q2edd10_aniso3redo","q2edd10redo","q2edd20redo","q2redo"]

#     run_ids     = ["2014","2015","2015","2014","2014","2014"]
#     runs = ["q2edd05redo","q2edd10_aniso1fixed","q2edd10_aniso3fixed","q2edd10redo","q2edd20redo","q2redo"]
#     run_ids     = ["2015","2015"]
#     runs = ["q2edd10_aniso1fixed","q2edd10_aniso3fixed"]

#     run_ids     = ["1055","1055"]
#     runs = ["q2_SN","q2_SN_slow"]
#     run_ids     = ["2020"]*2
#     runs = ["restest0m01_small","restest0m02_small"]
#     run_ids     = ["2020"]*3
#     runs = ["restest0m005_tiny","restest0m01_tiny","restest0m02_tiny"]
#     runs = ["restest0m02","restest0m04"]
#      run_ids     = ["1055","1055"]
#     runs = ["q2_SN","q2_SN_slow"]

#     run_ids     = ["2014","2015","2015","2014","2014","2014","2015","2015","1055","1055"]
#     runs = ["q2edd05redo","q2edd10_aniso1fixed","q2edd10_aniso3fixed","q2edd10redo","q2edd20redo","q2redo","q2edd10_aniso1fixed","q2edd10_aniso3fixed","q2_SN","q2_SN_slow"]
#     run_ids = ["2018"]
#     runs = ["treetest"]
#     run_ids = ["2019"]*5
#     runs = ["restest0m"+x for x in ["02","04","06","08","1"]]
#     run_ids = ["2020"]*5
#     runs = ["aniso_"+x for x in ["0","1","2","3","4"]]


#     runs = ["SF_test_high_rho_floor","SF_test_truelove","SF_test_high_rho","SN_test_m2_slowerlong"]
#     runs = ["SF_test_high_rho_floor_1000","SF_test_high_rho_floor_100","SF_test_high_rho_floor_30","SF_test_high_rho_floor","no_SNe_long"]
#     run_ids = ["2030"]*4+["3030"]
#     runs = ["SF_test_high_rho_floor_100","SF_test_high_rho_floor_30","SF_test_high_rho_floor_100_thinner_Q1","SF_test_high_rho_floor_30_thinner_Q1"]
    runs = ["SF_test_high_rho_floor_30_thinner_Q1_5_cut","SF_test_high_rho_floor_30_thinner_Q1_cut","SF_test_high_rho_floor_30_thinner_Q2_cut","SN_test_high_rho_floor_30_thinner_Q1_5_cut","SN_test_high_rho_floor_30_thinner_Q1_cut","SN_test_high_rho_floor_30_thinner_Q2_cut"]
    run_ids = ["2030"]*6
    save = False
    nskip = 10
    sfr_filename = "../figures/sfr_floor.png"
    cum_sf_filename = "../figures/cum_sf_floor.png"
    

    print(run_ids)
    print(runs)
    sfr_data = []
    for run_id,output_dir in zip(run_ids,runs):
        sfr_data+=[extract_sfr(run_id,output_dir,save=save,nskip=nskip)]
    
    if sfr_filename:
        fig = plt.figure()
        plt.xlabel("$t$ (Myr)")
        plt.ylabel("SFR (M$_\odot$/yr)")
        for run_name,(times,sf_cumulative,sfr) in zip(runs,sfr_data):
            plt.plot(times,sfr,label=run_name)
        plt.yscale('log')
        plt.legend(loc='best')
        plt.grid(which='both',linewidth=.5,color='grey')
        plt.savefig(sfr_filename,dpi=200)
        plt.close()
    
    if cum_sf_filename:
        fig = plt.figure()
        plt.xlabel("$t$ (Myr)")
        plt.ylabel("SF (M$_\odot$)")
        for run_name,(times,sf_cumulative,sfr) in zip(runs,sfr_data):
            plt.plot(times,sf_cumulative,label=run_name)
        plt.yscale('log')
        plt.legend(loc='best')
        plt.grid(which='both',linewidth=.5,color='grey')
        plt.savefig(cum_sf_filename,dpi=200)
        plt.close()

