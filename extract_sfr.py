import numpy as np
import gizmo_tools

# import matplotlib as mpl
# mpl.use('Agg')
# import matplotlib.pyplot as P

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

def extract_sfr(run_ids,runs):
    gizmoDir = gizmo_tools.getGizmoDir()
    for irun,(run_id,output_dir) in enumerate(zip(run_ids,runs)):
        print("Plotting: ",output_dir)
        gizmoDir = gizmo_tools.getGizmoDir(run_id)
        fullDir = gizmoDir+"/"+run_id+"/"+output_dir
        infilename = fullDir+"/info.txt"
        print("Loading and parsing")
        sf_data_str = np.loadtxt(infilename,delimiter=",",usecols=(1,4),dtype=str)
#         sf_data_str = np.loadtxt(LineSkipper(open(infilename)),delimiter=",",usecols=(1,4),dtype=str)
#         sf_data_str = np.loadtxt(LineSkipper(open(infilename)),delimiter=",",dtype=str)
#         exit()
#         nskip = 1000
        nskip = 10
#         timeskip = 0.01*1.e6/tunit
        timeskip = 0.001*1.e6/tunit
        
        print(sf_data_str.shape)
        
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
        
        if sf_data_str.shape[0]>3:
        
            times = np.loadtxt(sf_data_str[:,0],usecols=[1])
            sf_cumulative = np.loadtxt(sf_data_str[:,1],usecols=[1])
        
            times *= tunit
            sf_cumulative *= munit
        
    #         times = times[good_dt]
    #         sf_cumulative = sf_cumulative[good_dt]
            sfr = np.gradient(sf_cumulative)/np.gradient(times)
            
            times /= 1e6 # Myr for plot
        else:
            times = [0,0]
            sf_cumulative = [0,0]
            sfr = [0,0]
            print("No times")
        print("Dumping")
#         P.plot(times,sfr,label=runs[irun])
        
        output_data = np.vstack((times,sf_cumulative,sfr)).T
        np.savetxt("data/sf"+run_id+output_dir+".dat",output_data)
    
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
    run_ids = ["2020"]*5
    runs = ["aniso_"+x for x in ["0","1","2","3","4"]]
    print(run_ids)
    print(runs)
    extract_sfr(run_ids,runs)
    
    
#     fig = P.figure()

