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

if __name__ == '__main__':
    run_id = "2014"
    runs = ["q2edd05redo","q2edd10_aniso1redo","q2edd10_aniso3redo","q2edd10redo","q2edd20redo","q2redo"]
    
#     fig = P.figure()

    gizmoDir = gizmo_tools.getGizmoDir()
    for irun,output_dir in enumerate(runs):
        print("Plotting: ",output_dir)
        fullDir = gizmoDir+"/"+run_id+"/"+output_dir
        infilename = fullDir+"/info.txt"
        print("Loading and parsing")
        sf_data_str = np.loadtxt(infilename,delimiter=",",usecols=(1,4),dtype=str)
#         sf_data_str = np.loadtxt(LineSkipper(open(infilename)),delimiter=",",usecols=(1,4),dtype=str)
#         sf_data_str = np.loadtxt(LineSkipper(open(infilename)),delimiter=",",dtype=str)
#         exit()
        nskip = 1000
        sf_data_str = sf_data_str[::nskip,:]
        
        times = np.loadtxt(sf_data_str[:,0],usecols=[1])
        sf_cumulative = np.loadtxt(sf_data_str[:,1],usecols=[1])
        
        times *= tunit
        sf_cumulative *= munit

        sfr = np.gradient(sf_cumulative)/np.gradient(times)
        
        times /= 1e6 # Myr for plot
        print("Plotting")
#         P.plot(times,sfr,label=runs[irun])
        
        output_data = np.vstack((times,sf_cumulative,sfr)).T
        np.savetxt("data/sf"+run_id+output_dir+".dat",output_data)
    
#     P.xlabel("$t$ (Myr)")
#     P.ylabel("SFR (M$_\odot$/yr)")
#     P.legend(loc='best')
#     P.yscale('log')
#     P.savefig("../figures/sfr_q2.png",dpi=200)
#     P.close()