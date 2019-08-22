import numpy as np
import matplotlib.pyplot as plt
# from scipy import stats
import sys
# from astropy import modeling
# from scipy import signal
import pandas as pd

run_id = "3032"
run_names = [
"longrun_medflow_settled_defaultaniso_polar",
"longrun_medflow_vesc_defaultaniso",
"longrun_medflow_vesc_defaultaniso_polar",
"longrun_weakflow_rapid_defaultaniso",
"longrun_weakflow_rapid_defaultaniso_polar",
"longrun_weakflow_settled_defaultaniso",
"longrun_weakflow_settled_defaultaniso_polar",
"longrun_weakflow_vesc_defaultaniso",
"longrun_weakflow_vesc_defaultaniso_polar",
"newflow_settled_thin_up",
"newflow_vesc_thin_45",
"newflow_vesc_thin_side",
"newflow_vesc_thin_up"]

nruns = len(run_names)
snap_str = "100"

nbins = 400
dv_cent = 20.
vmin=-1000.
vmax=-vmin
v_values = np.arange(vmin,vmax,(vmax-vmin)/nbins)

lines = ["co1","co2","hcn1","hcn2","h2_1","h2_2","h2_3"]
line_codes = ["line_"+line for line in lines]

nray_strip = 41
nray = nray_strip*2

# max_clip = 10**20

max_clip = 21.
min_clip = 0.

features = ["left","centre","right"]
line_properties = ["strength","width","v"]
line_columns = [x+"_"+y for x in features for y in line_properties]

# gaussian_model = modeling.models.Gaussian1D(1.e10,-200,50)+modeling.models.Gaussian1D(1.e10,200,50)+modeling.models.Gaussian1D(1.e20,0,10)
# gaussian_fit = modeling.fitting.LevMarLSQFitter()

# def find_peaks_stepthrough(y,v):
#     i = 0
#     
#     initial_peaks,peak_dict = signal.find_peaks(y,prominence=0.7) # i.e. factor of 5 minimum prominence
    

# def contiguous_regions(condition):
#     """Finds contiguous True regions of the boolean array "condition". Returns
#     a 2D array where the first column is the start index of the region and the
#     second column is the end index.
#     
#     From: https://stackoverflow.com/questions/4494404/find-large-number-of-consecutive-values-fulfilling-condition-in-a-numpy-array"""
# 
#     # Find the indicies of changes in "condition"
#     d = np.diff(condition)
#     idx, = d.nonzero() 
# 
#     # We need to start things after the change in "condition". Therefore, 
#     # we'll shift the index by 1 to the right.
#     idx += 1
# 
#     if condition[0]:
#         # If the start of condition is True prepend a 0
#         idx = np.r_[0, idx]
# 
#     if condition[-1]:
#         # If the end of condition is True, append the length of the array
#         idx = np.r_[idx, condition.size] # Edit
# 
#     # Reshape the result into two columns
#     idx.shape = (-1,2)
#     return idx

def line_info(y,slice,x):
    strength = np.trapz(y[slice],x[slice])
    peakwidth = strength/np.max(y[slice])
    mean_x = np.sum(y[slice]*x[slice])/np.sum(y[slice])
    return strength,peakwidth,mean_x

def line_info_preset_x(x):
    def f(y,slice):
        return line_info(y,slice,x)
    return f

line_info_v = line_info_preset_x(v_values)

def tophat(x, base_level, hat_level, hat_mid, hat_width):
    return np.where((hat_mid-hat_width/2. < x) & (x < hat_mid+hat_width/2.), hat_level, base_level)

def triangle(x, base_level, tri_peak, tri_mid, tri_width):
    return np.where((tri_mid-tri_width/2. < x) & (x < tri_mid+tri_width/2.), tri_peak-2*np.abs(x-tri_mid)*(tri_peak-base_level)/tri_width, base_level)

def three_peaks(x,left_level,left_peak,left_mid,left_width,cent_level,cent_peak,cent_mid,cent_width,right_level,right_peak,right_mid,right_width):
    return tophat(x,left_level,left_peak,left_mid,left_width)+tophat(x,right_level,right_peak,right_mid,right_width)+triangle(x,cent_level,cent_peak,cent_mid,cent_width)

for run_name in ["newflow_vesc_thin_side_line"]:
    for line in lines:
#     for line in ["co1"]:
        line_fig,line_ax = plt.subplots(1,2,figsize=(8.,24.))
#         fig,ax = plt.subplots(1,2,figsize=(12.,4.))
        ray_spectra = np.loadtxt("data/line_profs_{}_{}.dat".format(run_name,line))
#         mask = np.ones(ray_spectra.shape,dtype=np.bool)
#         mask[197:204,:] = False
#         spectra_masked_centre = np.ma.masked_array(ray_spectra,mask=mask)
#         ray_spectra[197:204,:] = 0. # CUT OUT THE DISC
#         print(np.percentile(ray_spectra,100.*(1.-1./nray)),np.percentile(ray_spectra,100.*(1.-1./nray_strip)),np.ceil(np.log10(np.percentile(ray_spectra,100.*(1.-1./nray)))))
#         max_clip = np.ceil(np.log10(np.percentile(ray_spectra,100.*(1.-1./nbins))))
#         max_clip = np.log10(np.max(spectra_masked_centre))
        max_clip = np.log10(np.max(ray_spectra))
#         max_clip = np.max(ray_spectra)
    
        
        df = [pd.DataFrame(columns=line_columns),pd.DataFrame(columns=line_columns)]
        for iray in range(nray_strip):
#         for iray in [5]:
#             ax.plot(v_values,np.clip(ray_spectra[:,iray],0,10**13)/np.max(ray_spectra[:,iray])-1.2*iray)
            for subplot_index in range(2):
#             for subplot_index in [0]:
                offset = subplot_index*nray_strip
                ray_spectrum = ray_spectra[:,iray+offset]
                log_spectrum = np.log10(ray_spectrum)
#                 log_spectrum = ray_spectra[:,iray+offset]
                line_ax[subplot_index].plot(v_values,np.clip(log_spectrum,min_clip,max_clip)/(max_clip-min_clip)-1.2*iray)
#                 ax[subplot_index].plot(v_values,log_spectrum)
#                 ax[subplot_index].plot(v_values,ray_spectrum)
# Try this "algorithm": remove central peak, get mean & weighted standard deviation of left, and of right?
#                 threshold = min_clip
#                 above_threshold = log_spectrum>min_clip
#                 threshold_regions = contiguous_regions(above_threshold)
                central_region = [-dv_cent,dv_cent]
#                 for threshold_region in threshold_regions:
#                     v_region = v_values[threshold_region]
#                     if v_region[0]<0. and v_region[1]>0.:
#                         central_region = v_region
#                 if central_region is None:
#                     central_slice = np.zeros(v_values.shape,dtype=np.bool)
#                 else:
                central_slice = (v_values>=central_region[0]) & (v_values<=central_region[1])

                left_slice = (v_values<0.) & ~central_slice
                right_slice = (v_values>0.) & ~central_slice

                # find central peak (if any)
                left_peakwidth,mean_left,right_peakwidth,mean_right,central_peakwidth,mean_central = [0.]*6

                if np.sum(left_slice)>0.:
                    left_line = line_info_v(ray_spectrum,left_slice)
                if np.sum(right_slice)>0.:
                    right_line = line_info_v(ray_spectrum,right_slice)
                if np.sum(central_slice)>0.:
                    central_line = line_info_v(ray_spectrum,central_slice)
                
                df[subplot_index]=df[subplot_index].append({x:y for x,y in zip(line_columns,left_line+central_line+right_line)},ignore_index=True)

#                 print(subplot_index,left_line,central_line,right_line)

                
#         ax[0].set_xlim([-500,500])
#         ax[1].set_xlim([-500,500])
        line_ax[0].set_xlabel(r"$v$ (km/s)")
        line_ax[1].set_xlabel(r"$v$ (km/s)")
#         ax[0].set_ylim([0.,None])
#         ax[1].set_ylim([0.,None])
        line_fig.savefig("pics/line_profs_{}_{}.pdf".format(run_name,line))
        plt.close('all')
        
        
        summary_fig,summary_ax = plt.subplots(2,3,figsize=(9.,6.))
        for subplot_index in range(2):
            for feature in ["left","centre","right"]:
                summary_ax[subplot_index,0].plot(np.log10(df[subplot_index][feature+"_strength"]),label=feature)
                summary_ax[subplot_index,1].plot(df[subplot_index][feature+"_v"],label=feature)
                summary_ax[subplot_index,2].plot(df[subplot_index][feature+"_width"],label=feature)
            summary_ax[subplot_index,0].set_ylim(bottom=0.)
        summary_ax[0,0].legend(loc='best')
        summary_fig.savefig("pics/line_summary_{}_{}.pdf".format(run_name,line))
        plt.close('all')
        
