import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys
import gizmo_tools
import itertools

testplot = False

run_id = "3032"
run_names = [
"longrun_medflow_vesc_defaultaniso_polar",
"newflow_vesc_thin_45"
]

run_parameters = gizmo_tools.load_run_parameters("3032")
gizmo_tools.run_parameters_names(run_parameters)

# run_parameters = {x:y for x,y in run_parameters.items() if x in run_names}


nruns = len(run_names)
snap_str = "200"

nbins = 400
dv_cent = 20.
vmin=-200.
vmax=-vmin
v_values = np.arange(vmin,vmax,(vmax-vmin)/nbins)

# lines = ["co1"]
lines = ["co1","co2","hcn1","hcn2","h2_1","h2_2","h2_3"]
# lines_selected = lines

line_order = [0,1,2,3,4,6,5]
lines_selected = [lines[x] for x in line_order]

# for poster:
# run_names = ["longrun_medflow_vesc_defaultaniso_polar"]
# lines_selected = ["co1","h2_3"]

line_codes = ["line_"+line for line in lines]


line_names_formatted = {}

for line in lines:
    words = gizmo_tools.lineNamesFull[line].rsplit(' ',2)
    formatted_name = " ".join(words[:-2])+"\n"+" ".join(words[-2:])
    line_names_formatted[line] = formatted_name

nray_strip = 41
nray = nray_strip*2

x_values = np.linspace(-20,20,nray_strip)

flux_norm = 1.e32

scale = .4

# angles = [0.,10.,20.]
angles = [10.]

m_smbh = 1.e6
smbh_soft = 0.01
m_hernquist = 1.e9
a_hernquist = 250.


def v_circ_smbh_hernquist(m_smbh,smbh_soft,m_hernquist,a_hernquist,r):
    G_msun_km2s2_pc = 0.004302052
    x = r/a_hernquist
    r2 = r**2
    m =   m_smbh*r2/(r2+smbh_soft**2) + m_hernquist * (x/(1.+x))**2
    v_circ = np.sqrt(G_msun_km2s2_pc * m/r)
    return v_circ
    

# GM = m_smbh*G_msun_km2s2_pc


v_circ = v_circ_smbh_hernquist(m_smbh,smbh_soft,m_hernquist,a_hernquist,np.abs(x_values))
v_circ_los = np.zeros_like(v_circ)

nx = len(lines_selected)*2+1
ny = len(run_names)*len(angles)

gap_size = 24#16

fw_inches=2*scale*nx
# fh_inches=3.*scale*ny # paper
fh_inches=4*scale*ny # paper
line_fig,line_ax = plt.subplots(ny,nx,figsize=(fw_inches,fh_inches),sharex=True,squeeze=False,gridspec_kw = {'width_ratios':[16]*len(lines_selected)+[gap_size]+[16]*len(lines_selected)})#,sharey=True)

# fw_inches=3.*scale*nx
# fh_inches=6.*scale*ny # poster
# line_fig,line_ax = plt.subplots(ny,nx,figsize=(fw_inches,fh_inches),sharex=True,squeeze=False)#,sharey=True)
# line_fig,line_ax = plt.subplots(ny,nx,figsize=(fw_inches,fh_inches),sharex=True,squeeze=False,gridspec_kw = {'width_ratios':[16]*len(lines_selected)+[8]+[16]*len(lines_selected)})#,sharey=True)

for irun,run_name in enumerate(run_names):
    ix=0
    for suffix in ["_ext"]:
        for iangle,rotate_phi in enumerate(angles):
            v_circ_los[x_values>0.] = v_circ[x_values>0.] * np.cos(rotate_phi) 
            v_circ_los[x_values<=0.]=-v_circ[x_values<=0.] * np.cos(rotate_phi)
            
            for iline,line in enumerate(lines_selected):
                
                
#                 file_name = "data/line_profs{}_{}_line_{}_{}deg.dat".format(suffix,run_name,line,rotate_phi)
                file_name = "data/line_profs{}_{}_line_{}_{}deg_{}.dat".format(suffix,run_name,line,rotate_phi,snap_str)
                if testplot:
                    vmaps = np.random.random((2,nbins//5,nray_strip//5))*flux_norm
                    print("testplot: ",file_name)
                else:
                    try:
                        ray_spectra = np.loadtxt(file_name)
                        print(file_name," loaded")
                    except OSError:
                        print("file not found - skipping")
                        continue

                    vmaps = np.zeros((2,nbins,nray_strip))

                    for iray in range(nray_strip):
                        for subplot_index in range(2):
                            offset = subplot_index*nray_strip
                            ray_spectrum = ray_spectra[:,iray+offset]

                            vmaps[subplot_index,:,iray] = ray_spectrum


                for subplot_index in range(2):
                    ix = iline+len(lines_selected)*subplot_index+subplot_index
                    iy = len(angles)*irun+iangle
                    
                    
                    lines_mappable=line_ax[iy,ix].imshow(        np.flip(vmaps[subplot_index,:,:])/flux_norm
                                                                ,norm=mpl.colors.LogNorm()
#                                                                 ,norm=mpl.colors.SymLogNorm(linthresh=1.e15)
                                                                ,vmin=1.e4/flux_norm
                                                                ,interpolation='nearest'
                                                                ,origin='lower'
                                                                ,extent=[x_values[0],x_values[-1],v_values[0],v_values[-1]]
                                                                ,aspect='auto'
                                                                ,cmap='viridis'
                                                                )
#                     if subplot_index==1:
#                         line_ax[iy,ix].plot(x_values,v_circ_los,lw=.5)
#                         line_ax[iy,ix].plot(x_values,v_circ_los*np.sqrt(2.),lw=.5)
                    
#                     line_ax[iy,ix].contour( np.flip(vmaps[subplot_index,:,:],1)/flux_norm
#                                             ,levels=np.logspace(4,32,5)/flux_norm
#                                             ,colors='white'
#                                             ,alpha=1
#                                             ,extent=[x_values[0],x_values[-1],v_values[0],v_values[-1]]
#                                             ,linewidths=.5)

# for iy,run_name in enumerate(run_names):
#     line_ax[iy,-1].yaxis.set_label_position("right")
#     line_ax[iy,-1].set_ylabel(run_parameters[run_name]['name'])

for iy in range(ny):
    line_ax[iy,0].set_ylabel(r"$v$ (km/s)")
    line_ax[iy,len(lines_selected)+1].set_ylabel(r"$v$ (km/s)")
    line_ax[iy,len(lines_selected)].set_visible(False)

for ix in itertools.chain(range(1,len(lines_selected)),range(len(lines_selected)+2,len(lines_selected)*2+1)):
    for iy in range(ny):
        line_ax[iy,ix].yaxis.set_visible(False)


for irun in range(len(run_names)):
    for iangle,angle in enumerate(angles):
        line_ax[iangle+len(angles)*irun,len(lines_selected)*2-1].yaxis.set_label_position("right")
        line_ax[iangle+len(angles)*irun,len(lines_selected)*2-1].set_ylabel(r"$\theta={:.0f}^\circ$".format(angle)+"\n{}".format(run_parameters[run_names[0]]['name']))

#     line_ax[iangle+len(angles),len(lines_selected)*2-1].yaxis.set_label_position("right")
#     line_ax[iangle+len(angles),len(lines_selected)*2-1].set_ylabel(r"$\theta={:.0f}^\circ$".format(angle)+"\n{}".format(run_parameters[run_names[1]]['name']))
    
# for ix in range(len(angles)):
for iline,line in enumerate(lines_selected):
    line_ax[0,iline].set_title(line_names_formatted[line],size='xx-small')
    line_ax[0,iline+len(lines_selected)+1].set_title(line_names_formatted[line],size='xx-small')

#     line_ax[0,ix].set_title("Vertical slit")
#     line_ax[0,ix+3].set_title("Horizontal slit")
#     line_ax[0,ix].set_title(r"$\theta={:.0f}^\circ$".format(angles[ix]))
#     line_ax[0,ix+len(angles)].set_title(r"$\theta={:.0f}^\circ$".format(angles[ix]))

    line_ax[-1,iline].set_xlabel(r"$z$ (pc)")
    line_ax[-1,iline+len(lines_selected)+1].set_xlabel(r"$x$ (pc)")

    line_ax[-1,iline].set_xticks([-20,-10,0,10])
    line_ax[-1,iline+len(lines_selected)+1].set_xticks([-20,-10,0,10])
    
    for ax in [line_ax[-1,iline],line_ax[-1,iline+len(lines_selected)+1]]:
        plt.sca(ax)
        plt.xticks(rotation=90)

#                 line_ax[1].set_ylabel(r"$v$ (km/s)")
#                 line_ax[1].set_xlabel(r"$x$ (pc)")
#                 line_ax[1].set_title("Horizontal slit")

# line_fig.subplots_adjust(wspace=0, hspace=0)
# line_fig.set_constrained_layout_pads(w_pad=0.,h_pad=0.,hspace=0.,wspace=0.)
# line_fig.tight_layout(w_pad=0.,h_pad=0.,pad=0.)

# PAPER
line_fig.subplots_adjust(bottom=0.25, top=0.9, left=0.1, right=0.8,
                    wspace=0., hspace=0.)

cb_ax = line_fig.add_axes([0.83, 0.25, 0.01, 0.65])

# POSTER
# line_fig.subplots_adjust(bottom=0.2, top=0.85, left=0.1, right=0.8,
#                     wspace=0., hspace=0.)
# cb_ax = line_fig.add_axes([0.83, 0.1, 0.02, 0.8])


cb = line_fig.colorbar(lines_mappable,cax=cb_ax)
cb.set_label('Flux (normalised)')

L = vmaps.shape[1]
pos = line_ax[0,0].get_position()
lpixx = (L/(pos.y1-pos.y0))
one_pixel_per_velbin_dpi = int(np.floor(lpixx/fh_inches))
print("dpi=",one_pixel_per_velbin_dpi)
line_fig.savefig("../figures/line_profs_collated_bigun_{}.png".format(snap_str),dpi=one_pixel_per_velbin_dpi)
# line_fig.savefig("../figures/line_profs_collated_bigun_poster_{}.pdf".format(snap_str),dpi=one_pixel_per_velbin_dpi)
plt.close('all')

