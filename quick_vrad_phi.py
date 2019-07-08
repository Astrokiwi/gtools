import matplotlib.pyplot as plt
import numpy as np
import gizmo_tools
from scipy import stats


run_id="3001"
run_dir="a2_e01"
snap_str="100"

phi_bin_edges = np.linspace(0.,20.,21)
phi_centres = (phi_bin_edges[:-1]+phi_bin_edges[1:])/2.

r_bin_edges = np.linspace(.7,1.5,5)
r_bin_centres = (r_bin_edges[:-1]+r_bin_edges[1:])/2.
n_r_bins = r_bin_centres.size

t,data=gizmo_tools.load_gizmo_pandas(run_id,run_dir,snap_str,["Masses","Coordinates","Velocities"])
print("loaded time=",t)
gizmo_tools.calculate_phi(data)
gizmo_tools.calculate_vrad(data)
data["rad3d"] = np.sqrt(data["Coordinates_x"]**2+data["Coordinates_y"]**2+data["Coordinates_z"]**2)
binstats = stats.binned_statistic_2d(data["phi"], data["rad3d"], data["vrad"], statistic='mean',bins=[phi_bin_edges,r_bin_edges])
bin_map = binstats[0]


fig,sp = plt.subplots(4,n_r_bins,figsize=(4.*n_r_bins,16.),sharex='col',sharey='row')


binstats = stats.binned_statistic_2d(data["phi"], data["rad3d"], data["vrad"], statistic='mean',bins=[phi_bin_edges,r_bin_edges])
bin_map = binstats[0]
for ix in range(n_r_bins):
    sp[0,ix].plot(phi_centres,bin_map[:,ix])
    
oneparsec_vel = bin_map[:,1]

binstats = stats.binned_statistic_2d(data["phi"], data["rad3d"], np.abs(data["Velocities_z"]), statistic='mean',bins=[phi_bin_edges,r_bin_edges])
bin_map = binstats[0]
for ix in range(n_r_bins):
    sp[1,ix].plot(phi_centres,bin_map[:,ix])
    
oneparsec_vz = bin_map[:,1]

binstats = stats.binned_statistic_2d(data["phi"], data["rad3d"], data["vrad"]*data["Masses"], statistic='sum',bins=[phi_bin_edges,r_bin_edges])
bin_map = binstats[0]
for ix in range(n_r_bins):
    sp[2,ix].plot(phi_centres,bin_map[:,ix])
    
oneparsec_mom = bin_map[:,1]

binstats = stats.binned_statistic_2d(data["phi"], data["rad3d"], data["Masses"], statistic='sum',bins=[phi_bin_edges,r_bin_edges])
bin_map = binstats[0]
for ix in range(n_r_bins):
    sp[3,ix].plot(phi_centres,bin_map[:,ix])

oneparsec_mass = bin_map[:,1]

fig.savefig("../figures/vrood.png",dpi=150)

mass_tot = np.nansum(oneparsec_mass)
vel_tot = np.nansum(oneparsec_mom)/mass_tot

oneparsec_mdot_profile = oneparsec_mass*(oneparsec_vel*1.0227e-6)/0.2 # n.b. converting velocity from km/s to pc/year

oneparsec_mdot = np.nansum(oneparsec_mdot_profile)
print(oneparsec_mdot,mass_tot,vel_tot)