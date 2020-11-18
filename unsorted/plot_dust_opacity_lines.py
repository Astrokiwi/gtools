import matplotlib.pyplot as plt
import gizmo_tools
import numpy as np

gizmo_tools.derive_opacities()


n_mu = 200

# fig,sp = plt.subplots(2,1,figsize=(6,6))
fig,sp = plt.subplots(figsize=(6,3))

ax = sp
log_mu = np.linspace(-7.7,7.3,n_mu)
mu = 10.**log_mu
ax.plot(mu,gizmo_tools.dust_opacity_function(mu))
ax.set_xlabel(r'$\lambda$ ($\mu$m)')
ax.set_ylabel(r'$\kappa_\lambda$ (cm$^{2}$ g$^{-1}$)')
ax.set_xscale('log')
ax.set_yscale('log')

# ax = sp[1]
# mu_local = 10.**np.linspace(0,2.95,n_mu)
# # mu_local = np.linspace(1,900,n_mu)
# ax.plot(mu_local,gizmo_tools.dust_opacity_function(mu_local))
# ax.set_yscale('log')
# ax.set_xscale('log')
# ax.set_xlabel(r'$\lambda$ ($\mu$m)')
# ax.set_ylabel(r'$\kappa_\lambda$ (cm$^{2}$ g$^{-1}$)')

# wavelengths_to_plot = wavelengths[2:] # skip the first two

for icol,(wavelength,line_name) in enumerate(
                                    zip(gizmo_tools.line_wavelengths,
                                        [gizmo_tools.lineNamesFull[x] for x in gizmo_tools.line_opacities])
                                    ):
    ax.axvline(wavelength,c=gizmo_tools.colors[icol],label=line_name,linewidth=1.)

# ax.legend(loc='best',fontsize='xx-small')
#     ax.text(wavelength,1.e5,line_name,rotation=90.)
fig.tight_layout()

fig.savefig("../figures/dust_opacity.pdf")
plt.close('all')