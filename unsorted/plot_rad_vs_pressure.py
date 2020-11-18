import matplotlib.pyplot as plt
import numpy as np

pc_unit = 1.e3
nbins = 100

if __name__ == '__main__':
    big_rp = np.loadtxt("data/rad_vs_pressure_3032_newflow_vesc_thin_45_100.txt",skiprows=1)
    small_rp = np.loadtxt("data/rad_vs_pressure_3001_a2_e01_100.txt",skiprows=1)
    
    max_acc = np.max((np.max(big_rp[:,1]),np.max(small_rp[:,1])))
#     big_rp[:,1]/=max_acc
#     small_rp[:,1]/=max_acc
    
    small_hist,xedges,yedges = np.histogram2d(np.log10(small_rp[:,0]*pc_unit),np.log10(small_rp[:,1]/max_acc),range=[[-1,2.5],[-11,0]],bins=nbins)
    big_hist,xedges,yedges = np.histogram2d(np.log10(big_rp[:,0]*pc_unit),np.log10(big_rp[:,1]/max_acc),range=[[-1,2.5],[-11,0]],bins=nbins)
    small_hist = np.arcsinh(small_hist.T)
    big_hist = np.arcsinh(big_hist.T)
    
    allmax = np.max((np.max(small_hist[np.isfinite(small_hist)]),(np.max(big_hist[np.isfinite(big_hist)]))))
    allmin = np.min((np.min(small_hist[np.isfinite(small_hist)]),(np.min(big_hist[np.isfinite(big_hist)]))))
    
    small_hist=(small_hist-allmin)/(allmax-allmin)
    big_hist=(big_hist-allmin)/(allmax-allmin)
    
    
    im = np.ones((nbins,nbins,3))
    im[:,:,0]=1.-big_hist
    im[:,:,1]=1.-small_hist-big_hist
    im[:,:,2]=1.-small_hist
#     im-=np.min(im[np.isfinite(im)])
#     im/=np.max(im[np.isfinite(im)])


    fig,sp = plt.subplots()
    ax = sp
    X,Y=np.meshgrid(xedges,yedges)
    
    ax.imshow(np.clip(im,0,1),extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],origin='lower',aspect=(xedges[-1]-xedges[-0])/(yedges[-1]-yedges[0]))
    ax.set_xlabel(r'$\log r$ (pc)')
    ax.set_ylabel(r'$\log a_r$ (normalised)')
#     ax.scatter(big_rp[:,0]/pc_unit,big_rp[:,1])
#     ax.scatter(small_rp[:,0]/pc_unit,small_rp[:,1])
#     ax.set_xscale('log')
#     ax.set_yscale('log')
    fig.tight_layout()
    fig.savefig("../figures/rad_vs_pressure.png",dpi=200)
    
    plt.close('all')