import k2jansky
import numpy as np


def create_gaussian_src(src_size,src_std,src_flux_jansky):
    # create a source with amplitude of 1 k
    x,y = np.meshgrid(np.linspace(-(src_size//2),src_size//2,src_size),np.linspace(-(src_size//2),src_size//2,src_size))
    mean_x=0
    mean_y=0
    src=np.exp(-(x-mean_x)**2/(2*src_std**2)) * np.exp(-(y-mean_y)**2/(2*src_std**2))
    flux_kelvin_of_1k_amplitude=np.sum(src)
    flux_jansky_of_1k_amplitude=k2jansky.k2jansky('freq',1300*10**6,16.5,1,flux_kelvin_of_1k_amplitude)
    # add the src_flux_jansky factor
    src=(src_flux_jansky/flux_jansky_of_1k_amplitude)*np.exp(-(x-mean_x)**2/(2*src_std**2)) * np.exp(-(y-mean_y)**2/(2*src_std**2))
    return src