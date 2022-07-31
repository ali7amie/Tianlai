import matplotlib.pyplot as plt

from photutils.aperture import CircularAperture

import numpy as np

import convert_map_index


def ploting_test(map,coor_center,bar_center,simulated_coor_center):

    
    detected_pixels_center = CircularAperture(coor_center, r=0.15)

    barycenter_center=CircularAperture(bar_center, r=0.4)

    simulated_src=CircularAperture(simulated_coor_center, r=0.3)

    
    plt.figure()
    
    plt.subplot(121)
    plt.imshow(map,extent=(-np.shape(map)[0]/2, np.shape(map)[0]/2, -np.shape(map)[0]/2, np.shape(map)[0]/2))
    plt.colorbar()
    detected_pixels_center.plot(color='red', lw=1.5, alpha=0.5)
    barycenter_center.plot(color='orange', lw=3, alpha=0.9)

    plt.subplot(122)
    plt.imshow(map,extent=(-np.shape(map)[0]/2, np.shape(map)[0]/2, -np.shape(map)[0]/2, np.shape(map)[0]/2))    
    plt.colorbar()
    barycenter_center.plot(color='red', lw=3, alpha=0.9)
    simulated_src.plot(color='orange', lw=3, alpha=0.9)
    plt.show(block=False)


