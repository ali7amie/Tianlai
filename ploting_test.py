import matplotlib.pyplot as plt

from photutils.aperture import CircularAperture

import numpy as np


def ploting_test(map,detected_src_coor_upper_pixcorner,barycenter_list):
    
    detected_pixels= CircularAperture(np.transpose((detected_src_coor_upper_pixcorner[1],detected_src_coor_upper_pixcorner[0])), r=0.15)
    barycenter=CircularAperture(barycenter_list, r=0.3)
 
    plt.figure()
    plt.imshow(map)
    plt.colorbar()
    detected_pixels.plot(color='red', lw=1.5, alpha=0.5)
    barycenter.plot(color='orange', lw=3, alpha=0.5)
    plt.show()
      