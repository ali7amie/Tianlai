import matplotlib.pyplot as plt

from photutils.aperture import CircularAperture

import numpy as np

import convert_map_index


def ploting_test(map,detected_src_coor_upper_pixcorner,barycenter_list):

    coor_center=convert_map_index.convert_upper_to_center(detected_src_coor_upper_pixcorner,np.shape(map)[0],len(detected_src_coor_upper_pixcorner[0]))[4]
    bar_center=convert_map_index.convert_upper_to_center([barycenter_list[:,1],barycenter_list[:,0]],np.shape(map)[0],len(barycenter_list[0]))[4]
 
    
    detected_pixels = CircularAperture(np.transpose((detected_src_coor_upper_pixcorner[1],detected_src_coor_upper_pixcorner[0])), r=0.15)
    detected_pixels_center = CircularAperture(coor_center, r=0.15)

    barycenter_center=CircularAperture(bar_center, r=0.3)
    barycenter=CircularAperture(barycenter_list, r=0.3)

   #plot upper 

    
  #  plt.figure()
  #  plt.subplot(121)
  #  plt.title('upper_corner')
  #  plt.imshow(map)
  #  plt.colorbar()
  #  detected_pixels.plot(color='red', lw=1.5, alpha=0.5)
  #  barycenter.plot(color='orange', lw=3, alpha=0.5)
   
   
  #  plt.subplot(122)
  #  plt.title('center_corner')
  #  plt.imshow(map,extent=(-np.shape(map)[0]/2, np.shape(map)[0]/2, -np.shape(map)[0]/2, np.shape(map)[0]/2))
   # plt.colorbar()
   # detected_pixels_center.plot(color='red', lw=1.5, alpha=0.5)
   # barycenter_center.plot(color='orange', lw=3, alpha=0.5)
   # plt.show()


    
    plt.figure()
    plt.subplot(121)
    plt.imshow(map)
    plt.colorbar()
    detected_pixels.plot(color='red', lw=1.5, alpha=0.5)
    barycenter.plot(color='orange', lw=3, alpha=0.5)
    plt.subplot(122)
    plt.imshow(map,extent=(-np.shape(map)[0]/2, np.shape(map)[0]/2, -np.shape(map)[0]/2, np.shape(map)[0]/2))
    plt.colorbar()
    detected_pixels_center.plot(color='red', lw=1.5, alpha=0.5)
    barycenter_center.plot(color='orange', lw=3, alpha=0.5)
    plt.show()
    return (coor_center,bar_center)
      