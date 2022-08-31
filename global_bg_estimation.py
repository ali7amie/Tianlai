from astropy.stats import sigma_clip
from astropy.stats import sigma_clipped_stats
import numpy as np


def global_bg_estimation(rectmap):

    clipped_map=sigma_clip(rectmap,sigma=3)
    stats=sigma_clipped_stats(rectmap)#mean,median,std

    change_ratio=stats[2]/np.std(rectmap)

    if change_ratio >= 0.8:
        crouded=False     
        return (stats[1],stats[2],clipped_map,crouded)

        
    else:
        crouded=True
        return (2.5*stats[1]-1.5*stats[0],stats[2],clipped_map,crouded) #bg,std

#def global_bg_estimation(map):
 #   '''
  #  This function take a map and give its background parameters ( med?? clipping??)

#    parameters
 #   -----------
 #   map : 2D numpy array

  #  Return
    #-------
   # statistics : tuple
   #              (map_median,map_std)
   # '''
    #self.clipmap = sigma_clip(self.rectmap,sigma=s)
    #self.clipstd2 = np.std(self.clipmap)
    #self.clipmed2 = np.median(self.clipmap)
   # map_median = np.median(map)
   # map_std = np.std(map)
   # return (map_median,map_std)


