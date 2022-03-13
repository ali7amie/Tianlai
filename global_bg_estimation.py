from astropy.stats import sigma_clip
from astropy.stats import sigma_clipped_stats
import numpy as np

def global_bg_estimation(map):
    '''
    This function take a map and give its background parameters ( med?? clipping??)

    parameters
    -----------
    map : 2D numpy array

    Return
    -------
    statistics : tuple
                 (map_median,map_std)
    '''
    #self.clipmap = sigma_clip(self.rectmap,sigma=s)
    #self.clipstd2 = np.std(self.clipmap)
    #self.clipmed2 = np.median(self.clipmap)
    map_median = np.median(map)
    map_std = np.std(map)
    return (map_median,map_std)