from astropy.stats import sigma_clip
from astropy.stats import sigma_clipped_stats
import numpy as np

def global_bg_estimation(rectmap):
    #self.clipmap = sigma_clip(self.rectmap,sigma=s)
    #self.clipstd2 = np.std(self.clipmap)
    #self.clipmed2 = np.median(self.clipmap)
    map_median = np.median(rectmap)
    map_std = np.std(rectmap)
    return (map_median,map_std)