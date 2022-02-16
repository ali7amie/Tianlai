import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
import create_map
import global_bg_estimation
import detection_kernels
import sort_signal_to_noise

class src_finder:


    def __init__(self,use,fits_filename,map_size,map_resolution,projection_center,detection_kernels_size,n):
        if use == 'python':
            self.rectmap=create_map.create_map(100,10,5,1.5,1,0.001)
            self.detection_kernels_size=detection_kernels_size
            self.n=n
            self.map_size=np.shape(self.rectmap)[0]

        else:
            #create the instance, read a fits file, transform a spherical map into bumpy, extract a gnomonic projection.
            self.fits_filename = fits_filename
            self.map_resolution = map_resolution
            self.map_size = map_size
            self.detection_kernels_size=detection_kernels_size
            self.n=n
            self.projection_center = projection_center
            self.spherical_map=hp.read_map(self.fits_filename,dtype=None)
            self.rectmap=hp.gnomview(self.spherical_map,rot=[self.projection_center[0],self.projection_center[1]],reso=self.map_resolution,xsize=self.map_size,ysize=self.map_size,return_projected_map=True,no_plot=True)
    
        #global background estimation, with and without clipping   
        self.global_bg_stat=global_bg_estimation.global_bg_estimation(self.rectmap)
 
        # create the 3 detection kernels
        self.detection_tools=detection_kernels.create(self.rectmap,self.detection_kernels_size)
        # set the detection criteria
        self.detected_src_coor_upper_pixcorner = detection_kernels.set_detection_criteria(self.detection_tools[1],self.n,self.global_bg_stat[1])

        #compute the noise map of signal to noise (std of peripheric kernel)
        self.peripheric_std_map=sort_signal_to_noise.create_peripheric_std_map(self.rectmap,self.map_size,self.detection_kernels_size,self.detection_tools)
        self.signal_to_noise_results = sort_signal_to_noise.compute_signal_to_noise_results(self.detection_tools,self.peripheric_std_map,self.detected_src_coor_upper_pixcorner)
        
        
        

