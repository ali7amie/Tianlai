import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
import create_map
import global_bg_estimation
import detection_kernels
import sort_signal_to_noise
import hierarchical_merging
import compute_barycenter
import ploting_test
import compute_flux
import pixel2world



class src_finder:


    def __init__(self,use,fits_filename,map_size,map_resolution,projection_center,detection_kernels_size,n,aperture):
        if use == 'python':
            self.create=create_map.create_map(100,10,7,1.5,1,0.001)
            self.rectmap=self.create[0]
            self.catalog=self.create[1]
            self.detection_kernels_size=detection_kernels_size
            self.n=n
            self.map_size=np.shape(self.rectmap)[0]
            self.aperture=aperture
            self.projection_center=projection_center
            self.map_resolution=map_resolution

        else:
            #create the instance, read a fits file, transform a spherical map into bumpy, extract a gnomonic projection.
            self.fits_filename = fits_filename
            self.map_resolution = map_resolution
            self.map_size = map_size
            self.detection_kernels_size=detection_kernels_size
            self.n=n
            self.aperture=aperture
            self.projection_center = projection_center
            self.spherical_map=hp.read_map(self.fits_filename,dtype=None)
            self.rectmap=hp.gnomview(self.spherical_map,rot=[self.projection_center[0],self.projection_center[1]],reso=self.map_resolution,xsize=self.map_size,ysize=self.map_size,return_projected_map=True,no_plot=True)
    
        #global background estimation, with and without clipping   
        self.global_bg_stat=global_bg_estimation.global_bg_estimation(self.rectmap)
 
        # create the 3 detection kernels + statistics of the peripheric kernel
        self.detection_tools=detection_kernels.create(self.rectmap,self.detection_kernels_size)
        # set the detection criteria
        self.detected_src_coor_upper_pixcorner = detection_kernels.set_detection_criteria(self.detection_tools,self.n,self.global_bg_stat[1])

        #compute the noise map of signal to noise 
        self.signal_to_noise_results = sort_signal_to_noise.compute_signal_to_noise_results(self.detection_tools,self.detected_src_coor_upper_pixcorner,self.rectmap)

        #hierarchical merging of neighbours pixels to form sources
        self.all_agglomerations = hierarchical_merging.hierarchical_merging(self.signal_to_noise_results[2],self.aperture)

        #compyte barycenters 
        self.barycenter_list=compute_barycenter.compute_barycenter(self.signal_to_noise_results[1],self.all_agglomerations)



        # compute flux
        self.flux_list=compute_flux.create_integration_kernel(self.rectmap,self.detection_kernels_size,self.detection_tools,self.barycenter_list)

        # a ploting test
        self.plot=ploting_test.ploting_test(self.rectmap,self.detected_src_coor_upper_pixcorner,self.barycenter_list)

        #convert pixel to world coordinate
        self.barycenter_list_world=pixel2world.pixel2world(self.projection_center,self.map_resolution,self.barycenter_list)

        #final catalog
        self.final_catalog=np.transpose((self.barycenter_list_world[:,0],self.barycenter_list_world[:,1],self.flux_list[1]))


        
        
        

