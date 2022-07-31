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
import convert_map_index
import sys




class src_finder:

    '''

        This class take a map, and other additional attributes (4 or 5 in total) to create an instance of the src finder algorithm

        Attributes
        -----------
        use : str
              write 'python' if the input is a map is a simple numpy made with python
              write 'jskymap' if the map is puted in a FITS file and result from the map-making process, an additional argument should added the file_name
       
        fits_filename : str
                        FITS file name

        map_size : int
                   the length in pixels of a single side of the square map

        map_resolution : int
                         map resolution in arcminutes

        projection_center : tuple
                            (longitude, latitude) in degrees of the center of projection
                           
        detection_kernels_size : list
                                 [central kernel size, intermediate kernel size, peripherical kernel size] in pixels
        
        n : int
            threshold of detection in the condition : <central_kernel>   >   bg + n * std_bg

        aperture : list
                   [..,..] to be corrected

        Return
        --------
        ????????
                                 
    '''


    def __init__(self,use,detection_kernels_size,n,aperture,fits_filename='',map_size=80,map_resolution=18,projection_center=(90,180),noise=0.001,src_number=10,src_size=7,src_std=1.2,flux=1,freq=1300):

        self.use = use
        self.detection_kernels_size = detection_kernels_size
        self.n = n
        self.aperture = aperture
        self.map_size = map_size
        self.map_resolution = map_resolution
        self.projection_center = projection_center
        self.noise = noise
        self.src_number = src_number
        self.src_size = src_size
        self.src_std = src_std
        self.flux = flux
        self.freq = freq
            
        if self.use == 'python':
            self.simulated_map=create_map.create_map(self.map_size,self.src_number,self.src_size,self.src_std,self.flux,self.noise)
            self.rectmap=self.simulated_map[0]
            self.simulated_catalog=self.simulated_map[1]

        if self.use == 'JSkyMap':
            #create the instance, read a fits file, transform a spherical map into bumpy, extract a gnomonic projection.
            self.fits_filename = fits_filename
            self.spherical_map=hp.read_map(self.fits_filename,dtype=None)
            self.rectmap=hp.gnomview(self.spherical_map,rot=[self.projection_center[0],self.projection_center[1]],reso=self.map_resolution,xsize=self.map_size,ysize=self.map_size,return_projected_map=True,no_plot=True)
            #??????????????????? here you should import the txt file  catalog in case of jskymap
    
        #global background estimation, with and without clipping   
        self.global_bg_stat=global_bg_estimation.global_bg_estimation(self.rectmap)
 
        # create the 3 detection kernels + statistics of the peripheric kernel
        self.detection_tools=detection_kernels.create(self.rectmap,self.detection_kernels_size)

        # set the detection criteria
        self.detected_src_coor_upper_pixcorner = detection_kernels.set_detection_criteria(self.detection_tools,self.n,self.global_bg_stat[1])

        # resolving the no detection case
        if len(self.detected_src_coor_upper_pixcorner[0])!=0:
            self.no_detection=False
        else:
            self.no_detection=True
            
        

    def create_catalog(self):

        #compute the noise map of signal to noise 
        self.signal_to_noise_results = sort_signal_to_noise.compute_signal_to_noise_results(self.detection_tools,self.detected_src_coor_upper_pixcorner,self.rectmap)
    
        #hierarchical merging of neighbours pixels to form sources
        self.all_agglomerations = hierarchical_merging.hierarchical_merging(self.signal_to_noise_results[2],self.aperture)

        #compyte barycenters 
        self.barycenter_list=compute_barycenter.compute_barycenter(self.signal_to_noise_results[1],self.all_agglomerations)

        
            # compute flux
        self.flux_list=compute_flux.create_integration_kernel(self.rectmap,self.detection_kernels_size,self.detection_tools,self.barycenter_list)


        #convert upper to center
        self.simulated_coor_center=convert_map_index.convert_upper_to_center( [ self.simulated_catalog[:,0],self.simulated_catalog[:,1] ] , np.shape(self.rectmap)[0],len(self.simulated_catalog[:,0]) )[4]
        self.coor_center=convert_map_index.convert_upper_to_center(self.detected_src_coor_upper_pixcorner,np.shape(self.rectmap)[0],len(self.detected_src_coor_upper_pixcorner[0]))[4]
        self.bar_center=convert_map_index.convert_upper_to_center([self.barycenter_list[:,1],self.barycenter_list[:,0]],np.shape(self.rectmap)[0],len(self.barycenter_list[0]))[4]



        #convert pixel to world coordinate
        self.barycenter_list_world=pixel2world.pixel2world(self.projection_center,self.map_resolution,self.bar_center)
        self.simulated_src_world=pixel2world.pixel2world(self.projection_center,self.map_resolution,self.simulated_coor_center)

        #final catalog
        self.final_catalog=np.transpose((self.barycenter_list_world[:,0],self.barycenter_list_world[:,1],self.flux_list[1]))
        a=np.argsort(self.final_catalog[:,0])
        self.final_catalog=self.final_catalog[a]
        self.final_catalog=np.column_stack( ( range(0,len(self.final_catalog)) , self.final_catalog ))

        self.simulated_final_catalog=np.transpose((self.simulated_src_world[:,0],self.simulated_src_world[:,1],self.simulated_catalog[:,2]))
        b=np.argsort(self.simulated_final_catalog[:,0])
        self.simulated_final_catalog=self.simulated_final_catalog[b]
        self.simulated_final_catalog=np.column_stack( ( range(0,len(self.simulated_final_catalog)) , self.simulated_final_catalog ))



    def plot(self):
        ploting_test.ploting_test(self.rectmap,self.coor_center,self.bar_center,self.simulated_coor_center)


        
        
        

