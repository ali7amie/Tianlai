
import src_finder
import numpy as np
import cross_matching
import matplotlib.pyplot as plt

class post_processing:



    def __init__(self , use , detection_kernels_size , aperture , freq_list , noise_list , flux_list , n_list, iteration_number , fits_filename='' , map_size=80 , map_resolution=18 , projection_center=(90,180) , src_number=10 , src_size=7 , src_std=1.2, matching_aperture=0.5 , max_distance=0.08):


        self.use = use
        self.detection_kernels_size = detection_kernels_size
        self.n_list = n_list
        self.aperture = aperture
        self.map_size = map_size
        self.map_resolution = map_resolution
        self.projection_center = projection_center
        self.noise_list = noise_list
        self.src_number = src_number
        self.src_size = src_size
        self.src_std = src_std
        self.flux_list = flux_list
        self.freq_list = freq_list
        self.iteration_number = iteration_number
        self.matching_aperture = matching_aperture
        self.max_distance = max_distance


        self.maps=np.zeros((len(freq_list),len(self.noise_list),len(self.n_list),len(self.flux_list),self.iteration_number),dtype=tuple)#dtype=src_finder.src_finder)

        self.median_efficiency=np.zeros((len(freq_list),len(self.noise_list),len(self.n_list),len(self.flux_list)))
        self.std_efficiency=np.zeros((len(freq_list),len(self.noise_list),len(self.n_list),len(self.flux_list)))
        self.median_spurious=np.zeros((len(freq_list),len(self.noise_list),len(self.n_list),len(self.flux_list)))
        self.std_spurious=np.zeros((len(freq_list),len(self.noise_list),len(self.n_list),len(self.flux_list)))




        for i in range(0,len(self.freq_list)):
            for j in range(0,len(self.noise_list)):
                for k in range(0,len(self.n_list)):
                    for l in range(0,len(self.flux_list)):

                        efficiency_vector=np.zeros(self.iteration_number)
                        spurious_vector=np.zeros(self.iteration_number)

                        for m in range(0,self.iteration_number):
                    
                            
                            self.each_map=src_finder.src_finder('python' , [3,5,7] , self.n_list[k] , self.aperture , '' , self.map_size , self.map_resolution , self.projection_center , self.noise_list[j], self.src_number , self.src_size , self.src_std , self.flux_list[l] , self.freq_list[i]) #instantiation of the detector on each map

                            self.matching_results=cross_matching.cross_matching(self.each_map.final_catalog,self.each_map.simulated_final_catalog,self.matching_aperture,self.max_distance)
                            self.each_efficiency=len(self.matching_results[0])/len(self.each_map.simulated_final_catalog)
                            self.each_spurious=(len(self.each_map.final_catalog)-len(self.matching_results[0]))/len(self.each_map.final_catalog)

                            efficiency_vector[m]=self.each_efficiency
                            spurious_vector[m]=self.each_spurious
                        
                            self.maps[i][j][k][l][m]=[self.each_map,self.matching_results,self.each_efficiency,self.each_spurious]  # complete data for each map
                            print('*')  

                        self.median_efficiency[i][j][k][l]=np.median(efficiency_vector)
                        self.std_efficiency[i][j][k][k]=np.std(efficiency_vector)
                        self.median_spurious[i][j][k][l]=np.median(spurious_vector)
                        self.std_spurious[i][j][k][l]=np.std(spurious_vector)

       

                        
  
        for i in range(0,len(self.freq_list)):
            for j in range(0,len(self.noise_list)):
                plt.figure()
                for k in range(0,len(self.n_list)):
                    plt.plot(self.flux_list,self.median_efficiency[i][j][k])
                plt.legend(loc='lower right')
                plt.show()

        
