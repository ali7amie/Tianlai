
import src_finder
import numpy as np
import cross_matching
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import k2jansky
import math


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
        self.efficiency=np.zeros((len(freq_list),len(self.noise_list),len(self.n_list),len(self.flux_list),self.iteration_number))


        s=0
        ss=len(self.freq_list)*len(self.noise_list)*len(self.n_list)*len(self.flux_list)*self.iteration_number

        for i in range(0,len(self.freq_list)):
            for j in range(0,len(self.noise_list)):
                for k in range(0,len(self.n_list)):
                    for l in range(0,len(self.flux_list)):

                        self.efficiency_vector=np.zeros(self.iteration_number)
                        self.spurious_vector=np.zeros(self.iteration_number)

                        for m in range(0,self.iteration_number):
                    
                            s=s+1
                            self.each_map=src_finder.src_finder('python' , [3,5,7] , self.n_list[k] , self.aperture , '' , self.map_size , self.map_resolution , self.projection_center , self.noise_list[j], self.src_number , self.src_size , self.src_std , self.flux_list[l] , self.freq_list[i]) #instantiation of the detector on each map

                            
                            if self.each_map.no_detection==False:
                                self.each_map.create_catalog()
                                self.matching_results=cross_matching.cross_matching(self.each_map.final_catalog,self.each_map.simulated_final_catalog,self.matching_aperture,self.max_distance)
                                self.each_efficiency=len(self.matching_results[0])/len(self.each_map.simulated_final_catalog)
                                self.each_spurious=(len(self.each_map.final_catalog)-len(self.matching_results[0]))/len(self.each_map.final_catalog)
                                print(ss-s)
 
                            if self.each_map.no_detection==True:
                                print('no detection')
                                self.each_efficiency=0
                                self.each_spurious=0

                            self.efficiency[i][j][k][l][m]=self.each_efficiency

                            self.efficiency_vector[m]=self.each_efficiency
                            self.spurious_vector[m]=self.each_spurious
                        
                            self.maps[i][j][k][l][m]=[self.each_map,self.matching_results,self.each_efficiency,self.each_spurious]  # complete data for each map

                        self.median_efficiency[i][j][k][l]=np.median(self.efficiency_vector)
                        self.std_efficiency[i][j][k][l]=np.std(self.efficiency_vector)
                        self.median_spurious[i][j][k][l]=np.median(self.spurious_vector)
                        self.std_spurious[i][j][k][l]=np.std(self.spurious_vector)

       

# plot efficiency/spurious/flux to characterise the algorithm                       
  
        for i in range(0,len(self.freq_list)):
            for j in range(0,len(self.noise_list)):
                fig, ax = plt.subplots(constrained_layout=True)
                for k in range(0,len(self.n_list)):
                    #plt.style.use('seaborn')
                    plt.title('Efficiency as a function of Flux | Freq={}'.format(self.freq_list[i]))
                    x=self.flux_list
                    y=self.median_efficiency[i][j][k]
                    bar=self.std_efficiency[i][j][k]

                    ax.plot(x,y,label='n={}'.format(self.n_list[k]))
                    ax.set_xlabel('Flux density [Jy]')
                    ax.set_ylabel('Efficiency')
                    ax.legend(loc='lower right')
                    ax.grid(True)
                    
                    ax2=ax.twiny()          
                    ax2.set_xlabel('S/N')
                    ax2.set_xticks((1/(2*math.pi*self.src_std**2)) * k2jansky.jansky2k('freq',1300*1000000,16,0,self.flux_list)/self.noise_list[j])
                    ax.errorbar(x,y, yerr=bar, ls='none', fmt='o', elinewidth=1, capsize=1, barsabove=False)
                    

                    ax3=ax.twiny()
                    ax3.spines['top'].set_position(('outward', 40))
                    ax3.set_xticks( (1/(2*math.pi*self.src_std**2)) * k2jansky.jansky2k('freq',1300*1000000,16,0,self.flux_list))
                    ax3.set_xlabel('Amplitude [mK]')

                    plt.show(block=False)

# plot efficiency boxplot
  
        for i in range(0,len(self.freq_list)):
            for j in range(0,len(self.noise_list)):
                fig, ax = plt.subplots(constrained_layout=True)
                for k in range(0,len(self.n_list)):
                    plt.boxplot(np.transpose(self.efficiency[i][j][k]) )
                plt.legend(loc='lower right')
                plt.show(block=False)

