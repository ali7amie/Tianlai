import numpy as np
import scipy.ndimage as nd
import k2jansky


def compute_flux(detection_kernels_size, barycenter_list, conv_maps, rectmap, freq, global_stat, max_baseline):

    ''' This function take the the convolution maps, the kernels size, the list of barycenters and the frequency, to compute a list of fluxes density in Kelvin and in Jansky
    
    parameters:
    -------------

    detection_kernels_size: tuple
                            (size of central kernel, size of intermediate kernel, size of peripherical kernel)

    barycenter_list:   tuples list     (2 x source number) 
                       each line present (vertical coor, horizontal coor) 
                       
    conv_maps: tuple
              (avg central map (ndarray), avg intermediate map (ndarray), avg peripherical map (ndarray), median peripherical map (ndarray), std peripherical map (ndarray))
                                      
    freq: float
          the frequency of signal in hertz


    Return:
    -------   
    flux_list:   tuples    (list,list) 
                       list of fluxes in kelvin and in Jansky                
    '''
    flux_integration_kernel=np.ones((detection_kernels_size[2],detection_kernels_size[2]))
    flux_map=nd.convolve(rectmap,flux_integration_kernel)-(detection_kernels_size[2]**2)*conv_maps[3]  # what to do with negative median - use farer and thiner ring
    #coor=np.transpose((barycenter_list[0],barycenter_list[1]))
    coor=( barycenter_list[:,0] , barycenter_list[:,1] )
    amp_list=rectmap[coor]- global_stat[0]
    flux_list_k=flux_map[coor]
    flux_list_jansky=k2jansky.k2jansky(use='freq',freq=freq,max_baseline=max_baseline,beam_surface=274*10**(-7),density_k=flux_list_k)
    
    return (flux_list_k,flux_list_jansky,amp_list)
