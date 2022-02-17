import scipy.ndimage as nd
import numpy as np


def pad_with(vector, pad_width, iaxis, kwargs):
    pad_value = kwargs.get('padder', 1)
    vector[:pad_width[0]] = pad_value
    vector[-pad_width[1]:] = pad_value


def create(map,detection_kernels_size):
    
    a=1/detection_kernels_size[0]**2
    central_detection_kernel=a*np.ones((detection_kernels_size[0],detection_kernels_size[0]))
    avg_central_detection_kernel=nd.convolve(map,central_detection_kernel)


    b=1/(detection_kernels_size[1]**2-detection_kernels_size[0]**2)
    intermediate_detection_kernel=b*np.pad(0*central_detection_kernel,int((detection_kernels_size[1]-detection_kernels_size[0])/2),pad_with)
    avg_intermediate_detection_kernel=nd.convolve(map,intermediate_detection_kernel)


    c=1/(detection_kernels_size[2]**2-detection_kernels_size[1]**2)
    peripheric_detection_kernel=c*np.pad(0*intermediate_detection_kernel,int((detection_kernels_size[2]-detection_kernels_size[1])/2),pad_with)
    avg_peripheric_detection_kernel=nd.convolve(map,peripheric_detection_kernel)
 
    return ([central_detection_kernel,intermediate_detection_kernel,peripheric_detection_kernel],[avg_central_detection_kernel,avg_intermediate_detection_kernel,avg_peripheric_detection_kernel])

def set_detection_criteria(avg_detection_kernels_list,n,global_bg_std):

    #detection criteria/condition
        
    threshold = avg_detection_kernels_list[2] + n * global_bg_std # above local background (could be median) + n * global sigma
    first_condition = avg_detection_kernels_list[0] >= threshold 
    detected_src_coor_upper_pixcorner = np.where( first_condition )
    
    return detected_src_coor_upper_pixcorner
