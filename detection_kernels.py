import scipy.ndimage as nd
import numpy as np
import numpy.ma as ma


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

    
    peripheric_median_map=np.zeros_like(map)
    peripheric_std_map = np.zeros_like(map)
    mask=np.array(peripheric_detection_kernel==0)
    for i in range(int(detection_kernels_size[2]/2),map.shape[0]-int(detection_kernels_size[2]/2)):
        for j in range(int(detection_kernels_size[2]/2),map.shape[1]-int(detection_kernels_size[2]/2)):
            if (detection_kernels_size[2]/2).is_integer()==False:
                kernel=map[i-int(detection_kernels_size[2]/2):i+(int(detection_kernels_size[2]/2)+1),j-int(detection_kernels_size[2]/2):j+(int(detection_kernels_size[2]/2)+1)]
                    
            else:
                    
                kernel=map[i-int(detection_kernels_size[2]/2):i+int(detection_kernels_size[2]/2),j-int(detection_kernels_size[2]/2):j+int(detection_kernels_size[2]/2)]
                    
            masked_kernel=ma.masked_array(kernel,mask)
            peripheric_std_map[i][j]=np.ma.std(masked_kernel)
            peripheric_median_map[i][j]=np.ma.median(masked_kernel)
    

 
    return ( [ central_detection_kernel,intermediate_detection_kernel,peripheric_detection_kernel ] , [avg_central_detection_kernel , avg_intermediate_detection_kernel,[ avg_peripheric_detection_kernel , peripheric_median_map , peripheric_std_map ]  ] )

def set_detection_criteria(detection_tools,n,global_bg_std):

    #detection criteria/condition
        
    threshold = detection_tools[1][2][1] + n * global_bg_std # above local background + n * global sigma
    #threshold = detection_tools[1][2][1] + n * detection_tools[1][2][2] # above local background + n * local sigma
    first_condition = detection_tools[1][0] >= threshold # the average of central kernel or the value itself
    detected_src_coor_upper_pixcorner = np.where( first_condition )
    
    return detected_src_coor_upper_pixcorner
