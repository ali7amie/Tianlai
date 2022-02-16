import numpy as np
import numpy.ma as ma


# the signal to noise ratio here is the ratio of the average of central kernel and the std of the peripheric kernel, so this function compute the std and later the ratio


def create_peripheric_std_map(map,map_size,detection_kernels_size,detection_tools):

    # computing std map

    peripheric_std_map = np.zeros((map_size,map_size))
    mask=np.array(detection_tools[0][2]==0)
    for i in range(int(detection_kernels_size[2]/2),map.shape[0]-int(detection_kernels_size[2]/2)):
        for j in range(int(detection_kernels_size[2]/2),map.shape[1]-int(detection_kernels_size[2]/2)):
            if (detection_kernels_size[2]/2).is_integer()==False:
                kernel=map[i-int(detection_kernels_size[2]/2):i+(int(detection_kernels_size[2]/2)+1),j-int(detection_kernels_size[2]/2):j+(int(detection_kernels_size[2]/2)+1)]
                    
            else:
                    
                kernel=map[i-int(detection_kernels_size[2]/2):i+int(detection_kernels_size[2]/2),j-int(detection_kernels_size[2]/2):j+int(detection_kernels_size[2]/2)]
                    
            masked_kernel=ma.masked_array(kernel,mask)
            std=np.std(masked_kernel)
            peripheric_std_map[i][j]=std/np.sqrt(detection_kernels_size[2]**2-detection_kernels_size[1]**2)
    return peripheric_std_map

def compute_signal_to_noise_results(detection_tools,peripheric_std_kernel,detected_src_coor_upper_pixcorner):
    
    #computing signal to noise ratio
    signal_to_noise_map = detection_tools[1][0]/peripheric_std_kernel
    signal_to_noise_list = signal_to_noise_map[detected_src_coor_upper_pixcorner]
    return (signal_to_noise_map,signal_to_noise_list)

    