import numpy as np
import scipy.ndimage as nd
import k2jansky

def create_integration_kernel(map,detection_kernels_size,detection_tools,barycenter_list):
    flux_integration_kernel=np.ones((detection_kernels_size[2],detection_kernels_size[2]))
    flux_map=nd.convolve(map,flux_integration_kernel)-(detection_kernels_size[2]**2)*detection_tools[1][2][1]  # what to do with negative median - use farer and thiner ring
    #coor=np.transpose((barycenter_list[0],barycenter_list[1]))
    coor=[barycenter_list[:,0],barycenter_list[:,1]]
    flux_list_k=flux_map[coor]
    flux_list_jansky=k2jansky.k2jansky('freq',1300*10**6,16.5,1,flux_list_k)
     
    return (flux_list_k,flux_list_jansky)





#import cv2
#import numpy as np
#a=np.array([[2,2,2,2,2,2,2,2,2,2],[2,2,2,2,2,2,2,2,2,2],[2,2,2,2,2,2,2,2,2,2],[2,2,2,2,2,2,2,2,2,2],[2,2,2,2,2,2,2,2,2,2],[2,2,2,2,2,2,2,2,2,2],[2,2,2,2,2,2,2,2,2,2],[2,2,2,2,2,2,2,2,2,2],[2,2,2,2,2,2,2,2,2,2],[2,2,2,2,2,2,2,2,2,2]])
#b=np.array([[1,1,1],[1,1,1],[1,1,1]])
#c=cv2.filter2D(img, -1, kernel, borderType = cv2.BORDER_ISOLATED)