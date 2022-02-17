import numpy as np


def convert_upper_to_center(srccoor_upper_pixcorner,sizef,lenght):
    """ the default coordinate system of a map is matrix like, where A[x][y], 
    srccoor_upper_pixcenter=srccoor_upper_pixcorner+np.full((2,np.shape(srccoor_upper_pixcorner)[1]),0.5)
    srccoor_lower_pixcenter=[srccoor_upper_pixcenter[0],sizef-srccoor_upper_pixcenter[1]]
    srccoor_center_pixcenter=srccoor_lower_pixcenter+np.full(np.shape(srccoor_upper_pixcenter),-sizef/2)
    transpose_srccoor_center_pixcenter=np.transpose((-srccoor_center_pixcenter[1],-srccoor_center_pixcenter[0]))
    return [srccoor_upper_pixcorner,srccoor_upper_pixcenter,srccoor_lower_pixcenter,srccoor_center_pixcenter,transpose_srccoor_center_pixcenter]

def convert_center_to_upper(transpose_srccoor_center_pixcenter,sizef,lenght):
    #transpose_srccoor_center_pixcenter=np.transpose((np.zeros(sizef),np.zeros(sizef))
    srccoor_center_pixcenter=np.zeros((2,lenght))
    srccoor_center_pixcenter[1]=-transpose_srccoor_center_pixcenter[:,0]
    srccoor_center_pixcenter[0]=-transpose_srccoor_center_pixcenter[:,1]
    srccoor_lower_pixcenter=srccoor_center_pixcenter- np.full(lenght,-sizef/2)
    srccoor_upper_pixcenter=np.zeros((2,lenght))
    srccoor_upper_pixcenter[0]=srccoor_lower_pixcenter[0]
    srccoor_upper_pixcenter[1]=sizef-srccoor_lower_pixcenter[1]
    transpose_srccoor_upper_pixcenter=np.transpose((srccoor_upper_pixcenter[0],srccoor_upper_pixcenter[1]))
    srccoor_upper_pixcorner=srccoor_upper_pixcenter - np.full((2,lenght),0.5)
    transpose_srccoor_upper_pixcorner=np.transpose((srccoor_upper_pixcorner[0],srccoor_upper_pixcorner[1]))
    transpose_srccoor_upper_pixcenter=np.int_(transpose_srccoor_upper_pixcenter)
    return transpose_srccoor_upper_pixcenter