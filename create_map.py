import numpy as np
import create_gaussian_src
from photutils.datasets import make_noise_image


def create_map(map_size,src_number,src_size,src_std,src_flux_jansky,noise_std):
    """ This function create a map as a 2D Numpy array, enter map_size,src_number,src_size,src_std,src_flux_jansky,noise_std"""

    #create empty map
    rectmap=np.zeros((map_size,map_size))

    #create a gaussian source
    src=create_gaussian_src.create_gaussian_src(src_size,src_std,src_flux_jansky)

    #add src_number source to with random position, but with 2*src_size distant from the boundary.
    x=np.random.random_integers(low=2*src_size,high=map_size-2*src_size,size=src_number)
    y=np.random.random_integers(low=2*src_size,high=map_size-2*src_size,size=src_number) #generate random position, note that x and y here represent a matrix like coordinate, so the number of column or absice is y, the number of row or ordonnée is x

    
    
    for i in range(0,len(x)):
            rectmap[y[i]-(src_size//2):y[i]+((src_size//2)+1),x[i]-(src_size//2):x[i]+((src_size//2)+1)]=src
    
    #add gaussian noise
    noise =  make_noise_image((map_size,map_size), distribution='gaussian', mean=0,stddev=noise_std)
    rectmap=rectmap+noise
    catalog=np.transpose((y,x,src_flux_jansky*np.ones_like(x)))#upper pix_center, pix center to index. pix corner ??? plot??    
    return (rectmap,catalog)
    
    
    