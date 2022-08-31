import k2jansky
import numpy as np

def create_gaussian_src(src_size,src_std,src_flux_jansky,freq,max_baseline):

    ''' 
    This function simulate a single source having gaussian shape, it take the size, the std, the flux in jansky, and the observed frequency. It gives a 2D NumPy array as a source 
    It use the function k2jansky located in the file k2jansky.py to convert the flux density from kelvin to jansky

    parameters:
    -------------
    src_size: int
                
    src_std: float

    src_flux_jansky: float

    freq: float

    max_baseline: float
                          
       
    Return:
    -------  
    src: tuple (2D NumPy array, float, float, float)
         (a single source array, its flux in kelvin, its flux in jansky, its amplitude in kelvin)
    '''
          
    # create a source with amplitude of 1 k
    x,y = np.meshgrid(np.linspace(-(src_size//2),src_size//2,src_size),np.linspace(-(src_size//2),src_size//2,src_size))
    mean_x=0
    mean_y=0
    src_1=np.exp(-(x-mean_x)**2/(2*src_std**2)) * np.exp(-(y-mean_y)**2/(2*src_std**2))

    # Compute the flux density of the source having 1 K amplitude
    flux_kelvin_1=np.sum(src_1)

    # Converting this flux to jansky
    flux_jansky_1=k2jansky.k2jansky(use='freq',freq=freq,max_baseline=max_baseline,beam_surface=274*10**(-7),density_k=flux_kelvin_1)

    
    # create sources having src_flux_jansky Jansky 
    src=src_flux_jansky * (1 / flux_jansky_1) * np.exp(-(x-mean_x)**2/(2*src_std**2)) * np.exp(-(y-mean_y)**2/(2*src_std**2))

    # compute its flux density in kelvin
    flux_kelvin= np.sum(src)

    #compute the amplitude
    amp_kelvin=np.max(src)

    #compute its flux density in jansky
    flux_jansky=k2jansky.k2jansky(use='freq',freq=freq,max_baseline=max_baseline,beam_surface=274*10**(-7),density_k=flux_kelvin)


    return (src, flux_kelvin, flux_jansky, amp_kelvin)






    
