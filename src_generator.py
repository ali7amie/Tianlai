import numpy as np
import matplotlib.pyplot as plt
import scipy.constants

from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from photutils.aperture import CircularAperture
import math

def convert_K_to_Jansky(T):
    FREQ=1300*10**6
    #FREQ=1420*10**6
    c=3*10**8
    wavelenght=c/FREQ
    maxbaseline=16.5
    #beam_surf=6*10**(-5)
    #beam_surf=(wavelenght/maxbaseline)**2#33(58 arcmin**2 in rad) boltzmann (J.k-1)
    reso=12
    beam_surf=(reso*(1/60)*(math.pi/180))**2
    W=(2*scipy.constants.Boltzmann*beam_surf*T)/wavelenght**2
    S=W/10**(-26)
    return S

def pad_with(vector, pad_width, iaxis, kwargs):
    pad_value = kwargs.get('padder', 0)
    vector[:pad_width[0]] = pad_value
    vector[-pad_width[1]:] = pad_value

#def generate_xy(low,high,radius,srcnum):
   #while np.sqrt(x**2+y**2)<radius:
      # x=np.random.random_integers(low=low,high=high,size=srcnum)
     #  y=np.random.random_integers(low=low,high=high,size=srcnum)
   #return x,y

class generate:

    srccoef=10 #coeficient of src
    sigmag=1.5 # sigma of gaussian src
    sigmanoise=0.01 # sigma of gaussian noise
    meannoise=0 #mean of gaussian noise
    size=70 # size of map related to interferometer FOV
    srcnum=10
    srcsize=5
    side=10
    
    

    
    def __init__(self,srccoef) : 
        self.srccoef=srccoef 

    
    def set_src_number(self,srcnum): 
        self.srcnum=srcnum
    def set_src_coef(self,srccoef): 
        self.srccoef=srccoef
    def set_map_size(self,size): 
        self.size=size
        
    def set_sigma_noise(self,sigmanoise): 
        self.sigmanoise=sigmanoise

    def set_offset(self,meannoise):
        self.meannoise=meannoise

    def set_src_sigma(self,sigmag):
        self.sigmag=sigmag

    def set_side(self,side):
        self.side=side
        
    def set_src_size(self,srcsize):
        self.srcsize=srcsize

    def set_wavelenght(self,wavelenght):
        self.wavelenght=wavelenght

    def set_beam_surface(self,beamsurf):
        self.beamsurf=beamsurf

    def generate_position(self):
        self.x=np.random.random_integers(low=self.side,high=self.size-self.side,size=self.srcnum)
        self.y=np.random.random_integers(low=self.side,high=self.size-self.side,size=self.srcnum)
        #self.x,self.y=generate_xy(0,self.size,self.radius,self.srcnum)
        #self.x=[31]
        #self.y=[37]
        
        self.rectmap=np.zeros((self.size,self.size))
        for i in range(0,len(self.x)):
            self.rectmap[self.y[i]][self.x[i]]=1

    
    def generate_src_calibrate(self):
        gx, gy = np.meshgrid(np.linspace(-2,2,self.srcsize), np.linspace(-2,2,self.srcsize))
        mux=0
        muy=0
        #self.g = (1/248.69403730623083)*self.srccoef*np.exp(-( (gx-mux)**2 / ( 2.0 * self.sigmag**2 ) ) ) * np.exp(-( (gy-muy)**2 / ( 2.0 * self.sigmag**2 ) ) )
        self.g = np.exp(-( (gx-mux)**2 / ( 2.0 * self.sigmag**2 ) ) ) * np.exp(-( (gy-muy)**2 / ( 2.0 * self.sigmag**2 ) ) )
        #self.g = self.srccoef*np.exp(-( (gx-mux)**2 / ( 2.0 * self.sigmag**2 ) ) ) * np.exp(-( (gy-muy)**2 / ( 2.0 * self.sigmag**2 ) ) )


    def conversion(self):
        self.sumgK=np.sum(self.g)
        #self.wavelenght=0.2307692
        #self.beam_surf=0.00028465#33(58 arcmin**2 in sterad) boltzmann (J.k-1)
        #self.Watt=2*scipy.constants.Boltzmann*self.sumg*self.beam_surf/self.wavelenght**2
        #self.Jansky=self.Watt/(10**(-26))

        #self.FREQ=1300*10**6
        #self.c=3*10**8
        #self.wavelenght=self.c/self.FREQ
        #self.maxbaseline=22.5
        #self.beam_surf=(self.wavelenght/self.maxbaseline)**2#33(58 arcmin**2 in rad) boltzmann (J.k-1)
        #self.Watt=(2*scipy.constants.Boltzmann*self.sumgK*self.beam_surf)/self.wavelenght**2
        #self.sumgJ=self.Watt/(10**(-26))
        self.sumgJ=convert_K_to_Jansky(self.sumgK)

    def generate_src(self):
        gx, gy = np.meshgrid(np.linspace(-2,2,self.srcsize), np.linspace(-2,2,self.srcsize))
        mux=0
        muy=0
        #self.g = (1/248.69403730623083)*self.srccoef*np.exp(-( (gx-mux)**2 / ( 2.0 * self.sigmag**2 ) ) ) * np.exp(-( (gy-muy)**2 / ( 2.0 * self.sigmag**2 ) ) )
        self.gg = (1/self.sumgJ)*self.srccoef*np.exp(-( (gx-mux)**2 / ( 2.0 * self.sigmag**2 ) ) ) * np.exp(-( (gy-muy)**2 / ( 2.0 * self.sigmag**2 ) ) )
        #self.g = self.srccoef*np.exp(-( (gx-mux)**2 / ( 2.0 * self.sigmag**2 ) ) ) * np.exp(-( (gy-muy)**2 / ( 2.0 * self.sigmag**2 ) ) )
        self.sumggK=np.sum(self.gg)
        self.sumggJ=convert_K_to_Jansky(self.sumggK)
        
    def generate_map(self):
        
        self.h=np.where(self.rectmap==1)
        self.h_upper_pixcenter=self.h+np.full((2,np.shape(self.h)[1]),0.5)
        self.mod=np.sqrt(self.h_upper_pixcenter[0]**2+self.h_upper_pixcenter[1]**2)
        self.h_lower_pixcenter=[self.h_upper_pixcenter[0],self.size-self.h_upper_pixcenter[1]]
        self.h_center_pixcenter=self.h_lower_pixcenter+np.full(np.shape(self.h_upper_pixcenter),-self.size/2)
        self.transpose_h_center_pixcenter=np.transpose((-self.h_center_pixcenter[1],-self.h_center_pixcenter[0]))
        #self.mod=np.sqrt((-self.h_center_pixcenter[1])**2+(-self.h_center_pixcenter[0])**2)
        print('src number ={}'.format(np.shape(self.h)[1]))
        for k in range(0,np.shape(self.h)[1]):
            i=self.h[0][k]
            j=self.h[1][k]
            self.rectmap[i][j]=0
            self.rectmap[i-2:i+3,j-2:j+3]=self.rectmap[i-2:i+3,j-2:j+3]+self.gg

    def generate_noise(self):
        from photutils.datasets import make_noise_image
        shape = (self.size,self.size)
        self.noise =  make_noise_image(shape, distribution='gaussian', mean=self.meannoise,stddev=self.sigmanoise)
        #self.noise = np.pad(self.noise,5,pad_with)
    

    def choose_map(self):
        #self.rectmap=self.noise
        self.rectmap=self.rectmap+self.noise
        #self.rectmap=np.pad(self.rectmap,10,pad_with)

    def plot(self):
        
        self.transh=np.transpose((self.h[1],self.h[0]))
        src =  CircularAperture(self.transpose_h_center_pixcenter, r=0.5)
        plt.figure()
        
        plt.imshow(self.rectmap,extent=[-self.size/2., self.size/2., -self.size/2., self.size/2. ])
        src.plot(color='red',lw=1.5, alpha=0.5)
        plt.colorbar()

    def generator(self):
        self.generate_position()
        self.generate_src_calibrate()
        self.conversion()
        self.generate_src()
        self.generate_map()
        self.generate_noise()
        self.choose_map()
        #self.plot()
        
