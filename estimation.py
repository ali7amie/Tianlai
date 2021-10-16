import extract_class
import importlib
importlib.reload(extract_class)
import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.constants
import statistics
import healpy as hp
import src_generator

def convert_Jansky_to_K(S):
    FREQ=1300*10**6
    #FREQ=1420*10**6
    c=3*(10**8)
    wavelenght=c/FREQ
    maxbaseline=16.5
    #beam_surf=6*10**(-5)
    #beam_surf= 1.21*10**(-5)# pixel of 12 arcmin  instrument
    reso=12
    beam_surf=(reso*(1/60)*(math.pi/180))**2
    #beam_surf=4*10**(-6)   #healpix 512
    #beam_surf=1.59*10**(-5) #healpix 256
    #beam_surf=(wavelenght/maxbaseline)**2#33(58 arcmin**2 in rad) boltzmann (J.k-1)
    T=( (10**(-26))*S*wavelenght**2)/(2*scipy.constants.Boltzmann*beam_surf)
    return T

# convert kelvin to Jansky
def convert_K_to_Jansky(T):
    FREQ=1300*10**6
    #FREQ=1420*10**6
    c=3*(10**8)
    wavelenght=c/FREQ
    maxbaseline=16.5
    #beam_surf=6*10**(-5)
    #beam_surf= 1.21*10**(-5)
    #beam_surf=(wavelenght/maxbaseline)**2# 1.94 10-4
    reso=12
    beam_surf=(reso*(1/60)*(math.pi/180))**2
    #beam_surf=4*10**(-6)
    #beam_surf=1.59*10**(-5)
    W=(2*scipy.constants.Boltzmann*beam_surf*T)/(wavelenght**2)
    S=W/10**(-26)
    return S

a=extract_class.extract()
a.set_filename('filtmap_noise_256_0.001.fits')
a.set_size(90)
a.set_detuse('skymaps')
a.set_n_sigma(5)
a.set_query(1)
a.set_k(0)
a.set_query_disk_radius(35)
a.set_freq(1300*10**6)
a.set_maxbaseline(16.5)
a.set_reso(12)
a.set_nker1(3)
a.set_nker2(5)
a.set_radius(np.sqrt(2*(2.5**2)))
a.set_nker3(7)
a.set_croping_radius(35)
a.extract()
a.detect()

c=extract_class.extract()
c.set_filename('filtmap_noise_256.fits')
c.set_size(90)
c.set_detuse('skymaps')
c.set_n_sigma(5)
c.set_query(1)
c.set_k(0)
c.set_query_disk_radius(35)
c.set_freq(1300*10**6)
c.set_maxbaseline(16.5)
c.set_reso(12)
c.set_nker1(3)
c.set_nker2(5)
c.set_radius(np.sqrt(2*(2.5**2)))
c.set_nker3(7)
c.extract()
c.detect()

b=src_generator.generate(5)
b.set_src_number(1)
b.set_map_size(90)
b.set_sigma_noise(0)
b.set_src_sigma(1.5)
b.set_src_size(5)
b.generator()




#plt.figure()
#plt.imshow(a.rectmap)
#plt.colorbar()
#plt.show()


#plt.figure()
#plt.imshow(b.rectmap)
#plt.colorbar()
#plt.show()

#plt.figure()
#plt.plot(a.rectmap[37])

#plt.figure()
#plt.plot(a.rectmap[42])

#plt.figure()
#plt.plot(a.rectmap[53])


%matplotlib
import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
import numpy.ma as ma
from photutils.aperture import CircularAperture
from photutils.aperture import RectangularAperture
from mpl_toolkits.mplot3d import Axes3D


nker1=3
nker2=5
nker3=7


realmap=a.rectmap
rectmap=b.rectmap
realmap2=c.rectmap



coor=[31,37]
src= CircularAperture(coor,r=0.25 )
ker1=RectangularAperture(coor,nker1,nker1)
ker2=RectangularAperture(coor,nker2,nker2)
ker3=RectangularAperture(coor,nker3,nker3)


plt.figure()
plt.subplot(131)
plt.imshow(realmap)
plt.title('real source 512')
src.plot(color='red',lw=1.5, alpha=0.5)
#ker1.plot(color='red',lw=1.5, alpha=0.5)
#ker2.plot(color='yellow',lw=1.5, alpha=0.5)
#ker3.plot(color='orange',lw=1.5, alpha=0.5)


plt.subplot(132)
plt.imshow(rectmap)
plt.title('Ideal guassian source')
src.plot(color='red',lw=1.5, alpha=0.5)
#ker1.plot(color='red',lw=1.5, alpha=0.5)
#ker2.plot(color='yellow',lw=1.5, alpha=0.5)
#ker3.plot(color='orange',lw=1.5, alpha=0.5)


plt.subplot(133)
plt.imshow(realmap2)
plt.title('realmap 256')
src.plot(color='red',lw=1.5, alpha=0.5)
#ker1.plot(color='red',lw=1.5, alpha=0.5)
#ker2.plot(color='yellow',lw=1.5, alpha=0.5)
#ker3.plot(color='orange',lw=1.5, alpha=0.5)
plt.legend(loc='lower right')
plt.show()


#x=np.linspace(20,50,30)
#X, Y = np.meshgrid(x, x)

#plt.figure()
#ax = plt.axes(projection='3d')
#ax.plot_wireframe(X, Y, rectmap[20:50,20:50], color='black')
#ax.set_title('surface');

#plt.figure()
#ax = plt.axes(projection='3d')
#ax.plot_wireframe(X, Y, realmap[20:50,20:50], color='black')
#ax.set_title('surface');


#plt.figure()
#ax = plt.axes(projection='3d')
#ax.plot_surface(X, Y, realmap[40:60,40:60], rstride=1, cstride=1,cmap='viridis', edgecolor='none')
#ax.plot_surface(X, Y, rectmap[40:60,40:60], rstride=1, cstride=1,cmap='binary', edgecolor='none')
#ax.plot_wireframe(X, Y, rectmap[20:50,20:50]-realmap[20:50,20:50], color='black')
#ax.set_title('surface');

plt.figure()
plt.subplot(121)
plt.plot(rectmap[37],label="gaussian")
plt.plot(realmap[37],label="real512")
plt.plot(realmap2[37],label="real256")
plt.title('horizontal comparison')
plt.legend(loc='lower right')
plt.subplot(122)
plt.plot(rectmap[:,31],label="gaussian")
plt.plot(realmap[:,31],label="real512")
plt.plot(realmap2[:,31],label="real256")
plt.title('vertical comparison')
plt.legend(loc='lower right')

