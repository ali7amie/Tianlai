%matplotlib
import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
import numpy.ma as ma
from photutils.aperture import CircularAperture
from photutils.aperture import RectangularAperture
from mpl_toolkits.mplot3d import Axes3D






#build map
#_____________parameters___________
srccoef=1 #coeficient of src
sigmag=0.9 # sigma of gaussian src
noisesigma=0.8 # sigma of gaussian noise
meannoise=0 #mean of gaussian noise
rayon=0.25 # radius of aperture around src on the map
radius=2.5# radius of src disk
ecart=8 #when comparing to a known list ecart=abs(coor output - coor input)
xsize=100 # size of map related to interferometer FOV
ysize=100
s=2# thresholding level of sigma clipping
n=5 # n*sigma is the level of thresholding
res=12 # in arcmin resolution of map ,should be inferior to healpix resolution(27 arcmin), and related to instrument resolution which is 58 arcmin
xs=90
ys=90
nker1=3
nker2=5
nker3=7

srccoef=1 #coeficient of src
sigmag=0.9 # sigma of gaussian src
sigmanoise=0.01 # sigma of gaussian noise
meannoise=0 #mean of gaussian noise
size=70 # size of map related to interferometer FOV
srcnum=10
srcsize=7
#____________real source_______________
realfull=hp.read_map('filtmap_src.fits')
#realfull=hp.read_map('fitmap_one.fits')

realmap=hp.gnomview(realfull,rot=[90,90],reso=res,xsize=xs,ysize=xs,return_projected_map=True)


gx, gy = np.meshgrid(np.linspace(-2,2,5), np.linspace(-2,2,5))
mux=0
muy=0
g = srccoef*np.exp(-( (gx-mux)**2 / ( 2.0 * sigmag**2 ) ) ) * np.exp(-( (gy-muy)**2 / ( 2.0 * sigmag**2 ) ) )

rectmap=np.zeros((90,90))
rectmap[58-2:58+3,37-2:37+3]=g


coor=[48,47]
src= CircularAperture(coor,r=0.25 )
ker1=RectangularAperture(coor,nker1,nker1)
ker2=RectangularAperture(coor,nker2,nker2)
ker3=RectangularAperture(coor,nker3,nker3)


plt.figure()
plt.subplot(121)
plt.imshow(realmap)
plt.title('real source')
src.plot(color='red',lw=1.5, alpha=0.5)
ker1.plot(color='red',lw=1.5, alpha=0.5)
ker2.plot(color='yellow',lw=1.5, alpha=0.5)
ker3.plot(color='orange',lw=1.5, alpha=0.5)


plt.subplot(122)
plt.imshow(rectmap)
plt.title('Ideal guassian source')
src.plot(color='red',lw=1.5, alpha=0.5)
ker1.plot(color='red',lw=1.5, alpha=0.5)
ker2.plot(color='yellow',lw=1.5, alpha=0.5)
ker3.plot(color='orange',lw=1.5, alpha=0.5)
plt.show()


x=np.linspace(40,60,20)
X, Y = np.meshgrid(x, x)

plt.figure()
ax = plt.axes(projection='3d')
ax.plot_wireframe(X, Y, rectmap[40:60,40:60], color='black')
ax.set_title('surface');

plt.figure()
ax = plt.axes(projection='3d')
ax.plot_wireframe(X, Y, realmap[40:60,40:60], color='black')
ax.set_title('surface');


plt.figure()
ax = plt.axes(projection='3d')
#ax.plot_surface(X, Y, realmap[40:60,40:60], rstride=1, cstride=1,cmap='viridis', edgecolor='none')
#ax.plot_surface(X, Y, rectmap[40:60,40:60], rstride=1, cstride=1,cmap='binary', edgecolor='none')
ax.plot_wireframe(X, Y, rectmap[40:60,40:60]-realmap[40:60,40:60], color='black')
ax.set_title('surface');

plt.figure()
plt.subplot(121)
plt.plot(rectmap[37,:],label="gaussian")
plt.plot(realmap[58,:],label="real")
plt.title('horizontal comparison')
plt.legend(loc='lower right')
plt.subplot(122)
plt.plot(rectmap[:,58],label="gaussian")
plt.plot(realmap[:,37],label="real")
plt.title('vertical comparison')
plt.legend(loc='lower right')
