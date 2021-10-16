setenv JSKYBASE /homeijclab/hamie/JSkyMap 
set JEXE=$JSKYBASE/Objs
set NCPDECPOINT=90,88,86,84
set FREQ=1300
set SETUPF = seven_setup.d
set NSIDE = 256
set LMAX = 700
set SRCFILE = noise.txt
set VISFILE = vis_noise_256_0.001.ppf
set OUTMAPF = map_noise_256_0.001.fits
set RECALMFILE = map_noise_256_0.001.fits
set FILTMAPFILE = filtmap_noise_256_0.001.fits
set NOISE = 0.001
set SIGMAELL = 450 


$JEXE/p4src2vis -setup $SETUPF -src $SRCFILE -freq $FREQ -ram 1 -healpix $NSIDE -pdec $NCPDECPOINT -ngb -out $VISFILE 
$JEXE/vis2map -freq $FREQ -lmax $LMAX -noise $NOISE -healpix $NSIDE  -psithr 0.02,0.001  -nthreads 8,6 -ngb -prt 2 -sharp  -mxprod blas,4 $SETUPF $VISFILE visarr_ $NCPDECPOINT $OUTMAPF 
$JEXE/filt2map -lmax $LMAX -healpix $NSIDE -gaussfilt $SIGMAELL -fmlow 1 $RECALMFILE $FILTMAPFILE 



import healpy as hp
fullmap=hp.read_map('filtmap_2ndec.fits')
rectmap=hp.gnomview(fullmap,rot=[180,90],reso=15,xsize=70,ysize=70,return_projected_map=True)

plt.figure()
plt.imshow(rectmap,vmax=2)



import extract_class
a=extract_class.extract()
a.set_filename('filtmap_noise_src.fits')
a.set_detuse('skymaps')
a.set_size(90)
a.extract()
a.detect()
