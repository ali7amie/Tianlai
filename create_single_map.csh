export JSKYBASE='/Users/narimanealmoussawi/Documents/new_tianlai/JSKYMAP/JSkyMap'
set JEXE = $JSKYBASE/Objs


set NCPDECPOINT=90,88,86,84
set FREQ=1300
set SETUPF = setup_tianlai.d
set NSIDE = 256
set LMAX = 700
set NOISE = 0.001

set SRCFILE = nvss_src_S1000_dec80.txt
set VISFILE = vis_ncp_{$NOISE}_{$FREQ}_{$NSIDE}_4dec.ppf
set OUTMAPF = map_ncp_{$NOISE}_{$FREQ}_{$NSIDE}_4dec.fits
set RECALMFILE = map_ncp_{$NOISE}_{$FREQ}_{$NSIDE}_4dec.fits
set FILTMAPFILE = filtmap_ncp_{$NOISE}_{$FREQ}_{$NSIDE}_4dec.fits
set SIGMAELL = 450


$JEXE/map2vis -sharp -setup $SETUPF -pdec $NCPDECPOINT -freq $FREQ -healpix $NSIDE -lmax $LMAX -nthreads 10,16 -mxprod blas,4 -ngb -in $INMAP -src $SRCFILE -out $VISFILE 
$JEXE/p4src2vis -setup $SETUPF -src $SRCFILE -freq $FREQ -noise $NOISE -ram 1 -healpix $NSIDE -pdec $NCPDECPOINT -ngb -out $VISFILE 
$JEXE/vis2map -freq $FREQ -lmax $LMAX -noise $NOISE -rdbeam $BEAMFILE -rdainv $AINV -healpix $NSIDE -psithr 0.02,0.001 -nthreads 8,6 -ngb -prt 2 -sharp -mxprod blas,4 $SETUPF $VISFILE visarr_ $NCPDECPOINT $OUTMAPF  
$JEXE/filt2map -lmax $LMAX -healpix $NSIDE -gaussfilt $SIGMAELL -fmlow 1 $RECALMFILE $FILTMAPFILE 
  
