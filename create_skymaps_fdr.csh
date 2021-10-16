#!/bin/tcsh
mkdir fdr
cp seven_setup.d extract_class.py skymaps_fdr.py nvss_src_S1000_dec80.txt haslam1300eq.fits fdr
cd /homeijclab/hamie/TstRun/neww/fdr_7dishes/all/fdr

setenv JSKYBASE /homeijclab/hamie/JSkyMap 
set JEXE=$JSKYBASE/Objs
set NCPDECPOINT=90,88,86,84
set FREQ=1300
set SETUPF = seven_setup.d
set NSIDE = 256
set LMAX = 700
set BEAMFILE = beam_autoon_{$FREQ}_4dec_double.ppf
set AINV = inv_autoon_{$FREQ}_4dec_double.ppf
set SRCFILE = nvss_src_S1000_dec80.txt
set INMAP = haslam1300eq.fits 
set VISFILE = vis_test_sav.ppf 
set OUTMAPF = map_test_sav.fits 

#save
$JEXE/map2vis -sharp -setup $SETUPF  -pdec $NCPDECPOINT -freq $FREQ -healpix $NSIDE -lmax $LMAX -nthreads 10,16 -mxprod blas,4 -ngb -in $INMAP -out $VISFILE -savbeam  $BEAMFILE

$JEXE/vis2map -freq $FREQ -lmax $LMAX -healpix $NSIDE -nthreads 8,6 -ngb -prt 2 -sharp  -rdbeam $BEAMFILE -saveainv $AINV -mxprod blas,4 $SETUPF $VISFILE visarr_ $NCPDECPOINT $OUTMAPF


#$JEXE/p4src2vis -setup $SETUPF -freq $FREQ -ram 1 -healpix $NSIDE  -pdec $NCPDECPOINT -ngb -src $SRCFILE -out $VISFILE 
#$JEXE/vis2map -freq $FREQ -lmax $LMAX -healpix $NSIDE -nthreads 8,6 -ngb -prt 2 -sharp  -savbeam $BEAMFILE -saveainv $AINV -mxprod blas,4 $SETUPF $VISFILE visarr_ $NCPDECPOINT $OUTMAPF

foreach NOISE ( 0.001 )
    foreach ITERATION ( 0 1 ) 
        set OUTMAPF = map_ncp_autoon_onlynoise_{$NOISE}_{$FREQ}_{$ITERATION}_4dec.fits
	set RECALMFILE=map_ncp_autoon_onlynoise_{$NOISE}_{$FREQ}_{$ITERATION}_4dec.fits
	set FILTMAPFILE=filtmap_ncp_autoon_onlynoise_{$NOISE}_{$FREQ}_{$ITERATION}_4dec.fits
	$JEXE/vis2map -freq $FREQ -lmax $LMAX  -zerovis -noise $NOISE -healpix $NSIDE -nthreads 8,6 -ngb -prt 2 -sharp -rdbeam $BEAMFILE -rdainv $AINV -mxprod blas,4 $SETUPF $VISFILE visarr_ $NCPDECPOINT $OUTMAPF 
	$JEXE/filt2map -lmax $LMAX -healpix $NSIDE  -fmlow 1 $RECALMFILE $FILTMAPFILE 
     end
end	
