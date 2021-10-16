#!/bin/tcsh
mkdir eff
cp srcgen.py a.py text_array_generator_efficiency.py post_process.py src_generator.py seven_setup.d extract_class.py nvss_src_S1000_dec80.txt haslam1300eq.fits empty.ppf eff

cd /homeijclab/hamie/TstRun/neww/allinone/eff

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
#set INMAP = haslam1300eq.fits 
set VISFILE = vis_test_sav.ppf 
set OUTMAPF = map_test_sav.fits 
set SIGMAELL = 450
set INMAP = empty.ppf

$JEXE/map2vis -sharp -setup $SETUPF -pdec $NCPDECPOINT -freq $FREQ -healpix $NSIDE -savbeam $BEAMFILE -lmax $LMAX -nthreads 10,16 -mxprod blas,4 -ngb -in $INMAP -src $SRCFILE -out $VISFILE 
$JEXE/vis2map -freq $FREQ -lmax $LMAX -healpix $NSIDE -psithr 0.02,0.001 -nthreads 8,6 -ngb -prt 2 -sharp  -rdbeam $BEAMFILE -saveainv $AINV -mxprod blas,4 $SETUPF $VISFILE visarr_ $NCPDECPOINT $OUTMAPF


python a.py

foreach NOISE ( 0.001 )
    foreach FLUX (0.0 0.025 0.05 0.075 0.1 0.125 0.15 0.175 0.2 0.225 0.25 0.275 0.3 )
	foreach ITERATION ( 0 1 2 3)
	    set SRCFILE = srctxt_{$FLUX}jansky_iteration{$ITERATION}.txt
	    set VISFILE = vis_ncp_autoon_{$NOISE}_{$FREQ}_{$FLUX}_{$ITERATION}_4dec.ppf
	    set OUTMAPF = map_ncp_autoon_{$NOISE}_{$FREQ}_{$FLUX}_{$ITERATION}_4dec.fits
	    set RECALMFILE=map_ncp_autoon_{$NOISE}_{$FREQ}_{$FLUX}_{$ITERATION}_4dec.fits
	    set FILTMAPFILE=filtmap_ncp_autoon_{$NOISE}_{$FREQ}_{$FLUX}_{$ITERATION}_4dec.fits
	    #$JEXE/map2vis -sharp -setup $SETUPF -pdec $NCPDECPOINT -freq $FREQ -healpix $NSIDE -lmax $LMAX -nthreads 10,16 -mxprod blas,4 -ngb -in $INMAP -src $SRCFILE -out $VISFILE 
	    #$JEXE/p4src2vis -setup $SETUPF -src $SRCFILE -freq $FREQ -noise $NOISE -ram 1 -healpix $NSIDE -pdec $NCPDECPOINT -ngb -out $VISFILE 
	    $JEXE/vis2map -freq $FREQ -lmax $LMAX -noise $NOISE -rdbeam $BEAMFILE -rdainv $AINV -healpix $NSIDE -psithr 0.02,0.001 -nthreads 8,6 -ngb -prt 2 -sharp -mxprod blas,4 $SETUPF $VISFILE visarr_ $NCPDECPOINT $OUTMAPF  
	    $JEXE/filt2map -lmax $LMAX -healpix $NSIDE -gaussfilt $SIGMAELL -fmlow 1 $RECALMFILE $FILTMAPFILE 
  	end
    end	
end

setenv JSKYBASE /homeijclab/hamie/JSkyMap 
set JEXE=$JSKYBASE/Objs
set NCPDECPOINT=90,88,86,84
set FREQ=1300
set SETUPF = seven_setup.d
set NSIDE = 256
set LMAX = 700
set SRCFILE = src.txt
set VISFILE = vis_one.ppf
set OUTMAPF = map_one.fits
set RECALMFILE = map_one.fits
set FILTMAPFILE = filtmap_one.fits
set NOISE = 0.001

$JEXE/p4src2vis -setup $SETUPF -src $SRCFILE  -freq $FREQ -ram 1 -healpix $NSIDE -pdec $NCPDECPOINT -ngb -out $VISFILE
$JEXE/vis2map -freq $FREQ -lmax $LMAX  -healpix $NSIDE -nthreads 8,6 -ngb -prt 2 -sharp -mxprod blas,4 $SETUPF $VISFILE visarr_ $NCPDECPOINT $OUTMAPF 
$JEXE/filt2map -lmax $LMAX -healpix $NSIDE  -fmlow 1 $RECALMFILE $FILTMAPFILE 
#$JEXE/vis2map -freq $FREQ -lmax $LMAX -healpix $NSIDE -nthreads 8,6 -ngb -prt 2 -sharp  -mxprod blas,4 $SETUPF $VISFILE visarr_ $NCPDECPOINT $OUTMAPF 

rm vis_one.ppf map_one.fits filtmap_one.fits
foreach NOISE ( 0.001 )
    foreach FLUX ( 0 0.25 0.5 0.75 1 1.25 1.5 1.75 2 2.25 2.5 2.75 3 3.25 3.5 3.75 4 4.25 4.5 4.75 5 5.25 5.5  5.75 6 6.25 6.5 6.75 7 7.25 7.5 7.75 8 8.25 8.5 8.75 9 9.25 9.5 9.75 )
	foreach ITERATION ( 0 1  )
	    set FILTMAPFILE=filtmap_ncp_autoon_{$NOISE}_{$FREQ}_{$FLUX}_{$ITERATION}_4dec.fits
	    if ($FILTMAPFILE=filtmap_ncp_autoon_{$NOISE}_{$FREQ}_{$FLUX}_{$ITERATION}_4dec.fits = FILTMAPFILE=filtmap_ncp_autoon_{$NOISE}_{$FREQ}_0_{$ITERATION}_4dec.fits)
	        mv $FILTMAPFILE filtmap_ncp_autoon_{$NOISE}_{$FREQ}_0.0_{$ITERATION}_4dec.fits
	    endif
  	end
    end	
end
#############################################################

#tes1
#set OUTMAPF = map_test1_sav.fits 
#$JEXE/vis2map -freq $FREQ -lmax $LMAX -healpix $NSIDE -nthreads 8,6 -ngb -prt 2 -sharp  -rdbeam $BEAMFILE -rdainv $AINV -mxprod blas,4 $SETUPF $VISFILE visarr_ $NCPDECPOINT $OUTMAPF

#test2
set VISFILE = vis_test2_sav.ppf 
set OUTMAPF = map_test2_sav.fits 
$JEXE/p4src2vis -setup $SETUPF -freq $FREQ -ram 1 -healpix $NSIDE  -pdec $NCPDECPOINT -ngb -src $SRCFILE -out $VISFILE 
#$JEXE/vis2map -freq $FREQ -lmax $LMAX -healpix $NSIDE -nthreads 8,6 -ngb -prt 2 -sharp  -rdbeam $BEAMFILE -rdainv $AINV -mxprod blas,4 $SETUPF $VISFILE visarr_ $NCPDECPOINT $OUTMAPF
$JEXE/vis2map -freq $FREQ -lmax $LMAX -healpix $NSIDE -nthreads 8,6 -ngb -prt 2 -sharp  -rdainv $AINV -mxprod blas,4 $SETUPF $VISFILE visarr_ $NCPDECPOINT $OUTMAPF



#create maps

python text_array_generator_efficiency.py

foreach NOISE ( 0.001 0.005 )
    foreach FLUX ( 1 11  )
	foreach ITERATION ( 0 1)
	    set SRCFILE = srctxt_{$FLUX}jansky_iteration{$ITERATION}.txt
	    set VISFILE = vis_ncp_autoon_{$NOISE}_{$FREQ}_{$FLUX}_{$ITERATION}_4dec.ppf
	    set OUTMAPF = map_ncp_autoon_{$NOISE}_{$FREQ}_{$FLUX}_{$ITERATION}_4dec.fits
	    set RECALMFILE=map_ncp_autoon_{$NOISE}_{$FREQ}_{$FLUX}_{$ITERATION}_4dec.fits
	    set FILTMAPFILE=filtmap_ncp_autoon_{$NOISE}_{$FREQ}_{$FLUX}_{$ITERATION}_4dec.fits
	    $JEXE/p4src2vis -setup $SETUPF -src $SRCFILE  -freq $FREQ -ram 1 -healpix $NSIDE -rdbeam $BEAMFILE -pdec $NCPDECPOINT -ngb -out $VISFILE 
	    $JEXE/vis2map -freq $FREQ -lmax $LMAX -noise $NOISE -healpix $NSIDE -nthreads 8,6 -ngb -prt 2 -sharp  -rdbeam $BEAMFILE -rdainv $AINV -mxprod blas,4 $SETUPF $VISFILE visarr_ $NCPDECPOINT $OUTMAPF 
	    $JEXE/filt2map -lmax $LMAX -healpix $NSIDE  -fmlow 1 $RECALMFILE $FILTMAPFILE 
  	end
    end	
end

#no saved beam
python text_array_generator_efficiency.py

foreach NOISE ( 0.001 0.005 )
    foreach FLUX ( 1 11  )
	foreach ITERATION ( 0 1)
	    set SRCFILE = srctxt_{$FLUX}jansky_iteration{$ITERATION}.txt
	    set VISFILE = vis_ncp_autoon_{$NOISE}_{$FREQ}_{$FLUX}_{$ITERATION}_4dec.ppf
	    set OUTMAPF = map_ncp_autoon_{$NOISE}_{$FREQ}_{$FLUX}_{$ITERATION}_4dec.fits
	    set RECALMFILE=map_ncp_autoon_{$NOISE}_{$FREQ}_{$FLUX}_{$ITERATION}_4dec.fits
	    set FILTMAPFILE=filtmap_ncp_autoon_{$NOISE}_{$FREQ}_{$FLUX}_{$ITERATION}_4dec.fits
	    $JEXE/p4src2vis -setup $SETUPF -src $SRCFILE  -freq $FREQ -ram 1 -healpix $NSIDE -pdec $NCPDECPOINT -ngb -out $VISFILE 
	    $JEXE/vis2map -freq $FREQ -lmax $LMAX -noise $NOISE -healpix $NSIDE -nthreads 8,6 -ngb -prt 2 -sharp  -mxprod blas,4 $SETUPF $VISFILE visarr_ $NCPDECPOINT $OUTMAPF 
	    $JEXE/filt2map -lmax $LMAX -healpix $NSIDE  -fmlow 1 $RECALMFILE $FILTMAPFILE 
  	end
    end	
end

#testing loops without beam
python text_array_generator_efficiency.py
foreach NOISE ( 0.001 0.005 )
    foreach FLUX ( 1 11 )
	foreach ITERATION ( 0 1)
	    set SRCFILE = srctxt_{$FLUX}jansky_iteration{$ITERATION}.txt
	    set VISFILE = vis_ncp_autoon_{$NOISE}_{$FREQ}_{$FLUX}_{$ITERATION}_4dec.ppf
	    set OUTMAPF = map_ncp_autoon_{$NOISE}_{$FREQ}_{$FLUX}_{$ITERATION}_4dec.fits
	    set RECALMFILE=map_ncp_autoon_{$NOISE}_{$FREQ}_{$FLUX}_{$ITERATION}_4dec.fits
	    set FILTMAPFILE=filtmap_ncp_autoon_{$NOISE}_{$FREQ}_{$FLUX}_{$ITERATION}_4dec.fits
	    $JEXE/p4src2vis -setup $SETUPF -src $SRCFILE  -freq $FREQ -ram 1 -healpix $NSIDE -pdec $NCPDECPOINT -ngb -out $VISFILE 
	    $JEXE/vis2map -freq $FREQ -lmax $LMAX -noise $NOISE -healpix $NSIDE -nthreads 8,6 -ngb -prt 2 -sharp  -mxprod blas,4 $SETUPF $VISFILE visarr_ $NCPDECPOINT $OUTMAPF 
	    $JEXE/filt2map -lmax $LMAX -healpix $NSIDE  -fmlow 1 $RECALMFILE $FILTMAPFILE 
	end
    end
end

#testing loop with beam
python text_array_generator_efficiency.py
set NOISE = 0.001
set FLUX = 11

foreach ITERATION ( 0 1)
   set SRCFILE = srctxt_{$FLUX}jansky_iteration{$ITERATION}.txt
   set VISFILE = vis_ncp_autoon_{$NOISE}_{$FREQ}_{$FLUX}_{$ITERATION}_4dec.ppf
   set OUTMAPF = map_ncp_autoon_{$NOISE}_{$FREQ}_{$FLUX}_{$ITERATION}_4dec.fits
   set RECALMFILE=map_ncp_autoon_{$NOISE}_{$FREQ}_{$FLUX}_{$ITERATION}_4dec.fits
   set FILTMAPFILE=filtmap_ncp_autoon_{$NOISE}_{$FREQ}_{$FLUX}_{$ITERATION}_4dec.fits
   $JEXE/p4src2vis -setup $SETUPF -src $SRCFILE  -freq $FREQ -ram 1 -healpix $NSIDE -rdbeam $BEAMFILE -pdec $NCPDECPOINT -ngb -out $VISFILE 
   $JEXE/vis2map -freq $FREQ -lmax $LMAX -noise $NOISE -healpix $NSIDE -nthreads 8,6 -ngb -prt 2 -sharp  -rdbeam $BEAMFILE -rdainv $AINV -mxprod blas,4 $SETUPF $VISFILE visarr_ $NCPDECPOINT $OUTMAPF 
   $JEXE/filt2map -lmax $LMAX -healpix $NSIDE  -fmlow 1 $RECALMFILE $FILTMAPFILE 
end


#######################
set FREQ = 1300
foreach NOISE ( 0.001 0.005 0.01 )
    foreach FLUX ( 0.0 0.025 0.05 0.075 0.1 0.125 0.15 0.175 0.2 0.225 0.25 0.275 0.3  )
	foreach ITERATION ( 0 1 2 3 4 5 )
	    set SRCFILE2 = central_{$FLUX}_{$ITERATION}.txt
	    set SRCFILE = srctxt_{$FLUX}jansky_iteration{$ITERATION}.txt
	    set VISFILE = vis_ncp_autoon_{$NOISE}_{$FREQ}_{$FLUX}_{$ITERATION}_4dec.ppf
	    set OUTMAPF = map_ncp_autoon_{$NOISE}_{$FREQ}_{$FLUX}_{$ITERATION}_4dec.fits
	    set RECALMFILE=map_ncp_autoon_{$NOISE}_{$FREQ}_{$FLUX}_{$ITERATION}_4dec.fits
	    set FILTMAPFILE=filtmap_ncp_autoon_{$NOISE}_{$FREQ}_{$FLUX}_{$ITERATION}_4dec.fits
	    set SRCFILE3 = xyoffcentral_{$FLUX}_{$ITERATION}.txt
	    set SRCFILE4 = offcentral_{$FLUX}_{$ITERATION}.txt
	    rm $SRCFILE3 $SRCFILE4 $SRCFILE2 $SRCFILE $VISFILE $OUTMAPF $RECALMFILE $FILTMAPFILE
  	end
    end	
end



# test if it take the same beam file
set BEAMFILE = beam_autoon_{$FREQ}_4dec_double.ppf
set AINV = inv_autoon_{$FREQ}_4dec_double.ppf

python text_array_generator_efficiency.py

foreach NOISE ( 0.001 0.005 )
    foreach FLUX ( 1 11  )
	foreach ITERATION ( 0 1)
	    $AINV
	    $BEAMFILE
  	end
    end	
end
