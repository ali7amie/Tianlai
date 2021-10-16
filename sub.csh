#!/bin/tcsh
mkdir eff
cp nvss.py post_process.py a.py srcgen.py src_generator.py srcgennvss.py seven_setup.d extract_class.py nvss_src_S1000_dec80.txt haslam1300eq.fits src.txt empty.ppf eff
cd /homeijclab/hamie/TstRun/neww/allinone/eff

setenv JSKYBASE /homeijclab/hamie/JSkyMap 
set JEXE=$JSKYBASE/Objs
set NCPDECPOINT=90,88,86,84
set FREQ=1300
set SETUPF = seven_setup.d
set NSIDE = 256
set LMAX = 700
#set INMAP = haslam1300eq.fits 
set INMAP = empty.ppf
set SIGMAELL = 450



#save beam and ainv (OK)

#foreach FREQ (1298 1300 1302)
foreach FREQ ( 1300 )
    set VISFILE = testvis_ncp_autoon_{$FREQ}.ppf
    set OUTMAPF = testmap_ncp_autoon_{$FREQ}.fits
    set BEAMFILE = beam_autoon_{$FREQ}.ppf
    set AINV = inv_autoon_{$FREQ}.ppf
    set SRCFILE = src.txt
    $JEXE/map2vis -sharp -setup $SETUPF -pdec $NCPDECPOINT -freq $FREQ -healpix $NSIDE -savbeam $BEAMFILE -lmax $LMAX -nthreads 10,16 -mxprod blas,4 -ngb -in $INMAP -src $SRCFILE -out $VISFILE 
    $JEXE/vis2map -freq $FREQ -lmax $LMAX -healpix $NSIDE -nthreads 8,6 -ngb -psithr 0.02,0.001 -prt 2 -sharp  -rdbeam $BEAMFILE -saveainv $AINV -mxprod blas,4 $SETUPF $VISFILE visarr_ $NCPDECPOINT $OUTMAPF
end

#create txtfile for central map
#python nvss.py 
python a.py

#create central maps
set FREQ = 1300

foreach NOISE ( 0.001 )
    foreach FLUX (0.0 0.025 0.05 0.075 0.1 0.125 0.15)
	foreach ITERATION ( 0 1 2 3 4 5 6 7 )
	    set BEAMFILE = beam_autoon_{$FREQ}.ppf
	    set AINV = inv_autoon_{$FREQ}.ppf
	    #set SRCFILE = central_{$FLUX}_{$ITERATION}.txt
	    set SRCFILE = srctxt_{$FLUX}jansky_iteration{$ITERATION}.txt
	    set VISFILE = vis_ncp_autoon_{$NOISE}_{$FREQ}_{$FLUX}_{$ITERATION}_4dec.ppf
	    set OUTMAPF = map_ncp_autoon_{$NOISE}_{$FREQ}_{$FLUX}_{$ITERATION}_4dec.fits
	    set RECALMFILE=map_ncp_autoon_{$NOISE}_{$FREQ}_{$FLUX}_{$ITERATION}_4dec.fits
	    set FILTMAPFILE=filtmap_ncp_autoon_{$NOISE}_{$FREQ}_{$FLUX}_{$ITERATION}_4dec.fits
	    $JEXE/map2vis -sharp -setup $SETUPF -pdec $NCPDECPOINT -freq $FREQ -healpix $NSIDE -lmax $LMAX -rdbeam $BEAMFILE -nthreads 10,16 -mxprod blas,4 -ngb -in $INMAP -src $SRCFILE -out $VISFILE 
	    $JEXE/vis2map -freq $FREQ -lmax $LMAX -noise $NOISE -healpix $NSIDE -psithr 0.02,0.001 -nthreads 8,6 -ngb -prt 2 -sharp -rdbeam $BEAMFILE -rdainv $AINV -mxprod blas,4 $SETUPF $VISFILE visarr_ $NCPDECPOINT $OUTMAPF 
	    $JEXE/filt2map -lmax $LMAX -healpix $NSIDE -gaussfilt $SIGMAELL -fmlow 1 $RECALMFILE $FILTMAPFILE 
	end
    end
end

#_____________________________________________________
#!/bin/tcsh
mkdir 3freq3
cp nvss.py post_process.py a.py srcgen.py src_generator.py srcgennvss.py seven_setup.d extract_class.py nvss_src_S1000_dec80.txt haslam1300eq.fits src.txt empty.ppf 3freq3
cd /homeijclab/hamie/TstRun/neww/allinone/3freq3

setenv JSKYBASE /homeijclab/hamie/JSkyMap 
set JEXE=$JSKYBASE/Objs
set NCPDECPOINT=90,88,86,84
set FREQ=1300
set SETUPF = seven_setup.d
set NSIDE = 256
set LMAX = 700
set INMAP = haslam1300eq.fits 
set SIGMAELL = 450



#save beam and ainv (OK)

foreach FREQ (1298 1300 1302)
    set VISFILE = testvis_ncp_autoon_{$FREQ}.ppf
    set OUTMAPF = testmap_ncp_autoon_{$FREQ}.fits
    set BEAMFILE = beam_autoon_{$FREQ}.ppf
    set AINV = inv_autoon_{$FREQ}.ppf
    set SRCFILE = src.txt
    $JEXE/map2vis -sharp -setup $SETUPF -pdec $NCPDECPOINT -freq $FREQ -healpix $NSIDE -savbeam $BEAMFILE -lmax $LMAX -nthreads 10,16 -mxprod blas,4 -ngb -in $INMAP -src $SRCFILE -out $VISFILE 
    $JEXE/vis2map -freq $FREQ -lmax $LMAX -healpix $NSIDE -nthreads 8,6 -ngb -psithr 0.02,0.001 -prt 2 -sharp  -rdbeam $BEAMFILE -saveainv $AINV -mxprod blas,4 $SETUPF $VISFILE visarr_ $NCPDECPOINT $OUTMAPF
end

#create txtfile for central map
python nvss.py 


#create central maps
set FREQ = 1300

foreach NOISE ( 0.001 )
    foreach FLUX (0.0 0.025 0.05 0.075 0.1 0.125 )
	foreach ITERATION ( 0 1 2 3 4 )
	    set BEAMFILE = beam_autoon_{$FREQ}.ppf
	    set AINV = inv_autoon_{$FREQ}.ppf
	    set SRCFILE = central_{$FLUX}_{$ITERATION}.txt
	    #set SRCFILE = srctxt_{$FLUX}jansky_iteration{$ITERATION}.txt
	    set VISFILE = vis_ncp_autoon_{$NOISE}_{$FREQ}_{$FLUX}_{$ITERATION}_4dec.ppf
	    set OUTMAPF = map_ncp_autoon_{$NOISE}_{$FREQ}_{$FLUX}_{$ITERATION}_4dec.fits
	    set RECALMFILE=map_ncp_autoon_{$NOISE}_{$FREQ}_{$FLUX}_{$ITERATION}_4dec.fits
	    set FILTMAPFILE=filtmap_ncp_autoon_{$NOISE}_{$FREQ}_{$FLUX}_{$ITERATION}_4dec.fits
	    $JEXE/map2vis -sharp -setup $SETUPF -pdec $NCPDECPOINT -freq $FREQ -healpix $NSIDE -lmax $LMAX -rdbeam $BEAMFILE -nthreads 10,16 -mxprod blas,4 -ngb -in $INMAP -src $SRCFILE -out $VISFILE 
	    $JEXE/vis2map -freq $FREQ -lmax $LMAX -noise $NOISE -healpix $NSIDE -psithr 0.02,0.001 -nthreads 8,6 -ngb -prt 2 -sharp -rdbeam $BEAMFILE -rdainv $AINV -mxprod blas,4 $SETUPF $VISFILE visarr_ $NCPDECPOINT $OUTMAPF 
	    $JEXE/filt2map -lmax $LMAX -healpix $NSIDE -gaussfilt $SIGMAELL -fmlow 1 $RECALMFILE $FILTMAPFILE 
	end
    end
end

#create side maps

foreach NOISE ( 0.001 )
    foreach FREQ ( 1298 1302)
        set BEAMFILE = beam_autoon_{$FREQ}.ppf
	set AINV = inv_autoon_{$FREQ}.ppf
	set SRCFILE = nvss_src_S1000_dec80.txt
	set VISFILE = vis_ncp_autoon_{$NOISE}_{$FREQ}_4dec.ppf
	set OUTMAPF = map_ncp_autoon_{$NOISE}_{$FREQ}_4dec.fits
        set RECALMFILE=map_ncp_autoon_{$NOISE}_{$FREQ}_4dec.fits
	set FILTMAPFILE=filtmap_ncp_autoon_{$NOISE}_{$FREQ}_4dec.fits
	$JEXE/map2vis -sharp -setup $SETUPF  -pdec $NCPDECPOINT -freq $FREQ -healpix $NSIDE -rdbeam $BEAMFILE -lmax $LMAX -nthreads 10,16 -mxprod blas,4 -ngb -in $INMAP -src $SRCFILE -out $VISFILE 
	$JEXE/vis2map -freq $FREQ -lmax $LMAX -rdbeam $BEAMFILE -rdainv $AINV -noise $NOISE -healpix $NSIDE -nthreads 8,6 -ngb -prt 2 -sharp -mxprod blas,4 $SETUPF $VISFILE visarr_ $NCPDECPOINT $OUTMAPF  
	$JEXE/filt2map -lmax $LMAX -healpix $NSIDE -fmlow 1 $RECALMFILE $FILTMAPFILE 
    end
end

#_______________________________________________________________________________________________
#test beam and ainv (OK)

foreach FREQ (1298 1300 1302)
    set VISFILE = testvis_ncp_autoon2_{$FREQ}.ppf
    set OUTMAPF = testmap_ncp_autoon2_{$FREQ}.fits
    set BEAMFILE = beam_autoon_{$FREQ}.ppf
    set AINV = inv_autoon_{$FREQ}.ppf
    set SRCFILE = src.txt
    $JEXE/p4src2vis -setup $SETUPF -freq $FREQ -ram 1 -healpix $NSIDE  -rdbeam $BEAMFILE -pdec $NCPDECPOINT -ngb -src $SRCFILE -out $VISFILE 
    $JEXE/vis2map -freq $FREQ -lmax $LMAX -healpix $NSIDE -nthreads 8,6 -ngb -prt 2 -sharp -rdainv $AINV -mxprod blas,4 $SETUPF $VISFILE visarr_ $NCPDECPOINT $OUTMAPF
end

#create txtfile for central map
python nvss.py 

#create central maps
set FREQ = 1300

foreach NOISE ( 0.001 )
    foreach FLUX (0.0 0.025 0.05 0.075 0.1 0.125 0.15 0.175 0.2 0.225 0.25 0.275 0.3)
	foreach ITERATION ( 0 )
	    set BEAMFILE = beam_autoon_{$FREQ}.ppf
	    set AINV = inv_autoon_{$FREQ}.ppf
	    set SRCFILE = central_{$NOISE}_{$FLUX}_{$ITERATION}.txt
	    set VISFILE = vis_ncp_autoon_{$NOISE}_{$FREQ}_{$FLUX}_{$ITERATION}_4dec.ppf
	    set OUTMAPF = map_ncp_autoon_{$NOISE}_{$FREQ}_{$FLUX}_{$ITERATION}_4dec.fits
	    set RECALMFILE=map_ncp_autoon_{$NOISE}_{$FREQ}_{$FLUX}_{$ITERATION}_4dec.fits
	    set FILTMAPFILE=filtmap_ncp_autoon_{$NOISE}_{$FREQ}_{$FLUX}_{$ITERATION}_4dec.fits
	    $JEXE/map2vis -sharp -setup $SETUPF -pdec $NCPDECPOINT -freq $FREQ -healpix $NSIDE -lmax $LMAX -nthreads 10,16 -mxprod blas,4 -ngb -in $INMAP -src $SRCFILE -out $VISFILE 
	    $JEXE/vis2map -freq $FREQ -lmax $LMAX -noise $NOISE -healpix $NSIDE -nthreads 8,6 -ngb -prt 2 -sharp -rdainv $AINV -mxprod blas,4 $SETUPF $VISFILE visarr_ $NCPDECPOINT $OUTMAPF 
	    #$JEXE/filt2map -lmax $LMAX -healpix $NSIDE -fmlow 1 $RECALMFILE $FILTMAPFILE 
	end
    end
end



# try the loop (OK)
foreach FREQ (1298 1300 1302)
    set VISFILE = testvis_ncp_autoon_{$FREQ}.ppf
    set OUTMAPF = testmap_ncp_autoon_{$FREQ}.fits
    set BEAMFILE = beam_autoon_{$FREQ}.ppf
    set AINV = inv_autoon_{$FREQ}.ppf
    set SRCFILE = src.txt
    $JEXE/p4src2vis -setup $SETUPF -freq $FREQ -ram 1 -healpix $NSIDE  -pdec $NCPDECPOINT -ngb -src $SRCFILE -out $VISFILE 
    $JEXE/vis2map -freq $FREQ -lmax $LMAX -healpix $NSIDE -nthreads 8,6 -ngb -prt 2 -sharp -mxprod blas,4 $SETUPF $VISFILE visarr_ $NCPDECPOINT $OUTMAPF
end
