import post_process

a=post_process.process()
a.set_use('efficiency_pythmaps')
a.set_noise([0.001,0.005,0.01],[0.001,0.005,0.01])
a.set_sigmag(1.5)
a.noise()
a.set_flux(0,3.1*a.signoiseJden[2],2*a.signoiseJden[0])
#a.jandec=[0.5*a.signoiseJden[0],1.8*a.signoiseJden[0],2.5*a.signoiseJden[0],3.5*a.signoiseJden[0],5*a.signoiseJden[0]]
a.set_iteration(3)
a.flux()
#a.set_n(1,8,2)
a.n=[7]
a.run_efficiency_pythmaps()
a.statistics()
%matplotlib
a.set_nx2(0)
a.plot_all()
a.plot_eachn()

import post_process

a=post_process.process()
a.set_use('efficiency_pythmaps')
a.set_noise([0.001],[0.006])
a.set_sigmag(1.4)
a.noise()
a.set_flux(0,3.01*0.05,0.5*0.05)
a.set_iteration(8)
a.flux()
#a.set_n(1,8,2)
a.n=[5,7]
a.src_txt_name()
a.map_fits_name()
a.run_efficiency_skymaps()
a.statistics()
a.set_nx2(0)
%matplotlib
a.plot_all()
a.plot_eachn()


import post_process

a=post_process.process()
a.set_use('efficiency_pythmaps')
a.set_noise([0.001],[0.006])
a.set_freq_subtraction([1298,1300,1302])
a.set_sigmag(1.4)
a.noise()
#a.set_flux(0,2*0.1,0.25*0.1)
a.set_flux(0,2.51*0.05,0.5*0.05)
#a.jandec=[2]
a.set_iteration(5)
a.flux()
#a.set_n(1,8,4)
a.n=[5,7]
a.src_txt_name_subtraction()
a.map_fits_name_subtraction()
a.run_efficiency_subtraction()
a.statistics()
%matplotlib
a.set_nx2(0)
a.plot_all()
a.plot_eachn()
