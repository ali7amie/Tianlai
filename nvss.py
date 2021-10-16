import srcgennvss
import numpy as np
jandecmax=2.51*0.05
jandecstep=0.5*0.05
jandec=np.arange(0,jandecmax,jandecstep)
#jandec=[1,2]
iteration=5
signoise=[0.001]
for i in range(0,len(jandec)):
    for j in range(0,iteration):
        a=srcgennvss.generate()
        a.set_coef(jandec[i])
        a.set_srcfilename_nvss('central_{}_{}.txt'.format(np.round(jandec[i],3),j))
        a.set_srcfilename('offcentral_{}_{}.txt'.format(np.round(jandec[i],3),j))
        a.set_srcfilename2('xyoffcentral_{}_{}.txt'.format(np.round(jandec[i],3),j))
        a.run()
