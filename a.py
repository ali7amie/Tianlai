import numpy as np
import srcgen

iteration=8
#jandecmax=50*0.015
#jandecstep=5*0.015
#jandec=np.arange(0,jandecmax,jandecstep)
#jandec=np.arange(1,10*12,12)

#jandec=[0,2]
#jandec=[0.135]
#jandec=[0, 0.07575578, 0.15151155, 0.22726733, 0.3030231, 0.37877888, 0.45453465, 0.53029043, 0.6060462,  0.68180198]
gtextseach=[]
gtexts=[]
gtextsnameeach=[]
gtextsname=[]
jandec=np.arange(0,3.01*0.05,0.5*0.05)
xygtextsname=[]
xygtextsnameeach=[]
for i in range(0,len(jandec)):
    for j in range(0,iteration):

        gtextseach.append(srcgen.generate())
        gtextseach[j].set_srcfilename('srctxt_{}jansky_iteration{}.txt'.format(np.round(jandec[i],3),j))
        gtextseach[j].set_srcfilename2('xysrctxt_{}jansky_iteration{}.txt'.format(np.round(jandec[i],3),j))
        gtextsnameeach.append(gtextseach[j].srcfilename)
        xygtextsnameeach.append(gtextseach[j].srcfilenamexy)
        gtextseach[j].set_coef(jandec[i])
        gtextseach[j].set_number(10)
        gtextseach[j].run()
    gtexts.append(gtextseach)
    gtextsname.append(gtextsnameeach)
    xygtextsname.append(xygtextsnameeach)
    
    gtextseach=0
    gtextseach=[]
    gtextsnameeach=0
    gtextsnameeach=[]
    xygtextsnameeach=0
    xygtextsnameeach=[]
    