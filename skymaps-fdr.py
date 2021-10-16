import numpy as np
import extract_class
import importlib
import matplotlib.pyplot as plt
import math


# put the output of JSkyMap in an array of fits file name
iteration=2
FREQ=1300
#ndec=np.arange(0,8,3)
ndec=[2]
#signoise=[0.001,0.005, 0.01]
signoise=[0.001]
rate=[]
eachmaps=[]
maps=[]
rateeach=[]
finalrate=[]
finalmaps=[]

for noisex in range(0,len(signoise)):
    for nx in range(0,len(ndec)):
        for j in range(0,iteration):

            eachmaps.append(extract_class.extract())
            eachmaps[j].set_filename('filtmap_ncp_autoon_onlynoise_{}_{}_{}_4dec.fits'.format(signoise[noisex],FREQ,j))
            
            eachmaps[j].set_n_sigma(ndec[nx])
            
            eachmaps[j].extract()
            eachmaps[j].detect()
            eachmaps[j].false_detection_curve()
        
            #rateeach.append(np.shape(eachmaps[j].transpose_srccoor_center_pixcenter)[0]/eachmaps[j].size**2)
            rateeach.append(len(eachmaps[j].final_list)/eachmaps[j].square**2)
        rate.append(rateeach)
        maps.append(eachmaps)
        #ntab.append(ndec[i])
        rateeach=0
        rateeach=[]
        eachmaps=0
        eachmaps=[]

    finalrate.append(rate)
    finalmaps.append(maps)

    rate=0
    rate=[]
    maps=0
    maps=[]
    
    
    

meanrate=[]
stdrate=[]
meanfinalrate=[]
stdfinalrate=[]
for i in range(0,np.shape(finalrate)[0]):
    for j in range(0,np.shape(finalrate)[1]):
        meanrate.append(np.mean(finalrate[i][j]))
        stdrate.append(np.std(finalrate[i][j]))
    meanfinalrate.append(meanrate)
    stdfinalrate.append(stdrate)
    meanrate=0
    meanrate=[]
    stdrate=0
    stdrate=[]



%matplotlib
plt.figure()
for i in range(0,np.shape(finalrate)[0]):
    plt.errorbar(ndec,meanfinalrate[i],yerr=stdfinalrate[i],label='FDR sigma={}'.format(signoise[i]))#blue
    plt.xlabel('n')
    plt.ylabel('False detection rate')
    plt.title('False detection rate with only noise added to visibilities')

plt.legend(loc='upper right')
plt.show()
