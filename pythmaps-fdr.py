import numpy as np
import extract_class
import src_generator
import importlib
import matplotlib.pyplot as plt
importlib.reload(extract_class)
importlib.reload(src_generator)


signoise=[0.001,0.005,0.01]
iteration=4
gmap=src_generator.generate(10)
gmap.generator()


rate=[]
s=0
eachmaps=[]
maps=[]
ndec=np.arange(0,11,1)
rateeach=[]
finalrate=[]
finalmaps=[]

for noisex in range(0,len(signoise)):
    for i in range(0,len(ndec)):
        for j in range(0,iteration):
            gmap=src_generator.generate(10)
            gmap.set_sigma_noise(signoise[noisex])
            gmap.generator()

            eachmaps.append(extract_class.extract())
            eachmaps[j].set_n_sigma(ndec[i])
            eachmaps[j].set_k(0)
            eachmaps[j].set_square(55)
            eachmaps[j].set_query(0)
            
    
            eachmaps[j].rectmap=gmap.noise
            eachmaps[j].detect()
            #maps[i].plot()
            eachmaps[j].false_detection_curve()
        
            #rateeach.append(np.shape(eachmaps[j].transpose_srccoor_center_pixcenter)[0]/eachmaps[j].size**2)
            #rateeach.append(len(eachmaps[j].final_list)*49/eachmaps[j].square**2)
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
    plt.errorbar(ndec,meanfinalrate[i],yerr=stdfinalrate[i],label='FDR sigma={} K'.format(signoise[i]))#blue
    plt.xlabel('n')
    plt.ylabel('False detection rate')
    plt.title('False detection rate as a function of n using only noisy maps')

plt.legend(loc='upper right')
plt.show()


