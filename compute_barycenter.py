import numpy as np

def compute_barycenter(signal_to_noise_results,all_agglomerations): # signal to noise here is signal to noise [1]
    x_barycenter_list=[]
    y_barycenter_list=[]
    for i in range(0,len(all_agglomerations[:,1])):
        x_barycenter=np.sum(signal_to_noise_results[:,4][all_agglomerations[:,1][i]]*signal_to_noise_results[:,1][all_agglomerations[:,1][i]])/np.sum(signal_to_noise_results[:,1][all_agglomerations[:,1][i]])
        y_barycenter=np.sum(signal_to_noise_results[:,3][all_agglomerations[:,1][i]]*signal_to_noise_results[:,1][all_agglomerations[:,1][i]])/np.sum(signal_to_noise_results[:,1][all_agglomerations[:,1][i]])
        x_barycenter_list.append(x_barycenter)
        y_barycenter_list.append(y_barycenter)

    return np.transpose((np.round(y_barycenter_list,1).astype(int),np.round(x_barycenter_list,1).astype(int) ))
    