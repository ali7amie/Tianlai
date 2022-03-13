import numpy as np

def hierarchical_merging(table,aperture): #aperture is a list of (limit, aperture) limit for rectangular, aperture for circular
    all_agglomerations_index=[]
    for i in range(0,len(table[:,0])):
        if table[:,5][i]==0:
            pixels_index_in_each_src=[]
            for j in range(0,len(table[:,0])):
                if ( table[:,3][j] >= table[:,3][i]-aperture[0] ) and ( table[:,3][j] <= table[:,3][i]+aperture[0] ) and  ( table[:,4][j] >= table[:,4][i]-aperture[0] ) and ( table[:,4][j] <= table[:,4][i]+aperture[0] ):
                    offset=np.sqrt( ( table[:,3][i] - table[:,3][j] )**2 + (table[:,4][i] - table[:,4][j])**2 )
                    if offset<=aperture[1]:
                        pixels_index_in_each_src.append(int(table[:,0][j]))
                        table[:,5][j]=1
            all_agglomerations_index.append(pixels_index_in_each_src)
            
    return np.transpose(( np.arange(0,len(all_agglomerations_index)), all_agglomerations_index))
