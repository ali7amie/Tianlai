import numpy as np


def merging(aperture,sorted_raw_detection):
  
    id=sorted_raw_detection.index
    all_agglomerations_index=[]


    for i in id:

        if sorted_raw_detection['Flag'][i]==0:

            #targeting pixels of interest / limiting the searched area into a radius=aperture[0]/ to limit computations 

            c1 = sorted_raw_detection['vertical coor'] >= sorted_raw_detection['vertical coor'][i] - aperture[0] 
            c2 = sorted_raw_detection['vertical coor'] <= sorted_raw_detection['vertical coor'][i] + aperture[0] 
            c3 = sorted_raw_detection['horizontal coor'] >= sorted_raw_detection['horizontal coor'][i] - aperture[0]
            c4 = sorted_raw_detection['horizontal coor'] <= sorted_raw_detection['horizontal coor'][i] + aperture[0]

            targeted_pixels=sorted_raw_detection[(c1) & (c2) & (c3) & (c4)]

            #select pixels inside a radius=aperture[1]

            offset = np.sqrt( ( sorted_raw_detection['vertical coor'][i] - targeted_pixels['vertical coor'] )**2 + ( sorted_raw_detection['horizontal coor'][i] - targeted_pixels['horizontal coor'] )**2 )
            selected_pixels=targeted_pixels[offset<=aperture[1]].index
            all_agglomerations_index.append(selected_pixels)

            #assign 1 as a flag for selected pixels
            sorted_raw_detection.loc[targeted_pixels[offset<=aperture[1]].index]=sorted_raw_detection.loc[targeted_pixels[offset<=aperture[1]].index].assign(Flag=1)    

    return all_agglomerations_index



#def merging(aperture,sorted_raw_detection):


 #   ''' This function take an aperture, and the dataframe of all detections'pixels and give a list where each element represent itself a list of indexes contributing to a single sources

  #  parameters:
   # -------------
    #aperture: tuple (float,float)
     #          (3,2.5)  
    #sorted_raw_detection: dataframe 
     #                     all pixels verifying detection criteria sorted by their decreasing S/N ratio
       
    #Return:
    #-------  
    #all_agglomerations_index: list
     #                         a list where each element represent itself a list of indexes contributing to a single sources
    #'''
        

    #id=sorted_raw_detection.index
    #all_agglomerations_index=[]
    #for i in id:
     #   if sorted_raw_detection['Flag'][i]==0:
     #       pixels_index_in_each_src=[]
     #       for j in id:
     #           c1=sorted_raw_detection['vertical coor'][j] >= sorted_raw_detection['vertical coor'][i] - aperture[0] 
     #           c2=sorted_raw_detection['vertical coor'][j] <= sorted_raw_detection['vertical coor'][i] + aperture[0] 
     #           c3=sorted_raw_detection['horizontal coor'][j] >= sorted_raw_detection['horizontal coor'][i] - aperture[0]
     #           c4=sorted_raw_detection['horizontal coor'][j] <= sorted_raw_detection['horizontal coor'][i] + aperture[0]
     #           if (c1 and c2 and c3 and c4):
     #               offset = np.sqrt( ( sorted_raw_detection['vertical coor'][i] - sorted_raw_detection['vertical coor'][j] )**2 + ( sorted_raw_detection['horizontal coor'][i] - sorted_raw_detection['horizontal coor'][j] )**2 )
     #               if offset<=aperture[1]:
     #                   pixels_index_in_each_src.append(j)
     #                   sorted_raw_detection['Flag'][j]=1
     #       all_agglomerations_index.append(pixels_index_in_each_src)
    #return all_agglomerations_index
