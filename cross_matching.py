#import numpy as np


#def angular_dist(ra1,dec1,ra2,dec2):
 #   r1=np.radians(ra1)
  #  r2=np.radians(ra2)
   # d1=np.radians(dec1)
   # d2=np.radians(dec2)
  
   # a=(np.sin(np.abs(d1-d2)/2))**2
   # b = np.cos(d1)*np.cos(d2)*np.sin(np.abs(r1 - r2)/2)**2
   # d = 2*np.arcsin(np.sqrt(a + b))
  
   # return np.degrees(d)

#def cross_matching(detection_catalog,simulated_catalog,matching_aperture,max_distance):
 #   matches=[]
  #  no_matches=[]
   # for i in range(0,len(detection_catalog)):
    #    for j in range(0,len(simulated_catalog)):
     #       if ( detection_catalog[i][1] >= (simulated_catalog[j][1]-matching_aperture) ) and (detection_catalog[i][1] <= (simulated_catalog[j][1]+matching_aperture) ):
      #          offset=angular_dist(detection_catalog[i][1],detection_catalog[i][2],simulated_catalog[j][1],simulated_catalog[j][2])
       #         if (offset <= max_distance) and (detection_catalog[i][3]>=simulated_catalog[j][3]*0.7) and (detection_catalog[i][3]<= 1.3*simulated_catalog[j][3]):
        #            matches.append((detection_catalog[i][0],simulated_catalog[j][0],offset))
         #       else:
          #          no_matches.append(detection_catalog[i][0])

 #   return (np.array(matches),np.array(no_matches))



import numpy as np
import pandas as pd

def angular_dist(ra1,dec1,ra2,dec2):


    ''' This function take world coordinate of two sources and compute the angular distance between them 

    parameters:
    -------------
    ra1: float
        right ascension of the first source [deg]

    dec1: float
        declination of the second source [deg]
               
    ra2: float
        right ascension of the second source [deg]

    dec2: float
        declination of the second source [deg]

    Return:
    -------   
    angular_distance: float [deg]                
    '''  

    #convert inputs into radians
    r1=np.radians(ra1)
    r2=np.radians(ra2)
    d1=np.radians(dec1)
    d2=np.radians(dec2)
  
    #compute angular distance
    a=(np.sin(np.abs(d1-d2)/2))**2
    b = np.cos(d1)*np.cos(d2)*np.sin(np.abs(r1 - r2)/2)**2
    d = 2*np.arcsin(np.sqrt(a + b))

    #convert it to degrees
    angular_distance = np.degrees(d)
  
    return angular_distance

def cross_matching(detection_catalog,simulated_catalog,matching_aperture,max_distance):


    ''' This function two sources'catalogs and give a list of matched and no matched sources
        To limit the number of operation, the search for matched sources for each source is performed within the given aperture
        and later on if the separation distance between two sources is smaller than the threshold max distance, we consider a matching case
        Note that catalogs should be sorted by their declination

    parameters:
    -------------
    detection_catalog: dataframe 
                       (ra, dec, flux)
        
    simulated_catalog: dataframe
                       (ra, dec, flux)
               
    matching aperture: float
                       
    max_distance: float

    Return:
    -------   
    matching_results: tuple (dataframe, dataframe)
                       (dataframe of mathced sources, dataframe of no matched sources)
                                   
    '''  

    matches=[]
    no_matches=[]

    idd=detection_catalog.index
    ids=simulated_catalog.index
    for i in idd:
        for j in ids:

            #reduce the search area to a disk of radius matching_aperture around each source

            c1 = detection_catalog['dec[deg]'][i] >= (simulated_catalog['dec[deg]'][j] - matching_aperture) 
            c2 = detection_catalog['dec[deg]'][i] <= (simulated_catalog['dec[deg]'][j] + matching_aperture) 

            if c1 and c2:
                offset = angular_dist(detection_catalog['ra[deg]'][i], detection_catalog['dec[deg]'][i], simulated_catalog['ra[deg]'][j], simulated_catalog['dec[deg]'][j])
            
                #within this disk compute angular distance, and if it is smaller than max_distance, we consider a matching case

                c3 = offset <= max_distance
                c4 = detection_catalog['flux[Jy]'][i] >= 0.7 * simulated_catalog['flux[Jy]'][j]
                c5 = detection_catalog['flux[Jy]'][i] <= 1.3 * simulated_catalog['flux[Jy]'][j]
                if c3 and c4 and c5:

                #create a matching results dataframes  
                    matches.append((i, detection_catalog['ra[deg]'][i], detection_catalog['dec[deg]'][i],j , simulated_catalog['ra[deg]'][j],	simulated_catalog['dec[deg]'][j], offset))
                else:
                    no_matches.append((i, detection_catalog['ra[deg]'][i], detection_catalog['dec[deg]'][i]))
    #set up the output
    matches_dataframe = pd.DataFrame(matches,columns=['idx src from det', 'ra1[deg]', 'dec1[deg]','idx src from sim', 'ra2[deg]', 'dec2[deg]', 'offset'])                 
    no_matches_dataframe = pd.DataFrame(no_matches,columns=['idx src from det', 'ra1[deg]', 'dec1[deg]'])

    return (matches_dataframe, no_matches_dataframe)
