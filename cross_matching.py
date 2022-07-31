import numpy as np


def angular_dist(ra1,dec1,ra2,dec2):
    r1=np.radians(ra1)
    r2=np.radians(ra2)
    d1=np.radians(dec1)
    d2=np.radians(dec2)
  
    a=(np.sin(np.abs(d1-d2)/2))**2
    b = np.cos(d1)*np.cos(d2)*np.sin(np.abs(r1 - r2)/2)**2
    d = 2*np.arcsin(np.sqrt(a + b))
  
    return np.degrees(d)

def cross_matching(detection_catalog,simulated_catalog,matching_aperture,max_distance):
    matches=[]
    no_matches=[]
    for i in range(0,len(detection_catalog)):
        for j in range(0,len(simulated_catalog)):
            if ( detection_catalog[i][1] >= (simulated_catalog[j][1]-matching_aperture) ) and (detection_catalog[i][1] <= (simulated_catalog[j][1]+matching_aperture) ):
                offset=angular_dist(detection_catalog[i][1],detection_catalog[i][2],simulated_catalog[j][1],simulated_catalog[j][2])
                if (offset <= max_distance) and (detection_catalog[i][3]>=simulated_catalog[j][3]*0.7) and (detection_catalog[i][3]<= 1.3*simulated_catalog[j][3]):
                    matches.append((detection_catalog[i][0],simulated_catalog[j][0],offset))
                else:
                    no_matches.append(detection_catalog[i][0])

    return (np.array(matches),np.array(no_matches))