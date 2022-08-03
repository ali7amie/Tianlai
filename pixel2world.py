import math  
import numpy as np 

def pixel2world(projection_center,map_resolution,coor_center):


    ''' This function convert coordinate expressed in pixels on a rectangular map to world coordinate (ra,dec) on the sky 

    parameters:
    -------------
    projection_center: tuple (float, float) 
                       (longitude[deg],latitude[deg])
                
    map_resolution: int
                    resolution in arcminutes

    coor_center: list of tuples (float,float)
                                (vertical coor, horizontal coor)
                          
       
    Return:
    -------  
    coor_world: list of tuples (float,float)
                               (ra[deg],dec[deg])
  '''  

    #converting projection center coordinate from degrees to radians
    
    rarad_center=projection_center[0]*math.pi/180
    decrad_center=projection_center[1]*math.pi/180
    
    #setting the scale

    scale=60/map_resolution    # resolution in arcminutes

    #convert distances in pixels to distances in radians
  
    vertical_rad = coor_center[:,0]/(scale*180/math.pi)   #coor_center[:,0] = vertical coordinate (y) in pixels     | vertical rad is its value in radians
    horizontal_rad = coor_center[:,1]/(scale*180/math.pi)   #coor_center[:,1] = horizontal coordinate (x) in pixels | horizontal rad is its value in radians
    
    #setting up formulas

    D=np.arctan(np.sqrt((vertical_rad**2+horizontal_rad**2)))

    B=np.arctan2(-horizontal_rad,vertical_rad)  #arctan2

    xx=np.sin(decrad_center)*np.sin(D)*np.cos(B) + np.cos(decrad_center)*np.cos(D)

    yy=np.sin(D)*np.sin(B)

    #world coordinate in radians

    rarad=rarad_center+np.arctan2(yy,xx)
    decrad=np.arcsin(np.sin(decrad_center)*np.cos(D)-np.cos(decrad_center)*np.sin(D)*np.cos(B))

    #world coordinate in degrees

    radeg=rarad*180/math.pi
    decdeg=decrad*180/math.pi

    world_coor=np.transpose((radeg,decdeg))
    return world_coor

   # def convert_radeg_xpixels(self)
   # self.scale=60/self.reso
  #  self.A = np.cos(srcdec)*np.cos(srcra-radeg_center)
   # self.F = (scale*(180/math.pi))/(np.sin(decdeg_center)*sin(srcdec)+A*np.cos(decdeg_center))
   # self.srcy = -F * (np.cos(decdeg_center) * np.sin(srcdec) - A * np.sin(decdeg_center)) 
   # self.srcx = -F * np.cos(srcdec) * np.sin(srcra-radeg_center)

#second method if necessary
def convert2(self):
    self.rarad2=np.arctan2(self.xrad,self.yrad)
    self.thetarad=np.arctan(np.sqrt(self.xrad**2+self.yrad**2))
    self.thetaraddeg=self.thetarad*180/math.pi
    self.decrad2=(np.pi/2)-self.thetarad 
    self.rasrc2=self.rarad2*180/math.pi+180
    self.decsrc2=self.decrad2*180/math.pi
    self.srclist2=np.transpose((self.rasrc2,self.decsrc2))
