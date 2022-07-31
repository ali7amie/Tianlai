import math  
import numpy as np 

def pixel2world(projection_center,map_resolution,barycenter_list):

    scale=60/map_resolution    #resolution in arcmin 

    decrad_center=projection_center[1]*math.pi/180
    rarad_center=projection_center[0]*math.pi/180
       
    xrad=barycenter_list[:,0]/(scale*180/math.pi)
    yrad=barycenter_list[:,1]/(scale*180/math.pi)
    D=np.arctan(np.sqrt((xrad**2+yrad**2)))
    B=np.arctan2(-xrad,yrad)
    xx=np.sin(decrad_center)*np.sin(D)*np.cos(B) + np.cos(decrad_center)*np.cos(D)
    yy=np.sin(D)*np.sin(B)
        
    rarad=rarad_center+np.arctan2(yy,xx)
    decrad=np.arcsin(np.sin(decrad_center)*np.cos(D)-np.cos(decrad_center)*np.sin(D)*np.cos(B))
    radeg=rarad*180/math.pi
    decdeg=decrad*180/math.pi
    srclist=np.transpose((radeg,decdeg))
    return srclist

   # def convert_radeg_xpixels(self)
   # self.scale=60/self.reso
  #  self.A = np.cos(srcdec)*np.cos(srcra-radeg_center)
   # self.F = (scale*(180/math.pi))/(np.sin(decdeg_center)*sin(srcdec)+A*np.cos(decdeg_center))
   # self.srcy = -F * (np.cos(decdeg_center) * np.sin(srcdec) - A * np.sin(decdeg_center)) 
   # self.srcx = -F * np.cos(srcdec) * np.sin(srcra-radeg_center)

def convert2(self):
    self.rarad2=np.arctan2(self.xrad,self.yrad)
    self.thetarad=np.arctan(np.sqrt(self.xrad**2+self.yrad**2))
    self.thetaraddeg=self.thetarad*180/math.pi
    self.decrad2=(np.pi/2)-self.thetarad 
    self.rasrc2=self.rarad2*180/math.pi+180
    self.decsrc2=self.decrad2*180/math.pi
    self.srclist2=np.transpose((self.rasrc2,self.decsrc2))

    
