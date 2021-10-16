import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
import astropy.stats
from astropy.stats import sigma_clip
from astropy.stats import sigma_clipped_stats
import scipy
import scipy.ndimage as nd
import numpy.ma as ma
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from photutils.aperture import CircularAperture
from photutils.aperture import RectangularAperture
import math
import scipy.constants

def convert_upper_to_center(srccoor_upper_pixcorner,sizef,lenght):
    srccoor_upper_pixcenter=srccoor_upper_pixcorner+np.full((2,np.shape(srccoor_upper_pixcorner)[1]),0.5)
    srccoor_lower_pixcenter=[srccoor_upper_pixcenter[0],sizef-srccoor_upper_pixcenter[1]]
    srccoor_center_pixcenter=srccoor_lower_pixcenter+np.full(np.shape(srccoor_upper_pixcenter),-sizef/2)
    transpose_srccoor_center_pixcenter=np.transpose((-srccoor_center_pixcenter[1],-srccoor_center_pixcenter[0]))
    return [srccoor_upper_pixcorner,srccoor_upper_pixcenter,srccoor_lower_pixcenter,srccoor_center_pixcenter,transpose_srccoor_center_pixcenter]

def convert_center_to_upper(transpose_srccoor_center_pixcenter,sizef,lenght):
    #transpose_srccoor_center_pixcenter=np.transpose((np.zeros(sizef),np.zeros(sizef))
    srccoor_center_pixcenter=np.zeros((2,lenght))
    srccoor_center_pixcenter[1]=-transpose_srccoor_center_pixcenter[:,0]
    srccoor_center_pixcenter[0]=-transpose_srccoor_center_pixcenter[:,1]
    srccoor_lower_pixcenter=srccoor_center_pixcenter- np.full(lenght,-sizef/2)
    srccoor_upper_pixcenter=np.zeros((2,lenght))
    srccoor_upper_pixcenter[0]=srccoor_lower_pixcenter[0]
    srccoor_upper_pixcenter[1]=sizef-srccoor_lower_pixcenter[1]
    transpose_srccoor_upper_pixcenter=np.transpose((srccoor_upper_pixcenter[0],srccoor_upper_pixcenter[1]))
    srccoor_upper_pixcorner=srccoor_upper_pixcenter - np.full((2,lenght),0.5)
    transpose_srccoor_upper_pixcorner=np.transpose((srccoor_upper_pixcorner[0],srccoor_upper_pixcorner[1]))
    transpose_srccoor_upper_pixcenter=np.int_(transpose_srccoor_upper_pixcenter)
    return transpose_srccoor_upper_pixcenter
    
def create_circular_mask(h, w, center=None, radius=None):

    if center is None: # use the middle of the image
        center = (int(w/2), int(h/2))
    if radius is None: # use the smallest distance between the center and image walls
        radius = min(center[0], center[1], w-center[0], h-center[1])

    Y, X = np.ogrid[:h, :w]
    dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2)

    mask = dist_from_center <= radius
    return mask
    

def pad_with(vector, pad_width, iaxis, kwargs):
    pad_value = kwargs.get('padder', 1)
    vector[:pad_width[0]] = pad_value
    vector[-pad_width[1]:] = pad_value
def merge(lsts):
    sts = [set(l) for l in lsts]
    i = 0
    while i < len(sts):
        j = i+1
        while j < len(sts):
            if len(sts[i].intersection(sts[j])) > 0:
                sts[i] = sts[i].union(sts[j])
                sts.pop(j)
            else: j += 1                        
        i += 1
    lst = [list(s) for s in sts]
    return lst


class extract:
    filename=""
    txtfilename=""
    detuse=""
    size = 70 # pixels 
    reso = 12. # arcmin
    rectmap = np.zeros((size,size))
    clipmap=np.zeros((size,size))
    nker1=3
    nker2=5
    nker3=7
    nker4=3
    nker5=5
    nker6=7
    STD = np.zeros((size,size))
    avg1 = np.zeros((size,size))
    avg2 = np.zeros((size,size))
    avg3 = np.zeros((size,size))
    STD = np.zeros((size,size))
    SUM = np.zeros((size,size))
    n=5
    
    c=3*(10**8)
    
    
    k=0
    query_disk_radius=30
    square=55
    decdeg_center=90
    radeg_center=180
    query=0
    FREQ=1300*(10**6)
    maxbaseline=22.5
    radius=np.sqrt(2*2.5**2)
    croping_radius=query_disk_radius
    
    
    def __init__(self,filename="",size=70,reso=15) : 
        self.filename=filename 
        self.reso=reso
        self.size=size 
    def set_filename(self,fnam):
        self.filename=fnam

    def set_nker1(self,nker1): 
        self.nker1=nker1
    def set_nker2(self,nker2): 
        self.nker2=nker2
    def set_nker3(self,nker3): 
        self.nker3=nker3
    def set_nker4(self,nker4): 
        self.nker4=nker4
    def set_nker5(self,nker5): 
        self.nker5=nker5
    def set_nker6(self,nker6): 
        self.nker6=nker6
        
    def set_reso(self,reso): 
        self.reso=reso
    
    def set_maxbaseline(self,maxbaseline): 
        self.maxbaseline=maxbaseline

    def set_freq(self,FREQ): 
        self.FREQ=FREQ

    def set_size(self,size):
        self.size=size
        self.image = np.zeros((self.size,self.size))

    def set_n_sigma(self,n):
        self.n=n

    def set_radius(self,radius):
        self.radius=radius
    
    def set_k(self,k): 
        self.k=k

    def set_query(self,query): 
        self.query=query
        
    def set_detuse(self,detuse):
        self.detuse=detuse

    def set_query_disk_radius(self,query_disk_radius):
        self.query_disk_radius=query_disk_radius
    def set_croping_radius(self,croping_radius):
        self.croping_radius=croping_radius

    def set_square(self,square):
        self.square=square
        
    def set_center(self,decdeg_center,radeg_center):
        self.decdeg_center=decdeg_center
        self.radeg_center=radeg_center

    def extract(self):
        print(self.filename)
        fullmap=hp.read_map(self.filename,dtype=None)
        self.rectmap=hp.gnomview(fullmap,rot=[self.radeg_center,self.decdeg_center],reso=self.reso,xsize=self.size,ysize=self.size,return_projected_map=True,no_plot=True)
    
    def get_image(self):
        return(self.rectmap)

    
    def global_bg_clip(self,s=1.5):
        self.clipmap=sigma_clip(self.rectmap,sigma=s)
        self.clipstd2 = np.std(self.clipmap)
        self.clipmed2 = np.median(self.clipmap)

    def second_clipping(self,s=1.5):
        self.clipmap2=sigma_clip(self.clipmap,sigma=s)
        self.clipstd3 = np.std(self.clipmap2)
        #self.clipmed3 = np.median(self.clipmap2)

    def global_bg(self):
        self.clipstd=np.std(self.rectmap)
        self.clipmed=np.median(self.rectmap)

    def ker1map(self):
            b=1/self.nker1**2
            self.ker1=b*np.ones((self.nker1,self.nker1))
            self.avg1=nd.convolve(self.rectmap,self.ker1)
    def ker4map(self):
            b=1/self.nker4**2
            self.ker4=b*np.ones((self.nker4,self.nker4))
            self.avg4=nd.convolve(self.rectmap,self.ker4)

    
    def ker2map(self):
        d=1/(self.nker2**2-self.nker1**2)
        self.ker2=d*np.pad(0*self.ker1,int((self.nker2-self.nker1)/2),pad_with)
        self.avg2=nd.convolve(self.rectmap,self.ker2)
    def ker5map(self):
        d=1/(self.nker5**2-self.nker4**2)
        self.ker5=d*np.pad(0*self.ker4,int((self.nker5-self.nker4)/2),pad_with)
        self.avg5=nd.convolve(self.rectmap,self.ker5)

    def local_bg(self):
        c=1/(self.nker3**2-self.nker2**2)
        #ker1=np.ones((self.nker1,self.nker1))
        #ker2=np.pad(0*ker1,int((self.nker2-self.nker1)/2),pad_with)
        self.ker3=c*np.pad(0*self.ker2,int((self.nker3-self.nker2)/2),pad_with)
        self.avg3=nd.convolve(self.rectmap, self.ker3)
    
    def ker3map(self):
        self.STD = np.zeros((self.size,self.size))
        #ker1=np.ones((self.nker1,self.nker1))
        #ker2=np.pad(0*ker1,int((self.nker2-self.nker1)/2),pad_with)
        #ker33=np.pad(0*self.ker2,int((self.nker3-self.nker2)/2),pad_with)
        self.mask=np.array(self.ker3==0)
        for i in range(int(self.nker3/2),self.rectmap.shape[0]-int(self.nker3/2)):
            for j in range(int(self.nker3/2),self.rectmap.shape[1]-int(self.nker3/2)):
                if (self.nker3/2).is_integer()==False:
                    self.kerr=self.rectmap[i-int(self.nker3/2):i+(int(self.nker3/2)+1),j-int(self.nker3/2):j+(int(self.nker3/2)+1)]
                    
                else:
                    
                    self.kerr=self.rectmap[i-int(self.nker3/2):i+int(self.nker3/2),j-int(self.nker3/2):j+int(self.nker3/2)]
                    
                self.mkerr=ma.masked_array(self.kerr,self.mask)
                self.std=np.std(self.mkerr)
                self.STD[i][j]=self.std/np.sqrt(self.nker3**2-self.nker2**2)
        #self.STD
    def ker6map(self):
        c=1/(self.nker6**2-self.nker5**2)
        #ker1=np.ones((self.nker1,self.nker1))
        #ker2=np.pad(0*ker1,int((self.nker2-self.nker1)/2),pad_with)
        self.ker6=c*np.pad(0*self.ker5,int((self.nker6-self.nker5)/2),pad_with)
        self.localmedian = np.zeros((self.size,self.size))
        #ker1=np.ones((self.nker1,self.nker1))
        #ker2=np.pad(0*ker1,int((self.nker2-self.nker1)/2),pad_with)
        #ker33=np.pad(0*self.ker2,int((self.nker3-self.nker2)/2),pad_with)
        self.mask2=np.array(self.ker6==0)
        for i in range(int(self.nker6/2),self.rectmap.shape[0]-int(self.nker6/2)):
            for j in range(int(self.nker6/2),self.rectmap.shape[1]-int(self.nker6/2)):
                if (self.nker6/2).is_integer()==False:
                    self.kerr2=self.rectmap[i-int(self.nker6/2):i+(int(self.nker6/2)+1),j-int(self.nker6/2):j+(int(self.nker6/2)+1)]
                    
                else:
                    
                    self.kerr2=self.rectmap[i-int(self.nker6/2):i+int(self.nker6/2),j-int(self.nker6/2):j+int(self.nker6/2)]
                    
                self.mkerr2=ma.masked_array(self.kerr2,self.mask2)
                self.lmedian=ma.median(self.mkerr2)
                self.localmedian[i][j]=self.lmedian/np.sqrt(self.nker6**2-self.nker5**2)
    
    def addker2(self):
        #self.SUM=np.zeros((self.size,self.size))
        #for i in range(int(self.nker3/2),self.rectmap.shape[0]-int(self.nker3/2)):
            #for j in range(int(self.nker3/2),self.rectmap.shape[1]-int(self.nker3/2)):
                #if (self.nker2/2).is_integer()==False:
                    #self.sumi=np.sum(self.rectmap[i-int(self.nker2/2):i+(int(self.nker2/2)+1),j-int(self.nker2/2):j+(int(self.nker2/2)+1)])
                    
                #else:
                    
                    #self.sumi=np.sum(self.rectmap[i-int(self.nker2/2):i+int(self.nker2/2),j-int(self.nker2/2):j+int(self.nker2/2)])
                #self.SUM[i][j]=self.sumi
        self.SUM=np.zeros((self.size,self.size))
        for i in range(int(self.nker6/2),self.rectmap.shape[0]-int(self.nker6/2)):
            for j in range(int(self.nker6/2),self.rectmap.shape[1]-int(self.nker6/2)):
                if (self.nker4/2).is_integer()==False:
                    self.sumi=np.sum(self.rectmap[i-int(self.nker4/2):i+(int(self.nker4/2)+1),j-int(self.nker4/2):j+(int(self.nker4/2)+1)])
                    
                else:
                    
                    self.sumi=np.sum(self.rectmap[i-int(self.nker4/2):i+int(self.nker4/2),j-int(self.nker4/2):j+int(self.nker4/2)])
                self.SUM[i][j]=self.sumi
        self.FLUX=self.SUM-((self.nker4**2)*self.localmedian)
       
    def addker3(self):
        self.SUM=np.zeros((self.size,self.size))
        for i in range(int(self.nker3/2),self.rectmap.shape[0]-int(self.nker3/2)):
            for j in range(int(self.nker3/2),self.rectmap.shape[1]-int(self.nker3/2)):
                if (self.nker2/2).is_integer()==False:
                    self.sumi=np.sum(self.rectmap[i-int(self.nker2/2):i+(int(self.nker2/2)+1),j-int(self.nker2/2):j+(int(self.nker2/2)+1)])
                    
                else:
                    
                    self.sumi=np.sum(self.rectmap[i-int(self.nker2/2):i+int(self.nker2/2),j-int(self.nker2/2):j+int(self.nker2/2)])
                self.SUM[i][j]=self.sumi

        self.FLUX=self.SUM-((self.nker2**2)*self.avg3)

    def first_detection(self):
        #self.limit=self.clipmed+self.n*self.clipstd
        #self.limit=self.avg3+self.n*self.clipstd
        self.limit=self.avg3+self.n*self.STD
        self.binary_map=np.where((self.avg1>self.avg3)&(self.avg1>=self.limit)&(self.avg1>self.avg2)&(self.rectmap>0),1,0)
        self.srccoor_upper_pixcorner=np.where((self.avg1>self.avg3)&(self.avg1>=self.limit)&(self.avg1>self.avg2)&(self.rectmap>0))
        self.transpose_srccoor_upper_pixcorner=np.transpose((self.srccoor_upper_pixcorner[0],self.srccoor_upper_pixcorner[1]))
        self.srccoor_upper_pixcenter=self.srccoor_upper_pixcorner+np.full((2,np.shape(self.srccoor_upper_pixcorner)[1]),0.5)
        self.srccoor_lower_pixcenter=[self.srccoor_upper_pixcenter[0],self.size-self.srccoor_upper_pixcenter[1]]
        self.srccoor_center_pixcenter=self.srccoor_lower_pixcenter+np.full(np.shape(self.srccoor_upper_pixcenter),-self.size/2)
        self.transpose_srccoor_center_pixcenter=np.transpose((-self.srccoor_center_pixcenter[1],-self.srccoor_center_pixcenter[0]))

   
    def src_signal_to_noise(self):
        self.snpix=[]
        sn=self.avg1/self.STD
        self.transpose_srccoor_upper_pixcorner=np.transpose((self.srccoor_upper_pixcorner[0],self.srccoor_upper_pixcorner[1]))
        for i in range(0,len(self.transpose_srccoor_center_pixcenter)):
            sn2pix=sn[self.transpose_srccoor_upper_pixcorner[i,0]][self.transpose_srccoor_upper_pixcorner[i,1]]
            self.snpix.append(sn2pix)
            
        

    def table(self):
    
        self.tab2=[np.transpose(self.snpix),self.srccoor_center_pixcenter[0],self.srccoor_center_pixcenter[1]]
        
        self.indexminmaxsn=np.argsort(self.tab2[0])
        self.tab2[0]=np.sort(self.tab2[0])
        self.tab2[1]=self.tab2[1][self.indexminmaxsn]
        self.tab2[2]=self.tab2[2][self.indexminmaxsn]
        self.tab2[0]=np.flip(self.tab2[0])
        self.tab2[1]=np.flip(self.tab2[1])
        self.tab2[2]=np.flip(self.tab2[2])
        self.indexmaxminsn=np.flip(self.indexminmaxsn)
        self.tab=[self.indexmaxminsn,self.tab2[0],self.tab2[1],self.tab2[2]]
        self.transpose_tab=np.transpose(self.tab)


    
    def clustering(self):
        flag=np.zeros(np.shape(self.tab[0])[0])
        self.flag2=np.zeros(np.shape(self.tab[0])[0])
        self.clusters_index_list=[]
        for i in range(0,np.shape(self.tab[0])[0]):
            if flag[self.tab[0][i]]==0:
               ### print('************************************')
               ### print('pixel {} '.format(self.tab[0][i]))
               ### print('flag  {} '.format(flag[self.tab[0][i]]))
                tabb=np.where(flag==0)
                mod=np.sqrt( (self.tab[2][i]-self.tab[2])**2 + (self.tab[3][i]-self.tab[3])**2 )
                rad=np.where(mod<=self.radius)
                radd=self.tab[0][rad]
                raddd=np.isin(radd,tabb)
                radddd=radd[raddd]
                print('pixels inside the disk  {} '.format(radd))
                print(raddd)
                print('pixels with flag 0 inside the disk  {} '.format(radddd))
                ngtab=0
                ngtab=[]
                for l in range(0,np.shape(radddd)[0]):
                    mod2=np.sqrt( (self.tab[2][np.where(self.tab[0]==radddd[l])]-self.tab[2][np.where(np.isin(self.tab[0],radddd))] )**2+(self.tab[3][np.where(self.tab[0]==radddd[l])]-self.tab[3][np.where(np.isin(self.tab[0],radddd))])**2)
                    ng=np.where(mod2<=np.sqrt(2))
                    ngg=radddd[ng] 
                    print('pixel {}  neighbours {} '.format(radddd[l],ngg))
                    ngtab.append(ngg)
                    ngtab2=merge(ngtab)
                print('clustered pixels {}'.format(ngtab2))


                
                self.clusters_index_list.append(ngtab2) 
        
                for m in range(0,np.shape(ngtab2)[0]):
                    flag[ngtab2[m]]=1
                    self.flag2[ngtab2[m]]=np.shape(self.clusters_index_list)[0]
                for h in range(0,len(ngtab2)):
                    print('flags {}'.format(flag[ngtab2[h]]))
#
                #print(self.clusters_index_list)

    
    def filter_clusters_index_list(self):        
        for i in range(0,len(self.clusters_index_list)):
            if   (len(self.clusters_index_list[i])>1):
                print('this row have many agglomeration')
                print(self.clusters_index_list[i])
            if   (len(self.clusters_index_list[i])>1)  &  (len(self.clusters_index_list[i][0])>3):
#
                
                print('the size of first agglomeration is larger than 3 pixels')
                print(self.clusters_index_list[i])
#
                #for j in range(1,len(self.clusters_index_list[i])):
                for j in range(len(self.clusters_index_list[i])-1,1,-1):
                    print('i={} and j=[]'.format(i,j))
                    if len(self.clusters_index_list[i][j])<3:
#
                        
                        self.clusters_index_list[i].remove(self.clusters_index_list[i][j])

    def correct_clusters_index_list(self):
        for i in range(0,len(self.clusters_index_list)):
           self.clusters_index_list[i]=self.clusters_index_list[i][0]  

    def clusters_dimension(self):
          self.pixnum=[]
          for i in range(0,len(self.clusters_index_list)):
              self.pixnum.append(len(self.clusters_index_list[i]))

    def compute_barycenter(self):
        self.clusters_coor_list=np.ones_like(self.clusters_index_list)
        #self.clusters_coor_list=np.ones_like(self.clusters_index_list[0])
        
        
        for i in range(0,len(self.clusters_index_list)):
            self.clusters_coor_list[i]=self.transpose_srccoor_center_pixcenter[self.clusters_index_list[i]]

        coor2=np.ones_like(self.clusters_index_list)
        for i in range(0,len(self.clusters_index_list)):
            coor2[i]=self.transpose_srccoor_center_pixcenter[self.clusters_index_list[i]]

        self.clusters_intensity_list=coor2      
        for i in range(0,len(coor2)):
            for j in range(0,len(coor2[i])):
                self.clusters_intensity_list[i][j]=self.rectmap[self.transpose_srccoor_upper_pixcorner[self.clusters_index_list[i][j]] [0] ] [self.transpose_srccoor_upper_pixcorner[self.clusters_index_list[i][j]][1] ]

        corint=self.clusters_coor_list*self.clusters_intensity_list


        self.clusters_sumintensity_list=[]
        for i in range(0,len(self.clusters_intensity_list)):
            self.clusters_sumintensity_list.append(np.sum(self.clusters_intensity_list[i])/2)
        sumcorint=0
        sumcorint=[]
        for i in range(0,len(corint)):
            sumcorint.append(np.sum(corint[i],axis=0))


        sumcorint=np.array(sumcorint)
        self.clusters_barycenter_list=[sumcorint[:,0]/self.clusters_sumintensity_list,sumcorint[:,1]/self.clusters_sumintensity_list]
        self.transpose_clusters_barycenter_list=np.transpose((self.clusters_barycenter_list[0],self.clusters_barycenter_list[1]))
        self.fluxdensityK=self.FLUX[np.int_(self.transpose_clusters_barycenter_list)]
        self.fluxdensityK=[]
        self.localbg=[]
        self.barycentercoor_upper_pixcenter=convert_center_to_upper(self.transpose_clusters_barycenter_list,self.size,np.shape(self.transpose_clusters_barycenter_list)[0])
        for i in range(0,np.shape(self.transpose_clusters_barycenter_list)[0]):
            self.fluxdensityK.append(self.FLUX[self.barycentercoor_upper_pixcenter[:,0][i]][self.barycentercoor_upper_pixcenter[:,1][i]])
            if self.detuse=='skymaps':
                self.localbg.append(self.localmedian[self.barycentercoor_upper_pixcenter[:,0][i]][self.barycentercoor_upper_pixcenter[:,1][i]])
            if self.detuse=='pythmaps':
                self.localbg.append(self.avg3[self.barycentercoor_upper_pixcenter[:,0][i]][self.barycentercoor_upper_pixcenter[:,1][i]])
            
            


    def search_nearer_src(self):
        self.nearer_distance_list=[]
        self.nearer_index_list=[]
        for i in range(0,np.shape(self.clusters_barycenter_list)[1]):
            mod3=np.sqrt( (self.clusters_barycenter_list[0][i]-self.clusters_barycenter_list[0])**2 + (self.clusters_barycenter_list[1][i]-self.clusters_barycenter_list[1])**2 )
            mod3[np.where(mod3==0)]=1000
            self.nearer_distance_list.append(min(mod3))
            self.nearer_index_list.append(np.where(mod3==min(mod3))[0][0])
  
    def convert_xpixels_radeg(self):

        self.decrad_center=self.decdeg_center*math.pi/180
        self.rarad_center=self.radeg_center*math.pi/180
        self.scale=60/self.reso
        #self.xrad=self.final_list[:,1]/(self.scale*180/math.pi)
        #self.yrad=self.final_list[:,2]/(self.scale*180/math.pi)
        self.xrad=self.clusters_barycenter_list[0]/(self.scale*180/math.pi)
        self.yrad=self.clusters_barycenter_list[1]/(self.scale*180/math.pi)
        self.D=np.arctan(np.sqrt((self.xrad**2+self.yrad**2)))
        self.B=np.arctan2(-self.xrad,self.yrad)
        self.xx=np.sin(self.decrad_center)*np.sin(self.D)*np.cos(self.B) + np.cos(self.decrad_center)*np.cos(self.D)
        self.yy=np.sin(self.D)*np.sin(self.B)
        
        self.rarad=self.rarad_center+np.arctan2(self.yy,self.xx)
        self.decrad=np.arcsin(np.sin(self.decrad_center)*np.cos(self.D)-np.cos(self.decrad_center)*np.sin(self.D)*np.cos(self.B))
        self.radeg=self.rarad*180/math.pi
        self.decdeg=self.decrad*180/math.pi
        self.srclist=np.transpose((self.radeg,self.decdeg))

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

    

    def convert_kelvin_jansky(self):
        
        self.c=3*(10**8)
        self.wavelenght=self.c/self.FREQ
        
        #self.beam_surf=(self.wavelenght/self.maxbaseline)**2#33(58 arcmin**2 in rad) boltzmann (J.k-1)
        #self.clusters_sumintensity_list=np.array(self.clusters_sumintensity_list)
        self.beam_surf=(self.reso*(1/60)*(math.pi/180))**2
        self.fluxdensityK=np.asarray(self.fluxdensityK)
        #self.Watt=2*scipy.constants.Boltzmann*self.clusters_sumintensity_list*self.beam_surf/self.wavelenght**2
        self.Watt=(2*scipy.constants.Boltzmann*self.fluxdensityK*self.beam_surf)/(self.wavelenght**2)
        self.Jansky=self.Watt/(10**(-26))


  
    def final_result(self):
        
        #self.final_list=np.transpose(( range(0,np.shape(self.clusters_barycenter_list)[1]), self.clusters_barycenter_list[0], self.clusters_barycenter_list[1], self.clusters_sumintensity_list, self.pixnum, self.nearer_distance_list, self.nearer_index_list))
        #self.final_list=np.transpose(( range(0,np.shape(self.clusters_barycenter_list)[1]), self.clusters_barycenter_list[0], self.clusters_barycenter_list[1],self.radeg,self.decdeg, self.clusters_sumintensity_list,self.Jansky, self.pixnum, self.nearer_distance_list, self.nearer_index_list))
        self.final_list=np.transpose(( range(0,np.shape(self.clusters_barycenter_list)[1]), self.clusters_barycenter_list[0], self.clusters_barycenter_list[1],self.radeg,self.decdeg, self.fluxdensityK,self.Jansky, self.pixnum, self.nearer_distance_list, self.nearer_index_list,self.localbg))
        #self.final_list=np.around(self.final_list,4)
        print('*********************all detected sources***************')
        print(self.final_list)
        self.before_filtering_query=[]
        self.before_filtering_query=self.final_list

    
        
    def filter_tiny_src(self):
        self.final_list=self.final_list[np.where(self.final_list[:,7]>self.k)]
        self.final_list=np.around(self.final_list,7)
        print('*********************sources bigger than given size *********')
        print(self.final_list)
        self.before_query=[]
        self.before_query=  self.final_list
        
    
    def query_disk(self):
        self.final_list=self.final_list[np.where(np.sqrt(self.final_list[:,1]**2+self.final_list[:,2]**2)<self.query_disk_radius)] 
        print('***********************sources inside croped disk*************')
        print(self.final_list)
        
    def query_rect(self):
        self.final_list=self.final_list[np.where((abs(self.final_list[:,1])<self.square/2)&(abs(self.final_list[:,2])<self.square/2))] 
        print('***********************sources inside croped square*************')
        print(self.final_list)

    def bg_estimation(self):
        self.cropedmap=create_circular_mask(self.size,self.size,center=None,radius=self.croping_radius)*self.rectmap
        self.mapstd=np.std(self.cropedmap)
        self.clipcropedmap=sigma_clip(self.cropedmap,sigma=3)
        self.clipmapstd = np.std(self.clipcropedmap)
        
    
    def detect(self):
        self.global_bg_clip()
        self.second_clipping()
        self.global_bg()
        self.ker1map()
        self.ker2map()
        self.local_bg()
        self.ker3map()
        self.ker4map()
        self.ker5map()
        self.ker6map()
        if self.detuse=='skymaps':
            self.addker2()
        if self.detuse=='pythmaps':
            self.addker3()
        self.first_detection()
        self.src_signal_to_noise()
        self.table()
        self.clustering()
        self.filter_clusters_index_list()
        self.correct_clusters_index_list()
        self.clusters_dimension()
        self.compute_barycenter()
        self.search_nearer_src()
        self.convert_xpixels_radeg()
        self.convert2()
        self.convert_kelvin_jansky()
        self.final_result()
        self.filter_tiny_src()
        if self.detuse=='skymaps':
            self.query_disk()
        if self.detuse=='pythmaps':
            self.query_rect()
        self.bg_estimation()
        
        
    
    def plot(self):
        rayon=0.25

        plt.figure()
        plt.imshow(self.clipmap,extent=[-self.size/2., self.size/2., -self.size/2., self.size/2. ])
        plt.title('clipped map with med={} and std={}'.format(self.clipmed2,self.clipstd2))
        plt.colorbar()
        plt.show()



        
        plt.figure()
        plt.imshow(self.cropedmap,extent=[-self.size/2., self.size/2., -self.size/2., self.size/2. ])
        plt.title('sigma on map ={}'.format(self.mapstd))
        plt.colorbar()
        plt.show()

        
        plt.figure()
        plt.imshow(self.clipcropedmap,extent=[-self.size/2., self.size/2., -self.size/2., self.size/2. ])
        plt.title('clipped cropedmap sigma={}'.format(self.clipmapstd))
        plt.colorbar()
        plt.show()
        
        

        alldetectedsrc= CircularAperture(self.transpose_srccoor_center_pixcenter, r=0.5*rayon)
        plt.figure()
        plt.imshow(self.rectmap,extent=[-self.size/2., self.size/2., -self.size/2., self.size/2. ])
        alldetectedsrc.plot(color='red', lw=1.5, alpha=0.5)
        plt.title('all detected sources')
        plt.colorbar()
        plt.show()

        srcinradbarycenter =   CircularAperture(self.transpose_clusters_barycenter_list, r=(self.nker2/2)/0.25*rayon)
        plt.figure()
        plt.imshow(self.rectmap,extent=[-self.size/2., self.size/2., -self.size/2., self.size/2. ])
        alldetectedsrc.plot(color='red', lw=1.5, alpha=0.5)
        srcinradbarycenter.plot(color='yellow',lw=1.5, alpha=0.5)
        plt.title('clustering')
        plt.colorbar()
        plt.show()

        baricenters  = CircularAperture(self.transpose_clusters_barycenter_list, r=rayon)
        srcbeforefilterquery=CircularAperture(np.transpose((self.before_filtering_query[:,1],self.before_filtering_query[:,2])),0.5)#same as barycenter
        plt.figure()
        plt.imshow(self.rectmap,extent=[-self.size/2., self.size/2., -self.size/2., self.size/2. ])
        srcinradbarycenter.plot(color='yellow',lw=1.5, alpha=0.5)
        baricenters.plot(color='red', lw=4, alpha=4)
        alldetectedsrc.plot(color='red', lw=1.5, alpha=0.5)
        plt.title('barycenters')
        plt.colorbar()
        plt.show()

        srcbeforequery=CircularAperture(np.transpose((self.before_query[:,1],self.before_query[:,2])),0.5)
        srcinradbeforequery =   CircularAperture(np.transpose((self.before_query[:,1],self.before_query[:,2])), r=(self.nker2/2)/0.25*rayon)
        plt.figure()
        plt.imshow(self.rectmap,extent=[-self.size/2., self.size/2., -self.size/2., self.size/2. ])
        srcinradbeforequery.plot(color='yellow',lw=1.5, alpha=0.5)
        srcbeforequery.plot(color='red',lw=1.5, alpha=0.5)
        plt.title('filtering')
        plt.colorbar()
        plt.show()

        srcinrad =   CircularAperture(np.transpose((self.final_list[:,1],self.final_list[:,2])), r=(self.nker2/2)/0.25*rayon)
        finalsrc = CircularAperture(np.transpose((self.final_list[:,1],self.final_list[:,2])),0.5)
        centerradius= CircularAperture([0,0],self.query_disk_radius )
        centersquare=RectangularAperture([0,0],self.square,self.square)
        plt.figure()
        plt.imshow(self.rectmap,extent=[-self.size/2., self.size/2., -self.size/2., self.size/2. ])
        #centerradius.plot(color='white',lw=1.5, alpha=0.5)
        centersquare.plot(color='red',lw=1.5, alpha=0.5)
        srcinrad.plot(color='yellow',lw=1.5, alpha=0.5)
        finalsrc.plot(color='red',lw=1.5, alpha=0.5)
        #srcinside.plot(color='red',lw=1.5, alpha=0.5)
        plt.title('final results')
        plt.colorbar()
        plt.show()


        
        #inputsrc = CircularAperture(positions, r=rayon)
        
        #srcshape =  CircularAperture(transbarycenter, r=6*rayon)
        #srcinrad =   CircularAperture(self.transpose_clusters_barycenter_list, r=(self.nker2/2)/0.25*rayon)
        
        
        srcinside = RectangularAperture(np.transpose((self.final_list[:,1],self.final_list[:,2])), 1,1)
        finalker1 = RectangularAperture(np.transpose((self.final_list[:,1],self.final_list[:,2])), self.nker1,self.nker1)
        finalker2 = RectangularAperture(np.transpose((self.final_list[:,1],self.final_list[:,2])), self.nker2,self.nker2)
        finalker3 = RectangularAperture(np.transpose((self.final_list[:,1],self.final_list[:,2])), self.nker3,self.nker3)
        
        


        #rectmap
        plt.figure()
        plt.imshow(self.rectmap,extent=[-self.size/2., self.size/2., -self.size/2., self.size/2. ])
        #inputsrc.plot(color='black', lw=6, alpha=0.5)
        #alldetectedsrc.plot(color='blue', lw=1.5, alpha=0.5)
        #baricenters.plot(color='red', lw=1.5, alpha=0.5)
        #srcshape.plot(color='yellow', lw=1.5, alpha=0.5)
        #centerradius.plot(color='white',lw=1.5, alpha=0.5)
        #centersquare.plot(color='yellow',lw=1.5, alpha=0.5)
        centerradius.plot(color='white',lw=1.5, alpha=0.5)
        srcinrad.plot(color='white',lw=1.5, alpha=0.5)
        finalsrc.plot(color='red',lw=1.5, alpha=0.5)
        #srcinside.plot(color='red',lw=1.5, alpha=0.5)

        #finalker1.plot(color='red',lw=1.5, alpha=0.5)
        #finalker2.plot(color='red',lw=1.5, alpha=0.5)
        #finalker3.plot(color='white',lw=1.5, alpha=0.5)
        plt.title('map')
        plt.colorbar()
        plt.show()

        #plt.figure()
        #plt.imshow(self.avg1,extent=[-self.size/2., self.size/2., -self.size/2., self.size/2. ])
        #inputsrc.plot(color='black', lw=6, alpha=0.5)
        #alldetectedsrc.plot(color='blue', lw=1.5, alpha=0.5)
        #centerradius.plot(color='white',lw=1.5, alpha=0.5)
        #baricenters.plot(color='red', lw=1.5, alpha=0.5)
        #srcshape.plot(color='yellow', lw=1.5, alpha=0.5)
        #srcinrad.plot(color='white',lw=1.5, alpha=0.5)
        #plt.title('map')
        #plt.colorbar()
        #plt.show()


        #binary map
        plt.figure()
        plt.imshow(self.binary_map,extent=[-self.size/2., self.size/2., -self.size/2., self.size/2. ])
        #inputsrc.plot(color='black', lw=6, alpha=0.5)
        #alldetectedsrc.plot(color='white', lw=1.5, alpha=0.5)
        centerradius.plot(color='white',lw=1.5, alpha=0.5)
        #centersquare.plot(color='red',lw=1.5, alpha=0.5)
        #baricenters.plot(color='red', lw=1.5, alpha=0.5)
        #srcshape.plot(color='yellow', lw=1.5, alpha=0.5)
        srcinrad.plot(color='white',lw=1.5, alpha=0.5)
        #centersquare.plot(color='yellow',lw=1.5, alpha=0.5)
        finalsrc.plot(color='red',lw=1.5, alpha=0.5)
        
        plt.title('binary map')
        plt.colorbar()
        plt.show()
        
        
        srcpixel = RectangularAperture(np.transpose((self.final_list[:,1],self.final_list[:,2])), 1,1)
        srczone =  CircularAperture(np.transpose((self.final_list[:,1],self.final_list[:,2])), r=self.nker3/2)
    

        #rectmap after filtration
        #plt.figure()
        #plt.imshow(self.rectmap,extent=[-self.size/2., self.size/2., -self.size/2., self.size/2. ])
        #srcpixel.plot(color='red',lw=1.5, alpha=0.5)
        #srczone.plot(color='white',lw=1.5, alpha=0.5)
        #centerradius.plot(color='white',lw=1.5, alpha=0.5)
        #centersquare.plot(color='yellow',lw=1.5, alpha=0.5)
        #plt.title('map')
        #plt.colorbar()
        #plt.show()

        
   

        

    def set_txtfilename(self,txtfnam):
        self.txtfilename=txtfnam

    def export(self):
        txtfile=open(self.txtfilename,'w')
        np.savetxt(txtfile,self.final_list,fmt='%.4f',delimiter='     ',header='detected sources number:{} \n index      x            y           ra         dec          flux(K)     flux(Jansky)    size in pix   dist near src  index near src'.format(np.shape(self.final_list)[0]))
        txtfile.close()

    def false_detection_curve(self):
        print('n={}'.format(self.n))
        print('limit/noise in % ={}'.format(self.limit/self.clipstd))
        print('false detections number ={}'.format(np.shape(self.transpose_srccoor_center_pixcenter)[0]))
        print('false detection rate  ={}'.format(np.shape(self.transpose_srccoor_center_pixcenter)[0]/self.size**2))
       


        
