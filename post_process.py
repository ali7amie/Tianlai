import numpy as np
import srcgen # class which generate sources text file as input of JSkyMap. 
import src_generator # class to generate synthetic skymaps
import extract_class # the detection algorithm class
import importlib
import matplotlib.pyplot as plt
import scipy.constants
import math
import statistics
import datetime
import matplotlib.dates as mdates
from matplotlib.ticker import AutoMinorLocator
import healpy as hp



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


def convert_Jansky_to_K(S):
    FREQ=1300*10**6
    #FREQ=1420*10**6
    c=3*(10**8)
    wavelenght=c/FREQ
    maxbaseline=16.5
    reso=12
    #beam_surf=6*10**(-5)
    #beam_surf=(wavelenght/maxbaseline)**2#33(58 arcmin**2 in rad) boltzmann (J.k-1)
    beam_surf=(reso*(1/60)*(math.pi/180))**2
    T=( (10**(-26))*S*wavelenght**2)/(2*scipy.constants.Boltzmann*beam_surf)
    return T

        # convert kelvin to Jansky
def convert_K_to_Jansky(T):
    FREQ=1300*10**6
    #FREQ=1420*10**6
    c=3*(10**8)
    wavelenght=c/FREQ
    maxbaseline=16.5
    reso=12
    #beam_surf=6*10**(-5)
    #beam_surf=(wavelenght/maxbaseline)**2#33(58 arcmin**2 in rad) boltzmann (J.k-1)
    beam_surf=(reso*(1/60)*(math.pi/180))**2
    W=(2*scipy.constants.Boltzmann*beam_surf*T)/(wavelenght**2)
    S=W/10**(-26)
    return S

#_____________________________________________________________________________________________



class process:
    delta=0.2
    freq=1300
    FREQ=1300*10**6
    sigmag=1.5
    signoise1=[]
    signoise=[]
    signoiseKden=[]
    signoiseJden=[]
    n=[]
    
    iteration=4
    jandecKden=[]
    jandecKamp=[]
    jandecmKamp=[]
    b=[]
    nx2=1
    use=""
    freq_subtraction=[1298,1300,1302]
    

    def __init__(self,use=""):
        self.use=use

    def set_use(self,use):
        self.use=use

    def set_iteration(self,iteration):
        self.iteration=iteration
    
    def set_freq(self,freq):
        self.freq=freq

    def set_FREQ(self,FREQ):
        self.FREQ=FREQ

    def set_noise(self,signoise1,signoise): 
        self.signoise1=signoise1
        self.signoise=signoise

    def set_sigmag(self,sigmag):
        self.sigmag=sigmag

    def set_nx2(self,nx2):
        self.nx2=nx2

    def set_freq_subtraction(self,freq_subtraction):
        self.freq_subtraction=freq_subtraction

    def noise(self):
        for i in range(0,len(self.signoise)):
            self.signoiseKden.append(2*math.pi*self.signoise[i]*self.sigmag**2) # density in K of a source with S/N=1 ,assuming a sigmag.
        for i in range(0,len(self.signoise)):
            self.signoiseJden.append(convert_K_to_Jansky(self.signoiseKden[i]))# density in Jy

    def set_flux(self,jandecmin,jandecmax,jandecstep): 
        self.jandecmin=jandecmin
        self.jandecmax=jandecmax
        self.jandecstep=jandecstep
        self.jandec=np.arange(self.jandecmin,self.jandecmax,self.jandecstep)
    
    def set_b(self):
        self.b=b   # jandec but with .0

    def flux(self):
        for i in range(0,len(self.jandec)):
            self.jandecKden.append(convert_Jansky_to_K(self.jandec[i]))
        self.jandecKden=np.array(self.jandecKden)
        self.jandecKamp=(1/(2*math.pi*self.sigmag**2))*self.jandecKden
        self.jandecKamp=np.asarray(self.jandecKamp)
        self.jandecmKamp=(1000)*self.jandecKamp

    def set_n(self,nmin,nmax,ndec):
        self.nmin=nmin
        self.nmax=nmax
        self.ndec=ndec
        self.n=np.arange(self.nmin,self.nmax,self.ndec)
    
    
    def src_txt_name(self):
        self.gtextseach=[]
        self.gtexts=[]
        self.gtextsnameeach=[]
        self.gtextsname=[]

        for i in range(0,len(self.jandec)):
            for j in range(0,self.iteration):

                self.gtextseach.append(srcgen.generate())
                self.gtextseach[j].set_srcfilename('srctxt_{}jansky_iteration{}.txt'.format(np.round(self.jandec[i],3),j))
                self.gtextsnameeach.append(self.gtextseach[j].srcfilename)
                self.gtextseach[j].set_coef(self.jandec[i])
        
            self.gtexts.append(self.gtextseach)
            self.gtextsname.append(self.gtextsnameeach)
            self.gtextseach=0
            self.gtextseach=[]
            self.gtextsnameeach=0
            self.gtextsnameeach=[]


    def src_txt_name_subtraction(self):
        self.gtextseach=[]
        self.gtexts=[]
        self.gtextsnameeach=[]
        self.gtextsname=[]

        for i in range(0,len(self.jandec)):
            for j in range(0,self.iteration):
                self.gtextseach.append(srcgen.generate())
                self.gtextseach[j].set_srcfilename('central_{}_{}.txt'.format(np.round(self.jandec[i],3),j))
                self.gtextsnameeach.append(self.gtextseach[j].srcfilename)
                self.gtextseach[j].set_coef(self.jandec[i])
        
            self.gtexts.append(self.gtextseach)
            self.gtextsname.append(self.gtextsnameeach)
            self.gtextseach=0
            self.gtextseach=[]
            self.gtextsnameeach=0
            self.gtextsnameeach=[]

    def map_fits_name(self):
        self.gnameseach=[]
        self.gnames=[]
        self.finalgnames=[]
        for noisex in range(0,len(self.signoise)):
            for i in range(0,len(self.jandec)):
                for j in range(0,self.iteration):
                    self.gnameseach.append('filtmap_ncp_autoon_{}_{}_{}_{}_4dec.fits'.format(self.signoise1[noisex],self.freq,np.round(self.jandec[i],3),j))
                self.gnames.append(self.gnameseach)
                self.gnameseach=0
                self.gnameseach=[]
            self.finalgnames.append(self.gnames)
            self.gnames=0
            self.gnames=[]

    def map_fits_name_subtraction(self):
        self.gnameseach=[]
        self.gnames=[]
        self.finalgnames=[]
        self.finalfinalgnames=[]
        for noisex in range(0,len(self.signoise)):
            for i in range(0,len(self.jandec)):
                for j in range(0,self.iteration):
                    self.gnameseach.append('filtmap_ncp_autoon_{}_{}_{}_{}_4dec.fits'.format(self.signoise1[noisex],self.freq_subtraction[1],np.round(self.jandec[i],3),j))
                self.gnames.append(self.gnameseach)
                self.gnameseach=0
                self.gnameseach=[]
            self.finalgnames.append(self.gnames)
            self.gnames=0
            self.gnames=[]




    def run_efficiency_skymaps(self):
        #maps
        self.deachmaps=[]
        self.dmaps=[]
        self.finaldmaps=[]
        self.finalfinaldmaps=[]

        self.geachmaps=[]
        self.gmaps=[]
        self.finalgmaps=[]
        self.finalfinalgmaps=[]

        #????
        self.detectedsrceach=[]
        self.detectedsrc=[]
        self.finaldetectedsrc=[]
        self.finalfinaldetectedsrc=[]

        # spurious sources, normalized spurious sources
        self.spuriouseach=[]
        self.spurious=[]
        self.finalspurious=[]
        self.finalfinalspurious=[]


        self.normspuriouseach=[]
        self.normspurious=[]
        self.finalnormspurious=[]
        self.finalfinalnormspurious=[]


        self.normspuriouspixeach=[]
        self.normspuriouspix=[]
        self.finalnormspuriouspix=[]
        self.finalfinalnormspuriouspix=[]

        self.normspuriousdeteach=[]
        self.normspuriousdet=[]
        self.finalnormspuriousdet=[]
        self.finalfinalnormspuriousdet=[]

        #all detected coordinate in xy and radec
        self.detectedcooreach=[]
        self.detectedcoor=[]
        self.finaldetectedcoor=[]
        self.finalfinaldetectedcoor=[]

        self.detectedcoordegeach=[]
        self.detectedcoordeg=[]
        self.finaldetectedcoordeg=[]
        self.finalfinaldetectedcoordeg=[]

        #efficiency
        self.efficiencyeach=[]
        self.efficiency=[]
        self.finalefficiency=[]
        self.finalfinalefficiency=[]

        #detected coordinate from the list with flux
        self.truedetectioncooreacheach=[]
        self.truedetectioncooreach=[]
        self.truedetectioncoor=[]
        self.finaltruedetectioncoor=[]
        self.finalfinaltruedetectioncoor=[]

        self.truedetectionfluxeacheach=[]
        self.truedetectionfluxeach=[]
        self.truedetectionflux=[]
        self.finaltruedetectionflux=[]
        self.finalfinaltruedetectionflux=[]



        #?????
        self.truedetectioneach=[]
        self.truedetection=[]
        self.finaltruedetection=[]
        self.finalfinaltruedetection=[]

        #statistics (median,mean) for each map
        self.stateach=[]
        self.stat=[]
        self.finalstat=[]
        self.finalfinalstat=[]

        #flag
        self.finalflag=[]
        self.finalfinalflag=[]
        self.finalfinalfinalflag=[]

        #bilan
        self.bilaneach=[]
        self.bilan=[]
        self.finalbilan=[]
        self.finalfinalbilan=[]

        #???

        self.true=[]
        self.finaltrue=[]
        self.finalfinaltrue=[]

        # not detected src and normalized
        self.notdetectedeach=[]
        self.notdetected=[]
        self.finalnotdetected=[]
        self.finalfinalnotdetected=[]

        self.normnotdetectedeach=[]
        self.normnotdetected=[]
        self.finalnormnotdetected=[]
        self.finalfinalnormnotdetected=[]

        # all detection
        self.alldetection=[]
        self.alldetectioneach=[]
        self.finalalldetection=[]
        self.finalfinalalldetection=[]

        #textfiles
        self.txteach=[]
        self.txt=[]
        self.finaltxt=[]
        self.finalfinaltxt=[]

        self.txtarray=[]
        self.textarray=[]
        self.finaltextarray=[]
        self.finalfinaltextarray=[]

        self.txtnameeach=[]
        self.txtname=[]
        self.finaltxtname=[]
        self.finalfinaltxtname=[]

        self.xytxteach=[]
        self.xytxt=[]
        self.xyfinaltxt=[]
        self.xyfinalfinaltxt=[]

        self.xytxtarray=[]
        self.xytextarray=[]
        self.xyfinaltextarray=[]
        self.xyfinalfinaltextarray=[]


        self.xytxtnameeach=[]
        self.xytxtname=[]
        self.xyfinaltxtname=[]
        self.xyfinalfinaltxtname=[]



        self.detectedfluxeach=[]
        self.finaldif=[]
        self.finalfinaldif=[]
        self.finalfinalfinaldif=[]
        #generated/detected comparison algorithm
        

        for sigx in range(0,len(self.signoise)):
            for nx in range(0,len(self.n)):
                print('for a given false detection rate n={}'.format(self.n[nx]))
                for i in range(0,len(self.jandec)):
                    for j in range(0,self.iteration):
                        self.txtnameeach.append('srctxt_{}jansky_iteration{}.txt'.format(np.round(self.jandec[i],3),j))
                        self.xytxtnameeach.append('xysrctxt_{}jansky_iteration{}.txt'.format(np.round(self.jandec[i],3),j))
                        self.txteach.append(open(self.txtnameeach[j],'r'))
                        self.xytxteach.append(open(self.xytxtnameeach[j],'r'))
                        self.txtarray.append(np.loadtxt(self.txteach[j]))
                        self.xytxtarray.append(np.loadtxt(self.xytxteach[j]))
                        self.txteach[j].close()
                        self.xytxteach[j].close()
                        self.deachmaps.append(extract_class.extract())
                        self.deachmaps[j].set_filename(self.finalgnames[sigx][i][j])
                        self.deachmaps[j].set_size(90)
                        self.deachmaps[j].set_n_sigma(self.n[nx])
                        self.deachmaps[j].set_detuse('skymaps')
                        self.deachmaps[j].set_k(0)
                        self.deachmaps[j].set_query_disk_radius(35)
                        self.deachmaps[j].set_freq(1300*10**6)
                        self.deachmaps[j].set_maxbaseline(16.5)
                        self.deachmaps[j].set_reso(12)
                        self.deachmaps[j].set_nker1(3)
                        self.deachmaps[j].set_nker2(5)
                        self.deachmaps[j].set_radius(np.sqrt(2*(2.5**2)))
                        self.deachmaps[j].set_nker3(7)
                        self.deachmaps[j].set_nker4(4)
                        self.deachmaps[j].set_nker5(8)
                        self.deachmaps[j].set_nker6(10)
                        self.deachmaps[j].extract()
                        self.deachmaps[j].detect()


                        print('START START START START START START START ********************************************')
                #print('srcnumber={} noise= {} n={} flux={} iteration={} '.format(,signoise[sigx],n[nx],jandec[i],j ))
                
                        self.detectedcooreach.append(np.transpose((self.deachmaps[j].final_list[:,1],self.deachmaps[j].final_list[:,2],self.deachmaps[j].final_list[:,5],self.deachmaps[j].final_list[:,6],self.deachmaps[j].final_list[:,10])))
                        self.detectedcoordegeach.append(np.transpose((self.deachmaps[j].final_list[:,3],self.deachmaps[j].final_list[:,4],self.deachmaps[j].final_list[:,5],self.deachmaps[j].final_list[:,6],self.deachmaps[j].final_list[:,10])))
                        self.detectedfluxeach.append(self.deachmaps[j].final_list[:,6])

                        print('generated src')
                        print(self.txtarray[j])
                        print('all detected src')
                        print(self.detectedcooreach[j])
                        self.flag=np.zeros(np.shape(self.detectedcoordegeach[j])[0])

                        print('computing difference between each generated coor and all detected coor')
                        ####
                        self.transposetxtarray=np.transpose((self.txtarray[j][:,0],self.txtarray[j][:,1]))
                        self.xytransposetxtarray=np.transpose((self.xytxtarray[j][:,0],self.xytxtarray[j][:,1]))
                        for o in range(0,np.shape(self.transposetxtarray)[0]):
                            #difeach.append(transposetxtarray[o]-detectedcoordegeach[j][:,[0,1]])
                    
                    
                            self.dif=((180/math.pi)*np.arccos(np.sin(self.transposetxtarray[o][1]*(math.pi/180))*np.sin(self.detectedcoordegeach[j][:,1]*(math.pi/180))+np.cos(self.transposetxtarray[o][1]*(math.pi/180))*np.cos(self.detectedcoordegeach[j][:,1]*(math.pi/180))*np.cos(self.transposetxtarray[o][0]*(math.pi/180)-self.detectedcoordegeach[j][:,0]*(math.pi/180))))
                            self.dif2=self.xytransposetxtarray[o]-self.detectedcooreach[j][:,[0,1]]
                            print('dif')
                            print(self.dif)
                            print(np.where(abs(self.dif)<=self.delta))
                            #matchsrceach.append(detectedcoordegeach[j][np.where(angdist[o]<=delta)])
                            #matchsrcfluxeach.append(detectedfluxeach[j][np.where(angdist[o]<=delta)])
                            #matchsrcdisteach.append(angdist[o][np.where(angdist[o]<=delta)])

                            if len(np.where(abs(self.dif)<=self.delta)[0])!=0:        
                                self.truedetectioncooreacheach.append(self.detectedcooreach[j][np.where( abs(self.dif)<=self.delta )])
                                self.truedetectionfluxeacheach.append(self.detectedfluxeach[j][np.where( abs(self.dif)<=self.delta )])
                                #if len(np.where(abs(dif)<=delta)[0])==1:
                                self.true.append(self.detectedfluxeach[j][np.where( abs(self.dif)<=self.delta )])
                            self.flag[ np.where(abs(self.dif)<=self.delta) ] = self.flag[ np.where(abs(self.dif)<=self.delta) ] +1
                    
                            print('**** **** **** ** **** * * * ** * * ** * * * ** * * * ')
                            print('second flag')
                            print(self.flag)
                            print(self.detectedfluxeach[j][np.where( abs(self.dif)<=self.delta )])
                            print('* * * * ** * * * ** * * * * * * ** * * * * * ** * * * *')

                            s=0
                            ss=0
                            sss=0
                            for k in self.flag:
                                if k==1:
                                    s=s+1
                                if k>1:
                                    ss=ss+k
                                sss=s+ss
                            print('{} one detection/one src'.format(s))
                            print('{} one detection for  many sources'.format(ss))
                            print('total {}'.format(sss))

                        self.truedetectioncooreach.append(self.truedetectioncooreacheach)
                        self.truedetectionfluxeach.append(self.truedetectionfluxeacheach)
                        self.truedetectioncooreacheach=0
                        self.truedetectioncooreacheach=[]
                        self.truedetectionfluxeacheach=0
                        self.truedetectionfluxeacheach=[]
                        self.bilaneach.append([s,ss,sss])
                        self.detectedsrceach.append(sss)        
                        self.spuriouseach.append(len(self.flag)-sss)
                        self.truedetectioncooreacheach2=0
                        self.truedetectioncooreacheach2=[]
                        self.truedetectionfluxeacheach2=0
                        self.truedetectionfluxeacheach2=[]
                        self.notdetectedeach.append(len(self.txtarray[j])-sss)
    
                        self.alldetectioneach.append(len(self.flag))
                        if len(self.flag)!=0:
                            self.normspuriouseach.append((len(self.flag)-sss)/len(self.txtarray[j]))
                            self.normspuriouspixeach.append((len(self.flag)-sss)/(self.deachmaps[j].square**2))
                            self.normspuriousdeteach.append((len(self.flag)-sss)/len(self.flag))
                            self.normnotdetectedeach.append((len(self.txtarray[j])-sss)/len(self.txtarray[j]))
                    
                        else:
                            self.normspuriouseach.append(0)
                            self.normspuriousdeteach.append(0)
                            self.normnotdetectedeach.append(0)
                
                        self.efficiencyeach.append(sss/len(self.txtarray[j]))
                        #stateach.append(np.median(truedetectionfluxeach[j]))
                        for h in range(0,len(self.truedetectionfluxeach[j])):
                            if len(self.truedetectionfluxeach[j][h])>1:
                                self.truedetectionfluxeach[j][h]=np.mean(self.truedetectionfluxeach[j][h])
                        self.stateach.append(np.mean(self.truedetectionfluxeach[j]))
                        print('alldetection= {} ,True detection= {} ,spurious detections={}, not detected src={} '.format(len(self.flag),sss,len(self.flag)-sss,len(self.txtarray[j])-sss))
                        print('efficiency ={} normalized spurious detection number to all detection/number of generated sources={} normalized not detected source ={}'.format(self.efficiencyeach[j],self.normspuriouseach[j],self.normnotdetectedeach[j]))
                        print('END END END END END END END END END END END***************************************')
                         
                    self.truedetectioncoor.append(self.truedetectioncooreach)
                    self.stat.append(self.stateach)
                    self.truedetectionflux.append(self.truedetectionfluxeach)
                    self.truedetection.append(self.truedetectioneach)
                    self.finalflag.append(self.flag)
                    self.bilan.append(self.bilaneach)
                    self.detectedcoor.append(self.detectedcooreach) 
                    self.detectedcoordeg.append(self.detectedcoordegeach)
                    self.dmaps.append(self.deachmaps)
                    self.detectedsrc.append(self.detectedsrceach)
                    self.spurious.append(self.spuriouseach)
                    self.notdetected.append(self.notdetectedeach)
                    self.alldetection.append(self.alldetectioneach)
                    self.normspurious.append(self.normspuriouseach)
                    self.normspuriousdet.append(self.normspuriousdeteach)
                    self.normspuriouspix.append(self.normspuriouspixeach)
                    self.normnotdetected.append(self.normnotdetectedeach)
                    self.efficiency.append(self.efficiencyeach)


            
                    self.finalfinaldif.append(self.finaldif)
                    self.bilan.append(self.bilaneach)
            

                    self.txt.append(self.txteach)
            
                    self.txtname.append(self.txtnameeach) 
                    self.textarray.append(self.txtarray)

                    self.xytxt.append(self.xytxteach)
            
                    self.xytxtname.append(self.xytxtnameeach) 
                    self.xytextarray.append(self.xytxtarray)
            
            
            
                    self.deachmaps=0
                    self.deachmaps=[]
                    self.detectedcooreach=0
                    self.detectedcooreach=[]
            
                    self.spuriouseach=0
                    self.spuriouseach=[]
                    self.notdetectedeach=0
                    self.notdetectedeach=[]
                    self.alldetectioneach=0
                    self.alldetectioneach=[]
                    self.normspuriouseach=0
                    self.normspuriouseach=[]
                    self.normspuriousdeteach=0
                    self.normspuriousdeteach=[]
                    self.normspuriouspixeach=0
                    self.normspuriouspixeach=[]
                    self.normnotdetectedeach=0
                    self.normnotdetectedeach=[]
                    self.efficiencyeach=0
                    self.efficiencyeach=[]
                    self.bilaneach=0
                    self.bilaneach=[]
                    self.flag=0
                    self.flag=[]
                    self.truedetectioneach=0
                    self.truedetectioneach=[] 
                    self.truedetectioncooreach=0
                    self.truedetectioncooreach=[]
                    self.truedetectionfluxeach=0
                    self.truedetectionfluxeach=[]
                    self.detectedfluxeach=0
                    self.detectedfluxeach=[]
                    self.stateach=0
                    self.stateach=[]
            
                    self.bbilaneach=0
                    self.bbilaneach=[]

                    self.finalflag=0
                    self.finalflag=[]
                    self.finalangdist=0
                    self.finalangdist=[]
                    self.finaldif=0                           
                    self.finaldif=[]
                    self.matchsrc=0
                    self.matchsrc=[]
                    self.matchsrcdist=0
                    self.matchsrcdist=[]

                    self.txteach=0
                    self.txteach=[]
                    self.txtnameeach=0
                    self.txtnameeach=[]
                    self.deachmaps=0
                    self.deachmaps=[]
                    self.txtarray=0
                    self.txtarray=[]

                    self.xytxteach=0
                    self.xytxteach=[]
                    self.xytxtnameeach=0
                    self.xytxtnameeach=[]
                    self.xydeachmaps=0
                    self.xydeachmaps=[]
                    self.xytxtarray=0
                    self.xytxtarray=[]
            
           
                    self.detectedcoordegeach=0
                    self.detectedcoordegeach=[]
                    self.detectedsrceach=0
                    self.detectedsrceach=[]
                    self.spuriouseach=0
                    self.spuriouseach=[]
                    self.notdetectedeach=0
                    self.notdetectedeach=[]
                    self.alldetectioneach=0
                    self.alldetectioneach=[]
                    self.normspuriouseach=0
                    self.normspuriouseach=[]
                    self.normnotdetectedeach=0
                    self.normnotdetectedeach=[]
                    self.efficiencyeach=0
                    self.efficiencyeach=[]
                    self.bilaneach=0
                    self.bilaneach=[]

                self.finaltrue.append(self.true)
                self.finalstat.append(self.stat)
                self.finaltruedetectioncoor.append(self.truedetectioncoor)
                self.finaltruedetectionflux.append(self.truedetectionflux)
                self.finaltruedetection.append(self.truedetection)
                self.finalfinalflag.append(self.finalflag)
                self.finalbilan.append(self.bilan)
                self.finaldetectedcoor.append(self.detectedcoor)
                self.finaldetectedcoordeg.append(self.detectedcoordeg)
                self.finaldmaps.append(self.dmaps)
                self.finaldetectedsrc.append(self.detectedsrc)
                self.finalspurious.append(self.spurious)
                self.finalnotdetected.append(self.notdetected)
                self.finalalldetection.append(self.alldetection)
                self.finalnormspurious.append(self.normspurious)
                self.finalnormspuriousdet.append(self.normspuriousdet)
                self.finalnormspuriouspix.append(self.normspuriouspix)
                self.finalnormnotdetected.append(self.normnotdetected)
                self.finalefficiency.append(self.efficiency)
        

        
        
                self.finaltxt.append(self.txt)
                self.finaltxtname.append(self.txtname)
                self.finaltextarray.append(self.textarray)
                self.finaldetectedsrc.append(self.detectedsrc)

                self.xyfinaltxt.append(self.xytxt)
                self.xyfinaltxtname.append(self.xytxtname)
                self.xyfinaltextarray.append(self.xytextarray)
        

                self.bbilan=0
                self.bbilan=[]
        
                self.finalfinalflag=0
                self.finalfinalflag=[]
                self.finalfinalangdist=0
                self.finalfinalangdist=[]
                self.finalfinaldif=0                           
                self.finalfinaldif=[]
                self.finalmatchsrc=0
                self.finalmatchsrc=[]
                self.finalmatchsrcdist=0
                self.finalmatchsrcdist=[]
        
                self.txt=0
                self.txt=[]
                self.txtname=0
                self.txtname=[]
                self.xytxt=0
                self.xytxt=[]
                self.xytxtname=0
                self.xytxtname=[]
                self.txtfilesrcxy=0
                self.txtfilesrcxy=[]
                self.detectedcoor=0
                self.detectedcoor=[]
                self.detectedcoordeg=0
                self.detectedcoordeg=[]
                self.textarray=0
                self.textarray=[]
                self.xytextarray=0
                self.xytextarray=[]
                self.dmaps=0
                self.dmaps=[]
                self.detectedsrc=0
                self.detectedsrc=[]
                self.spurious=0
                self.spurious=[]
                self.notdetected=0
                self.notdetected=[]
                self.alldetection=0
                self.alldetection=[]
                self.normspurious=0
                self.normspurious=[]
                self.normspuriousdet=0
                self.normspuriousdet=[]
                self.normspuriouspix=0
                self.normspuriouspix=[]
                self.normnotdetected=0
                self.normnotdetected=[]
                self.efficiency=0
                self.efficiency=[]
                self.bilan=0
                self.bilan=[]
                self.detectedcoor=0
                self.detectedcoor=[]
        
                self.dmaps=0
                self.dmaps=[]
                self.detectedsrc=0
                self.detectedsrc=[]
                self.spurious=0
                self.spurious=[]
                self.notdetected=0
                self.notdetected=[]
                self.alldetection=0
                self.alldetection=[]
                self.finalflag=0
                self.finalflag=[]
                self.truedetection=0
                self.truedetection=[]
                self.truedetectioncoor=0
                self.truedetectioncoor=[]
                self.truedetectionflux=0
                self.truedetectionflux=[]
                self.stat=0
                self.stat=[]
                self.true=0
                self.true=[]

            self.finalfinalstat.append(self.finalstat)
            self.finalfinaltruedetectioncoor.append(self.finaltruedetectioncoor)
            self.finalfinaltruedetectionflux.append(self.finaltruedetectionflux)
            self.finalfinaltruedetection.append(self.finaltruedetection)
            self.finalfinalfinalflag.append(self.finalfinalflag)
            self.finalfinalbilan.append(self.finalbilan)
            self.finalfinaldetectedcoor.append(self.finaldetectedcoor)
            self.finalfinaldetectedcoordeg.append(self.finaldetectedcoordeg)
  
            self.finalfinaldmaps.append(self.finaldmaps)
            self.finalfinaldetectedsrc.append(self.finaldetectedsrc)
            self.finalfinalspurious.append(self.finalspurious)
            self.finalfinalnormspuriousdet.append(self.finalnormspuriousdet)
            self.finalfinalnormspuriouspix.append(self.finalnormspuriouspix)
            self.finalfinalnotdetected.append(self.finalnotdetected)
            self.finalfinalalldetection.append(self.finalalldetection)
            self.finalfinalnormspurious.append(self.finalnormspurious)
            self.finalfinalnormnotdetected.append(self.finalnormnotdetected)
            self.finalfinalefficiency.append(self.finalefficiency)  
            self.finalfinaltrue.append(self.finaltrue)

    

    
            self.finalfinaltxt.append(self.finaltxt)    
            self.finalfinaltxtname.append(self.finaltxtname) 
            self.finalfinaltextarray.append(self.finaltextarray)
            self.xyfinalfinaltxt.append(self.xyfinaltxt)    
            self.xyfinalfinaltxtname.append(self.xyfinaltxtname) 
            self.xyfinalfinaltextarray.append(self.xyfinaltextarray)

            self.finalfinalbilan.append(self.finalbilan)
   
            self.finalfinaldetectedsrc.append(self.finaldetectedsrc)
    
            self.finaltrue=0
            self.finaltrue=[]
            self.finalbbilan=0
            self.finalbbilan=[]
            self.finalfinalfinalflag=0
            self.finalfinalfinalflag=[]
            self.finalfinalfinalangdist=0
            self.finalfinalfinalangdist=[]
            self.finalfinalfinaldif=0                           
            self.finalfinalfinaldif=[]
            self.finalfinalmatchsrc=0
            self.finalfinalmatchsrc=[]
            self.finalfinalmatchsrcdist=0
            self.finalfinalmatchsrcdist=[] 
 
            self.finaltxt=0
            self.finaltxt=[]
            self.finaltxtname=0
            self.finaltxtname=[]
            self.xyfinaltxt=0
            self.xyfinaltxt=[]
            self.xyfinaltxtname=0
            self.xyfinaltxtname=[]
            self.finaldetectedcoor=0
            self.finaldetectedcoor=[]
            self.finaldetectedcoordeg=0
            self.finaldetectedcoordeg=[]
            self.finaltextarray=0
            self.finaltextarray=[]
            self.xyfinaltextarray=0
            self.xyfinaltextarray=[]
            self.finaldmaps=0
            self.finaldmaps=[]
            self.finaltxtfilesrcxy=0
            self.finaltxtfilesrcxy=[]
            self.finaldetectedsrc=0
            self.finaldetectedsrc=[]
            self.finalspurious=0
            self.finalspurious=[]
            self.finalnormspuriousdet=0
            self.finalnormspuriousdet=[]
            self.finalnormspuriouspix=0
            self.finalnormspuriouspix=[]
            self.finalnotdetected=0
            self.finalnotdetected=[]
            self.finalalldetection=0
            self.finalalldetection=[]
            self.finalnormspurious=0
            self.finalnormspurious=[]
            self.finalnormnotdetected=0
            self.finalnormnotdetected=[]
            self.finalefficiency=0
            self.finalefficiency=[]
            self.finalbilan=0
            self.finalbilan=[]
            self.finaltruedetection=0
            self.finaltruedetection=[]
            self.finalfinalflag=0
            self.finalfinalflag=[]
            self.finaltruedetectioncoor=0
            self.finaltruedetectioncoor=[]
            self.finaltruedetectionflux=0
            self.finaltruedetectionflux=[]
            self.finalstat=0
            self.finalstat=[]

    def run_efficiency_subtraction(self):
        #maps
        self.deachmaps=[]
        self.dmaps=[]
        self.finaldmaps=[]
        self.finalfinaldmaps=[]

        self.geachmaps=[]
        self.gmaps=[]
        self.finalgmaps=[]
        self.finalfinalgmaps=[]

        #????
        self.detectedsrceach=[]
        self.detectedsrc=[]
        self.finaldetectedsrc=[]
        self.finalfinaldetectedsrc=[]

        # spurious sources, normalized spurious sources
        self.spuriouseach=[]
        self.spurious=[]
        self.finalspurious=[]
        self.finalfinalspurious=[]


        self.normspuriouseach=[]
        self.normspurious=[]
        self.finalnormspurious=[]
        self.finalfinalnormspurious=[]


        self.normspuriouspixeach=[]
        self.normspuriouspix=[]
        self.finalnormspuriouspix=[]
        self.finalfinalnormspuriouspix=[]

        self.normspuriousdeteach=[]
        self.normspuriousdet=[]
        self.finalnormspuriousdet=[]
        self.finalfinalnormspuriousdet=[]

        #all detected coordinate in xy and radec
        self.detectedcooreach=[]
        self.detectedcoor=[]
        self.finaldetectedcoor=[]
        self.finalfinaldetectedcoor=[]

        self.detectedcoordegeach=[]
        self.detectedcoordeg=[]
        self.finaldetectedcoordeg=[]
        self.finalfinaldetectedcoordeg=[]

        #efficiency
        self.efficiencyeach=[]
        self.efficiency=[]
        self.finalefficiency=[]
        self.finalfinalefficiency=[]

        #detected coordinate from the list with flux
        self.truedetectioncooreacheach=[]
        self.truedetectioncooreach=[]
        self.truedetectioncoor=[]
        self.finaltruedetectioncoor=[]
        self.finalfinaltruedetectioncoor=[]

        self.truedetectionfluxeacheach=[]
        self.truedetectionfluxeach=[]
        self.truedetectionflux=[]
        self.finaltruedetectionflux=[]
        self.finalfinaltruedetectionflux=[]



        #?????
        self.truedetectioneach=[]
        self.truedetection=[]
        self.finaltruedetection=[]
        self.finalfinaltruedetection=[]

        #statistics (median,mean) for each map
        self.stateach=[]
        self.stat=[]
        self.finalstat=[]
        self.finalfinalstat=[]

        #flag
        self.finalflag=[]
        self.finalfinalflag=[]
        self.finalfinalfinalflag=[]

        #bilan
        self.bilaneach=[]
        self.bilan=[]
        self.finalbilan=[]
        self.finalfinalbilan=[]

        #???

        self.true=[]
        self.finaltrue=[]
        self.finalfinaltrue=[]

        # not detected src and normalized
        self.notdetectedeach=[]
        self.notdetected=[]
        self.finalnotdetected=[]
        self.finalfinalnotdetected=[]

        self.normnotdetectedeach=[]
        self.normnotdetected=[]
        self.finalnormnotdetected=[]
        self.finalfinalnormnotdetected=[]

        # all detection
        self.alldetection=[]
        self.alldetectioneach=[]
        self.finalalldetection=[]
        self.finalfinalalldetection=[]

        #textfiles
        self.txteach=[]
        self.txt=[]
        self.finaltxt=[]
        self.finalfinaltxt=[]

        self.txtarray=[]
        self.textarray=[]
        self.finaltextarray=[]
        self.finalfinaltextarray=[]

        self.txtnameeach=[]
        self.txtname=[]
        self.finaltxtname=[]
        self.finalfinaltxtname=[]

        self.xytxteach=[]
        self.xytxt=[]
        self.xyfinaltxt=[]
        self.xyfinalfinaltxt=[]

        self.xytxtarray=[]
        self.xytextarray=[]
        self.xyfinaltextarray=[]
        self.xyfinalfinaltextarray=[]


        self.xytxtnameeach=[]
        self.xytxtname=[]
        self.xyfinaltxtname=[]
        self.xyfinalfinaltxtname=[]



        self.detectedfluxeach=[]
        self.finaldif=[]
        self.finalfinaldif=[]
        self.finalfinalfinaldif=[]
        #generated/detected comparison algorithm
        

        for noisex in range(0,len(self.signoise)):
            for nx in range(0,len(self.n)):
                print('for a given false detection rate n={}'.format(self.n[nx]))
                for i in range(0,len(self.jandec)):
                    for j in range(0,self.iteration):
                        #if self.freq_subtraction[freqx]==1300:
                        self.txtnameeach.append('offcentral_{}_{}.txt'.format(np.round(self.jandec[i],3),j))
                        #else:
                        #self.txtnameeach.append('nvss_src_S1000_dec80.txt')
                        self.txteach.append(open(self.txtnameeach[j],'r'))
                        #self.xytxteach.append(open(self.xytxtnameeach[j],'r'))
                        self.txtarray.append(np.loadtxt(self.txteach[j]))
                        #self.xytxtarray.append(np.loadtxt(self.xytxteach[j]))
                        self.txteach[j].close()
                        #self.xytxteach[j].close()
                        self.deachmaps.append(extract_class.extract())
                        filenameleft='filtmap_ncp_autoon_{}_{}_4dec.fits'.format(self.signoise1[noisex],self.freq_subtraction[0])
                        filenamecenter='filtmap_ncp_autoon_{}_{}_{}_{}_4dec.fits'.format(self.signoise1[noisex],self.freq_subtraction[1],np.round(self.jandec[i],3),j)
                        filenameright='filtmap_ncp_autoon_{}_{}_4dec.fits'.format(self.signoise1[noisex],self.freq_subtraction[2])
                        self.fullmapleft=hp.read_map(filenameleft)
                        self.fullmapcenter=hp.read_map(filenamecenter)
                        self.fullmapright=hp.read_map(filenameright)
                        self.fullmap=self.fullmapcenter-0.5*(self.fullmapleft+self.fullmapright)
                        self.deachmaps[j].set_filename(self.finalgnames[noisex][i][j])
                        self.deachmaps[j].set_size(90)
                        self.deachmaps[j].set_n_sigma(self.n[nx])
                        self.deachmaps[j].set_detuse('skymaps')
                        self.deachmaps[j].set_k(0)
                        self.deachmaps[j].set_query_disk_radius(74)
                        self.deachmaps[j].set_freq(1300*10**6)
                        self.deachmaps[j].set_maxbaseline(16.5)
                        self.deachmaps[j].set_reso(12)
                        self.deachmaps[j].set_nker1(3)
                        self.deachmaps[j].set_nker2(5)
                        self.deachmaps[j].set_radius(np.sqrt(2*(2.5**2)))
                        self.deachmaps[j].set_nker3(7)
                        self.deachmaps[j].set_nker4(4)
                        self.deachmaps[j].set_nker5(8)
                        self.deachmaps[j].set_nker6(10)
                        self.deachmaps[j].extract()
                        rectmap=hp.gnomview(self.fullmap,rot=[self.deachmaps[j].radeg_center,self.deachmaps[j].decdeg_center],reso=self.deachmaps[j].reso,xsize=self.deachmaps[j].size,ysize=self.deachmaps[j].size,return_projected_map=True,no_plot=True)
                        self.deachmaps[j].rectmap=rectmap
                        self.deachmaps[j].detect()


                        print('START START START START START START START ********************************************')
                        #print('srcnumber={} noise= {} n={} flux={} iteration={} '.format(,signoise[sigx],n[nx],jandec[i],j ))
                
                        self.detectedcooreach.append(np.transpose((self.deachmaps[j].final_list[:,1],self.deachmaps[j].final_list[:,2],self.deachmaps[j].final_list[:,5],self.deachmaps[j].final_list[:,6],self.deachmaps[j].final_list[:,10],self.deachmaps[j].final_list[:,3],self.deachmaps[j].final_list[:,4])))
                        self.detectedcoordegeach.append(np.transpose((self.deachmaps[j].final_list[:,3],self.deachmaps[j].final_list[:,4],self.deachmaps[j].final_list[:,5],self.deachmaps[j].final_list[:,6],self.deachmaps[j].final_list[:,10])))
                        self.detectedfluxeach.append(self.deachmaps[j].final_list[:,6])

                        print('generated src')
                        print(self.txtarray[j])
                        print('all detected src')
                        print(self.detectedcooreach[j])
                        self.flag=np.zeros(np.shape(self.detectedcoordegeach[j])[0])
                            
                        print('computing difference between each generated coor and all detected coor')
                        ####
                        self.transposetxtarray=np.transpose((self.txtarray[j][:,0],self.txtarray[j][:,1]))
                        #self.xytransposetxtarray=np.transpose((self.xytxtarray[j][:,0],self.xytxtarray[j][:,1]))
                        for o in range(0,np.shape(self.transposetxtarray)[0]):
                            #difeach.append(transposetxtarray[o]-detectedcoordegeach[j][:,[0,1]])
                    
                    
                            self.dif=((180/math.pi)*np.arccos(np.sin(self.transposetxtarray[o][1]*(math.pi/180))*np.sin(self.detectedcoordegeach[j][:,1]*(math.pi/180))+np.cos(self.transposetxtarray[o][1]*(math.pi/180))*np.cos(self.detectedcoordegeach[j][:,1]*(math.pi/180))*np.cos(self.transposetxtarray[o][0]*(math.pi/180)-self.detectedcoordegeach[j][:,0]*(math.pi/180))))
                            #self.dif2=self.xytransposetxtarray[o]-self.detectedcooreach[j][:,[0,1]]
                            print('dif')
                            print(self.dif)
                            print(np.where(abs(self.dif)<=self.delta))
                            #matchsrceach.append(detectedcoordegeach[j][np.where(angdist[o]<=delta)])
                            #matchsrcfluxeach.append(detectedfluxeach[j][np.where(angdist[o]<=delta)])
                            #matchsrcdisteach.append(angdist[o][np.where(angdist[o]<=delta)])

                            if len(np.where(abs(self.dif)<=self.delta)[0])!=0:        
                                self.truedetectioncooreacheach.append(self.detectedcooreach[j][np.where( abs(self.dif)<=self.delta )])
                                self.truedetectionfluxeacheach.append(self.detectedfluxeach[j][np.where( abs(self.dif)<=self.delta )])
                                #if len(np.where(abs(dif)<=delta)[0])==1:
                                self.true.append(self.detectedfluxeach[j][np.where( abs(self.dif)<=self.delta )])
                            self.flag[ np.where(abs(self.dif)<=self.delta) ] = self.flag[ np.where(abs(self.dif)<=self.delta) ] +1
                    
                            print('**** **** **** ** **** * * * ** * * ** * * * ** * * * ')
                            print('second flag')
                            print(self.flag)
                            print(self.detectedfluxeach[j][np.where( abs(self.dif)<=self.delta )])
                            print('* * * * ** * * * ** * * * * * * ** * * * * * ** * * * *')

                            s=0
                            ss=0
                            sss=0
                            for k in self.flag:
                                if k==1:
                                    s=s+1
                                if k>1:
                                    ss=ss+k
                                sss=s+ss
                            print('{} one detection/one src'.format(s))
                            print('{} one detection for  many sources'.format(ss))
                            print('total {}'.format(sss))

                        self.truedetectioncooreach.append(self.truedetectioncooreacheach)
                        self.truedetectionfluxeach.append(self.truedetectionfluxeacheach)
                        self.truedetectioncooreacheach=0
                        self.truedetectioncooreacheach=[]
                        self.truedetectionfluxeacheach=0
                        self.truedetectionfluxeacheach=[]
                        self.bilaneach.append([s,ss,sss])
                        self.detectedsrceach.append(sss)        
                        self.spuriouseach.append(len(self.flag)-sss)
                        self.truedetectioncooreacheach2=0
                        self.truedetectioncooreacheach2=[]
                        self.truedetectionfluxeacheach2=0
                        self.truedetectionfluxeacheach2=[]
                        self.notdetectedeach.append(len(self.txtarray[j])-sss)
    
                        self.alldetectioneach.append(len(self.flag))
                        if len(self.flag)!=0:
                            self.normspuriouseach.append((len(self.flag)-sss)/len(self.txtarray[j]))
                            self.normspuriouspixeach.append((len(self.flag)-sss)/(self.deachmaps[j].square**2))
                            self.normspuriousdeteach.append((len(self.flag)-sss)/len(self.flag))
                            self.normnotdetectedeach.append((len(self.txtarray[j])-sss)/len(self.txtarray[j]))
                    
                        else:
                            self.normspuriouseach.append(0)
                            self.normspuriousdeteach.append(0)
                            self.normnotdetectedeach.append(0)
                
                        self.efficiencyeach.append(sss/len(self.txtarray[j]))
                        #stateach.append(np.median(truedetectionfluxeach[j]))
                        for h in range(0,len(self.truedetectionfluxeach[j])):
                            if len(self.truedetectionfluxeach[j][h])>1:
                                self.truedetectionfluxeach[j][h]=np.mean(self.truedetectionfluxeach[j][h])
                        self.stateach.append(np.mean(self.truedetectionfluxeach[j]))
                        print('alldetection= {} ,True detection= {} ,spurious detections={}, not detected src={} '.format(len(self.flag),sss,len(self.flag)-sss,len(self.txtarray[j])-sss))
                        print('efficiency ={} normalized spurious detection number to all detection/number of generated sources={} normalized not detected source ={}'.format(self.efficiencyeach[j],self.normspuriouseach[j],self.normnotdetectedeach[j]))
                        print('END END END END END END END END END END END***************************************')
                         
                    self.truedetectioncoor.append(self.truedetectioncooreach)
                    self.stat.append(self.stateach)
                    self.truedetectionflux.append(self.truedetectionfluxeach)
                    self.truedetection.append(self.truedetectioneach)
                    self.finalflag.append(self.flag)
                    self.bilan.append(self.bilaneach)
                    self.detectedcoor.append(self.detectedcooreach) 
                    self.detectedcoordeg.append(self.detectedcoordegeach)
                    self.dmaps.append(self.deachmaps)
                    self.detectedsrc.append(self.detectedsrceach)
                    self.spurious.append(self.spuriouseach)
                    self.notdetected.append(self.notdetectedeach)
                    self.alldetection.append(self.alldetectioneach)
                    self.normspurious.append(self.normspuriouseach)
                    self.normspuriousdet.append(self.normspuriousdeteach)
                    self.normspuriouspix.append(self.normspuriouspixeach)
                    self.normnotdetected.append(self.normnotdetectedeach)
                    self.efficiency.append(self.efficiencyeach)


            
                    self.finalfinaldif.append(self.finaldif)
                    self.bilan.append(self.bilaneach)
            

                    self.txt.append(self.txteach)
            
                    self.txtname.append(self.txtnameeach) 
                    self.textarray.append(self.txtarray)

                    self.xytxt.append(self.xytxteach)
            
                    self.xytxtname.append(self.xytxtnameeach) 
                    self.xytextarray.append(self.xytxtarray)
            
            
            
                    self.deachmaps=0
                    self.deachmaps=[]
                    self.detectedcooreach=0
                    self.detectedcooreach=[]
            
                    self.spuriouseach=0
                    self.spuriouseach=[]
                    self.notdetectedeach=0
                    self.notdetectedeach=[]
                    self.alldetectioneach=0
                    self.alldetectioneach=[]
                    self.normspuriouseach=0
                    self.normspuriouseach=[]
                    self.normspuriousdeteach=0
                    self.normspuriousdeteach=[]
                    self.normspuriouspixeach=0
                    self.normspuriouspixeach=[]
                    self.normnotdetectedeach=0
                    self.normnotdetectedeach=[]
                    self.efficiencyeach=0
                    self.efficiencyeach=[]
                    self.bilaneach=0
                    self.bilaneach=[]
                    self.flag=0
                    self.flag=[]
                    self.truedetectioneach=0
                    self.truedetectioneach=[] 
                    self.truedetectioncooreach=0
                    self.truedetectioncooreach=[]
                    self.truedetectionfluxeach=0
                    self.truedetectionfluxeach=[]
                    self.detectedfluxeach=0
                    self.detectedfluxeach=[]
                    self.stateach=0
                    self.stateach=[]
            
                    self.bbilaneach=0
                    self.bbilaneach=[]

                    self.finalflag=0
                    self.finalflag=[]
                    self.finalangdist=0
                    self.finalangdist=[]
                    self.finaldif=0                           
                    self.finaldif=[]
                    self.matchsrc=0
                    self.matchsrc=[]
                    self.matchsrcdist=0
                    self.matchsrcdist=[]

                    self.txteach=0
                    self.txteach=[]
                    self.txtnameeach=0
                    self.txtnameeach=[]
                    self.deachmaps=0
                    self.deachmaps=[]
                    self.txtarray=0
                    self.txtarray=[]

                    self.xytxteach=0
                    self.xytxteach=[]
                    self.xytxtnameeach=0
                    self.xytxtnameeach=[]
                    self.xydeachmaps=0
                    self.xydeachmaps=[]
                    self.xytxtarray=0
                    self.xytxtarray=[]
            
           
                    self.detectedcoordegeach=0
                    self.detectedcoordegeach=[]
                    self.detectedsrceach=0
                    self.detectedsrceach=[]
                    self.spuriouseach=0
                    self.spuriouseach=[]
                    self.notdetectedeach=0
                    self.notdetectedeach=[]
                    self.alldetectioneach=0
                    self.alldetectioneach=[]
                    self.normspuriouseach=0
                    self.normspuriouseach=[]
                    self.normnotdetectedeach=0
                    self.normnotdetectedeach=[]
                    self.efficiencyeach=0
                    self.efficiencyeach=[]
                    self.bilaneach=0
                    self.bilaneach=[]

                self.finaltrue.append(self.true)
                self.finalstat.append(self.stat)
                self.finaltruedetectioncoor.append(self.truedetectioncoor)
                self.finaltruedetectionflux.append(self.truedetectionflux)
                self.finaltruedetection.append(self.truedetection)
                self.finalfinalflag.append(self.finalflag)
                self.finalbilan.append(self.bilan)
                self.finaldetectedcoor.append(self.detectedcoor)
                self.finaldetectedcoordeg.append(self.detectedcoordeg)
                self.finaldmaps.append(self.dmaps)
                self.finaldetectedsrc.append(self.detectedsrc)
                self.finalspurious.append(self.spurious)
                self.finalnotdetected.append(self.notdetected)
                self.finalalldetection.append(self.alldetection)
                self.finalnormspurious.append(self.normspurious)
                self.finalnormspuriousdet.append(self.normspuriousdet)
                self.finalnormspuriouspix.append(self.normspuriouspix)
                self.finalnormnotdetected.append(self.normnotdetected)
                self.finalefficiency.append(self.efficiency)
        

        
        
                self.finaltxt.append(self.txt)
                self.finaltxtname.append(self.txtname)
                self.finaltextarray.append(self.textarray)
                self.finaldetectedsrc.append(self.detectedsrc)

                self.xyfinaltxt.append(self.xytxt)
                self.xyfinaltxtname.append(self.xytxtname)
                self.xyfinaltextarray.append(self.xytextarray)
        

                self.bbilan=0
                self.bbilan=[]
        
                self.finalfinalflag=0
                self.finalfinalflag=[]
                self.finalfinalangdist=0
                self.finalfinalangdist=[]
                self.finalfinaldif=0                           
                self.finalfinaldif=[]
                self.finalmatchsrc=0
                self.finalmatchsrc=[]
                self.finalmatchsrcdist=0
                self.finalmatchsrcdist=[]
        
                self.txt=0
                self.txt=[]
                self.txtname=0
                self.txtname=[]
                self.xytxt=0
                self.xytxt=[]
                self.xytxtname=0
                self.xytxtname=[]
                self.txtfilesrcxy=0
                self.txtfilesrcxy=[]
                self.detectedcoor=0
                self.detectedcoor=[]
                self.detectedcoordeg=0
                self.detectedcoordeg=[]
                self.textarray=0
                self.textarray=[]
                self.xytextarray=0
                self.xytextarray=[]
                self.dmaps=0
                self.dmaps=[]
                self.detectedsrc=0
                self.detectedsrc=[]
                self.spurious=0
                self.spurious=[]
                self.notdetected=0
                self.notdetected=[]
                self.alldetection=0
                self.alldetection=[]
                self.normspurious=0
                self.normspurious=[]
                self.normspuriousdet=0
                self.normspuriousdet=[]
                self.normspuriouspix=0
                self.normspuriouspix=[]
                self.normnotdetected=0
                self.normnotdetected=[]
                self.efficiency=0
                self.efficiency=[]
                self.bilan=0
                self.bilan=[]
                self.detectedcoor=0
                self.detectedcoor=[]
        
                self.dmaps=0
                self.dmaps=[]
                self.detectedsrc=0
                self.detectedsrc=[]
                self.spurious=0
                self.spurious=[]
                self.notdetected=0
                self.notdetected=[]
                self.alldetection=0
                self.alldetection=[]
                self.finalflag=0
                self.finalflag=[]
                self.truedetection=0
                self.truedetection=[]
                self.truedetectioncoor=0
                self.truedetectioncoor=[]
                self.truedetectionflux=0
                self.truedetectionflux=[]
                self.stat=0
                self.stat=[]
                self.true=0
                self.true=[]

            self.finalfinalstat.append(self.finalstat)
            self.finalfinaltruedetectioncoor.append(self.finaltruedetectioncoor)
            self.finalfinaltruedetectionflux.append(self.finaltruedetectionflux)
            self.finalfinaltruedetection.append(self.finaltruedetection)
            self.finalfinalfinalflag.append(self.finalfinalflag)
            self.finalfinalbilan.append(self.finalbilan)
            self.finalfinaldetectedcoor.append(self.finaldetectedcoor)
            self.finalfinaldetectedcoordeg.append(self.finaldetectedcoordeg)
  
            self.finalfinaldmaps.append(self.finaldmaps)
            self.finalfinaldetectedsrc.append(self.finaldetectedsrc)
            self.finalfinalspurious.append(self.finalspurious)
            self.finalfinalnormspuriousdet.append(self.finalnormspuriousdet)
            self.finalfinalnormspuriouspix.append(self.finalnormspuriouspix)
            self.finalfinalnotdetected.append(self.finalnotdetected)
            self.finalfinalalldetection.append(self.finalalldetection)
            self.finalfinalnormspurious.append(self.finalnormspurious)
            self.finalfinalnormnotdetected.append(self.finalnormnotdetected)
            self.finalfinalefficiency.append(self.finalefficiency)  
            self.finalfinaltrue.append(self.finaltrue)

    

    
            self.finalfinaltxt.append(self.finaltxt)    
            self.finalfinaltxtname.append(self.finaltxtname) 
            self.finalfinaltextarray.append(self.finaltextarray)
            self.xyfinalfinaltxt.append(self.xyfinaltxt)    
            self.xyfinalfinaltxtname.append(self.xyfinaltxtname) 
            self.xyfinalfinaltextarray.append(self.xyfinaltextarray)

            self.finalfinalbilan.append(self.finalbilan)
   
            self.finalfinaldetectedsrc.append(self.finaldetectedsrc)
    
            self.finaltrue=0
            self.finaltrue=[]
            self.finalbbilan=0
            self.finalbbilan=[]
            self.finalfinalfinalflag=0
            self.finalfinalfinalflag=[]
            self.finalfinalfinalangdist=0
            self.finalfinalfinalangdist=[]
            self.finalfinalfinaldif=0                           
            self.finalfinalfinaldif=[]
            self.finalfinalmatchsrc=0
            self.finalfinalmatchsrc=[]
            self.finalfinalmatchsrcdist=0
            self.finalfinalmatchsrcdist=[] 
 
            self.finaltxt=0
            self.finaltxt=[]
            self.finaltxtname=0
            self.finaltxtname=[]
            self.xyfinaltxt=0
            self.xyfinaltxt=[]
            self.xyfinaltxtname=0
            self.xyfinaltxtname=[]
            self.finaldetectedcoor=0
            self.finaldetectedcoor=[]
            self.finaldetectedcoordeg=0
            self.finaldetectedcoordeg=[]
            self.finaltextarray=0
            self.finaltextarray=[]
            self.xyfinaltextarray=0
            self.xyfinaltextarray=[]
            self.finaldmaps=0
            self.finaldmaps=[]
            self.finaltxtfilesrcxy=0
            self.finaltxtfilesrcxy=[]
            self.finaldetectedsrc=0
            self.finaldetectedsrc=[]
            self.finalspurious=0
            self.finalspurious=[]
            self.finalnormspuriousdet=0
            self.finalnormspuriousdet=[]
            self.finalnormspuriouspix=0
            self.finalnormspuriouspix=[]
            self.finalnotdetected=0
            self.finalnotdetected=[]
            self.finalalldetection=0
            self.finalalldetection=[]
            self.finalnormspurious=0
            self.finalnormspurious=[]
            self.finalnormnotdetected=0
            self.finalnormnotdetected=[]
            self.finalefficiency=0
            self.finalefficiency=[]
            self.finalbilan=0
            self.finalbilan=[]
            self.finaltruedetection=0
            self.finaltruedetection=[]
            self.finalfinalflag=0
            self.finalfinalflag=[]
            self.finaltruedetectioncoor=0
            self.finaltruedetectioncoor=[]
            self.finaltruedetectionflux=0
            self.finaltruedetectionflux=[]
            self.finalstat=0
            self.finalstat=[]




    def run_efficiency_pythmaps(self):
        #maps
        self.deachmaps=[]
        self.dmaps=[]
        self.finaldmaps=[]
        self.finalfinaldmaps=[]

        self.geachmaps=[]
        self.gmaps=[]
        self.finalgmaps=[]
        self.finalfinalgmaps=[]

        #????
        self.detectedsrceach=[]
        self.detectedsrc=[]
        self.finaldetectedsrc=[]
        self.finalfinaldetectedsrc=[]

        # spurious sources, normalized spurious sources
        self.spuriouseach=[]
        self.spurious=[]
        self.finalspurious=[]
        self.finalfinalspurious=[]


        self.normspuriouseach=[]
        self.normspurious=[]
        self.finalnormspurious=[]
        self.finalfinalnormspurious=[]


        self.normspuriouspixeach=[]
        self.normspuriouspix=[]
        self.finalnormspuriouspix=[]
        self.finalfinalnormspuriouspix=[]

        self.normspuriousdeteach=[]
        self.normspuriousdet=[]
        self.finalnormspuriousdet=[]
        self.finalfinalnormspuriousdet=[]

        #all detected coordinate in xy and radec
        self.detectedcooreach=[]
        self.detectedcoor=[]
        self.finaldetectedcoor=[]
        self.finalfinaldetectedcoor=[]

        self.detectedcoordegeach=[]
        self.detectedcoordeg=[]
        self.finaldetectedcoordeg=[]
        self.finalfinaldetectedcoordeg=[]

        #efficiency
        self.efficiencyeach=[]
        self.efficiency=[]
        self.finalefficiency=[]
        self.finalfinalefficiency=[]

        #detected coordinate from the list with flux
        self.truedetectioncooreacheach=[]
        self.truedetectioncooreach=[]
        self.truedetectioncoor=[]
        self.finaltruedetectioncoor=[]
        self.finalfinaltruedetectioncoor=[]

        self.truedetectionfluxeacheach=[]
        self.truedetectionfluxeach=[]
        self.truedetectionflux=[]
        self.finaltruedetectionflux=[]
        self.finalfinaltruedetectionflux=[]



        #?????
        self.truedetectioneach=[]
        self.truedetection=[]
        self.finaltruedetection=[]
        self.finalfinaltruedetection=[]

        #statistics (median,mean) for each map
        self.stateach=[]
        self.stat=[]
        self.finalstat=[]
        self.finalfinalstat=[]

        #flag
        self.finalflag=[]
        self.finalfinalflag=[]
        self.finalfinalfinalflag=[]

        #bilan
        self.bilaneach=[]
        self.bilan=[]
        self.finalbilan=[]
        self.finalfinalbilan=[]

        #???

        self.true=[]
        self.finaltrue=[]
        self.finalfinaltrue=[]

        # not detected src and normalized
        self.notdetectedeach=[]
        self.notdetected=[]
        self.finalnotdetected=[]
        self.finalfinalnotdetected=[]

        self.normnotdetectedeach=[]
        self.normnotdetected=[]
        self.finalnormnotdetected=[]
        self.finalfinalnormnotdetected=[]

        # all detection
        self.alldetection=[]
        self.alldetectioneach=[]
        self.finalalldetection=[]
        self.finalfinalalldetection=[]

        #textfiles
        self.txteach=[]
        self.txt=[]
        self.finaltxt=[]
        self.finalfinaltxt=[]

        self.txtarray=[]
        self.textarray=[]
        self.finaltextarray=[]
        self.finalfinaltextarray=[]

        self.txtnameeach=[]
        self.txtname=[]
        self.finaltxtname=[]
        self.finalfinaltxtname=[]

        self.xytxteach=[]
        self.xytxt=[]
        self.xyfinaltxt=[]
        self.xyfinalfinaltxt=[]

        self.xytxtarray=[]
        self.xytextarray=[]
        self.xyfinaltextarray=[]
        self.xyfinalfinaltextarray=[]


        self.xytxtnameeach=[]
        self.xytxtname=[]
        self.xyfinaltxtname=[]
        self.xyfinalfinaltxtname=[]



        self.detectedfluxeach=[]
        self.finaldif=[]
        self.finalfinaldif=[]
        self.finalfinalfinaldif=[]
        #generated/detected comparison algorithm
        

        for sigx in range(0,len(self.signoise)):
            for nx in range(0,len(self.n)):
                print('for a given false detection rate n={}'.format(self.n[nx]))
                for i in range(0,len(self.jandec)):
                    for j in range(0,self.iteration):
                        self.geachmaps.append(src_generator.generate(self.jandec[i]))
                        self.geachmaps[j].set_sigma_noise(self.signoise[sigx])
                        self.geachmaps[j].set_src_sigma(self.sigmag)
                        self.geachmaps[j].set_map_size(90)
                        self.geachmaps[j].set_side(10)
                        self.geachmaps[j].generator()
                        self.deachmaps.append(extract_class.extract())
                        self.deachmaps[j].set_reso(12)
                        self.deachmaps[j].set_size(90)
                        self.deachmaps[j].rectmap=self.geachmaps[j].rectmap
                        self.deachmaps[j].set_n_sigma(self.n[nx])
                        self.deachmaps[j].set_detuse('pythmaps')
                        self.deachmaps[j].set_radius(np.sqrt(2*(2.5**2)))
                        self.deachmaps[j].set_k(0)
                        self.deachmaps[j].set_square(74)
                        self.deachmaps[j].detect()
                        print('START START START START START START START ********************************************')
                        print('srcnumber={} noise= {} n={} flux={} iteration={} '.format(self.geachmaps[j].srcnum,self.signoise[sigx],self.n[nx],self.jandec[i],j ))
                
                        self.detectedcooreach.append(np.transpose((self.deachmaps[j].final_list[:,1],self.deachmaps[j].final_list[:,2],self.deachmaps[j].final_list[:,5],self.deachmaps[j].final_list[:,6])))
                        self.detectedfluxeach.append(self.deachmaps[j].final_list[:,6])
                        #detectedcooreach=np.transpose((deachmaps[j].final_list[:,1],deachmaps[j].final_list[:,2],deachmaps[j].final_list[:,5],deachmaps[j].final_list[:,6]))
                
                
                        print('generated src')
                        print(self.geachmaps[j].transpose_h_center_pixcenter)
                        print('all detected src')
                        print(self.detectedcooreach[j])
                        self.flag=np.zeros(np.shape(self.detectedcooreach[j])[0])
               
                        print('computing difference between each generated coor and all detected coor')
   
        
                
                        for o in range(0,np.shape(self.geachmaps[j].transpose_h_center_pixcenter)[0]):
                            #dif=geachmaps[j].transpose_h_center_pixcenter[o]-detectedcooreach
                    
                            self.dif=self.geachmaps[j].transpose_h_center_pixcenter[o]-self.detectedcooreach[j][:,[0,1]]
                            #dif=geachmaps[j].transpose_h_center_pixcenter[o]-detectedcooreach[:,[0,1]]
                            print('coor dif with {}'.format(self.geachmaps[j].transpose_h_center_pixcenter[o]))
                            #print(dif)
                            print('flag before')
                            print(self.flag)
                    
                            if len(np.where( (abs(self.dif[:,0])<=1) & (abs(self.dif[:,1])<=1))[0])!=0:        
                                self.truedetectioncooreacheach.append(self.detectedcooreach[j][np.where( (abs(self.dif[:,0])<=1) & (abs(self.dif[:,1])<=1))])
                                self.truedetectionfluxeacheach.append(self.detectedfluxeach[j][np.where( (abs(self.dif[:,0])<=1) & (abs(self.dif[:,1])<=1))])
                                #if len(np.where( (abs(dif[:,0])<=1) & (abs(dif[:,1])<=1))[0])==1:
                                self.true.append(self.detectedfluxeach[j][np.where( (abs(self.dif[:,0])<=1) & (abs(self.dif[:,1])<=1))])
                            self.flag[  np.where( (abs(self.dif[:,0])<=1) & (abs(self.dif[:,1])<=1) )  ]=self.flag[  np.where( (abs(self.dif[:,0])<=1) & (abs(self.dif[:,1])<=1) )  ]+1
                            print('flag after')
                            print(self.flag)
                            
                            s=0
                            ss=0
                            sss=0
                            for k in self.flag:
                                if k==1:
                                    s=s+1
                                if k>1:
                                    ss=ss+k
                                sss=s+ss
                            print('{} one detection/one src'.format(s))
                            print('{} one detection for  many sources'.format(ss))
                            print('total {}'.format(sss))

                        self.truedetectioncooreach.append(self.truedetectioncooreacheach)
                        self.truedetectionfluxeach.append(self.truedetectionfluxeacheach)
                        self.truedetectioncooreacheach=0
                        self.truedetectioncooreacheach=[]
                        self.truedetectionfluxeacheach=0
                        self.truedetectionfluxeacheach=[]
                        self.bilaneach.append([s,ss,sss])
                        self.detectedsrceach.append(sss)        
                        self.spuriouseach.append(len(self.flag)-sss)
                        self.truedetectioncooreacheach2=0
                        self.truedetectioncooreacheach2=[]
                        self.truedetectionfluxeacheach2=0
                        self.truedetectionfluxeacheach2=[]
                        self.notdetectedeach.append(self.geachmaps[j].srcnum-sss)
    
                        self.alldetectioneach.append(len(self.flag))
                        if len(self.flag)!=0:
                            self.normspuriouseach.append((len(self.flag)-sss)/self.geachmaps[j].srcnum)
                            self.normspuriouspixeach.append((len(self.flag)-sss)/(self.deachmaps[j].square**2))
                            self.normspuriousdeteach.append((len(self.flag)-sss)/len(self.flag))
                            self.normnotdetectedeach.append((self.geachmaps[j].srcnum-sss)/self.geachmaps[j].srcnum)
                    
                        else:
                            self.normspuriouseach.append(0)
                            self.normspuriousdeteach.append(0)
                            self.normnotdetectedeach.append(0)
                
                        self.efficiencyeach.append(sss/self.geachmaps[j].srcnum)
                        #stateach.append(np.median(truedetectionfluxeach[j]))
                        for h in range(0,len(self.truedetectionfluxeach[j])):
                            if len(self.truedetectionfluxeach[j][h])>1:
                                self.truedetectionfluxeach[j][h]=np.mean(self.truedetectionfluxeach[j][h])
                        self.stateach.append(np.mean(self.truedetectionfluxeach[j]))
                        print('alldetection= {} ,True detection= {} ,spurious detections={}, not detected src={} '.format(len(self.flag),sss,len(self.flag)-sss,self.geachmaps[j].srcnum-sss))
                        print('efficiency ={} normalized spurious detection number to all detection/number of generated sources={} normalized not detected source ={}'.format(self.efficiencyeach[j],self.normspuriouseach[j],self.normnotdetectedeach[j]))
                        print('END END END END END END END END END END END***************************************')
                         
                    self.truedetectioncoor.append(self.truedetectioncooreach)
                    self.stat.append(self.stateach)
                    self.truedetectionflux.append(self.truedetectionfluxeach)
                    self.truedetection.append(self.truedetectioneach)
                    self.finalflag.append(self.flag)
                    self.bilan.append(self.bilaneach)
                    self.detectedcoor.append(self.detectedcooreach) 
                    self.gmaps.append(self.geachmaps)
                    self.dmaps.append(self.deachmaps)
                    self.detectedsrc.append(self.detectedsrceach)
                    self.spurious.append(self.spuriouseach)
                    self.notdetected.append(self.notdetectedeach)
                    self.alldetection.append(self.alldetectioneach)
                    self.normspurious.append(self.normspuriouseach)
                    self.normspuriousdet.append(self.normspuriousdeteach)
                    self.normspuriouspix.append(self.normspuriouspixeach)
                    self.normnotdetected.append(self.normnotdetectedeach)
                    self.efficiency.append(self.efficiencyeach)
                    
                    self.geachmaps=0
                    self.geachmaps=[]
                    self.deachmaps=0
                    self.deachmaps=[]
                    self.detectedcooreach=0
                    self.detectedcooreach=[]
            
                    self.spuriouseach=0
                    self.spuriouseach=[]
                    self.notdetectedeach=0
                    self.notdetectedeach=[]
                    self.alldetectioneach=0
                    self.alldetectioneach=[]
                    self.normspuriouseach=0
                    self.normspuriouseach=[]
                    self.normspuriousdeteach=0
                    self.normspuriousdeteach=[]
                    self.normspuriouspixeach=0
                    self.normspuriouspixeach=[]
                    self.normnotdetectedeach=0
                    self.normnotdetectedeach=[]
                    self.efficiencyeach=0
                    self.efficiencyeach=[]
                    self.bilaneach=0
                    self.bilaneach=[]
                    self.flag=0
                    self.flag=[]
                    self.truedetectioneach=0
                    self.truedetectioneach=[] 
                    self.truedetectioncooreach=0
                    self.truedetectioncooreach=[]
                    self.truedetectionfluxeach=0
                    self.truedetectionfluxeach=[]
                    self.detectedfluxeach=0
                    self.detectedfluxeach=[]
                    self.stateach=0
                    self.stateach=[]
            
            

                self.finaltrue.append(self.true)
                self.finalstat.append(self.stat)
                self.finaltruedetectioncoor.append(self.truedetectioncoor)
                self.finaltruedetectionflux.append(self.truedetectionflux)
                self.finaltruedetection.append(self.truedetection)
                self.finalfinalflag.append(self.finalflag)
                self.finalbilan.append(self.bilan)
                self.finaldetectedcoor.append(self.detectedcoor)
                self.finaldetectedcoordeg.append(self.detectedcoordeg)
                self.finalgmaps.append(self.gmaps)
                self.finaldmaps.append(self.dmaps)
                self.finaldetectedsrc.append(self.detectedsrc)
                self.finalspurious.append(self.spurious)
                self.finalnotdetected.append(self.notdetected)
                self.finalalldetection.append(self.alldetection)
                self.finalnormspurious.append(self.normspurious)
                self.finalnormspuriousdet.append(self.normspuriousdet)
                self.finalnormspuriouspix.append(self.normspuriouspix)
                self.finalnormnotdetected.append(self.normnotdetected)
                self.finalefficiency.append(self.efficiency)
        
    

                self.bbilan=0
                self.bbilan=[]
        
                self.finalfinalflag=0
                self.finalfinalflag=[]
                self.finalfinalangdist=0
                self.finalfinalangdist=[]
                self.finalfinaldif=0                           
                self.finalfinaldif=[]
                self.finalmatchsrc=0
                self.finalmatchsrc=[]
                self.finalmatchsrcdist=0
                self.finalmatchsrcdist=[]
        
                self.detectedcoor=0
                self.detectedcoor=[]
                self.detectedcoordeg=0
                self.detectedcoordeg=[]
                self.textarray=0
                self.textarray=[]
                self.xytextarray=0
                self.xytextarray=[]
                self.gmaps=0
                self.gmaps=[]
                self.dmaps=0
                self.dmaps=[]
                self.detectedsrc=0
                self.detectedsrc=[]
                self.spurious=0
                self.spurious=[]
                self.notdetected=0
                self.notdetected=[]
                self.alldetection=0
                self.alldetection=[]
                self.normspurious=0
                self.normspurious=[]
                self.normspuriousdet=0
                self.normspuriousdet=[]
                self.normspuriouspix=0
                self.normspuriouspix=[]
                self.normnotdetected=0
                self.normnotdetected=[]
                self.efficiency=0
                self.efficiency=[]
                self.bilan=0
                self.bilan=[]
                self.detectedcoor=0
                self.detectedcoor=[]
        
                self.dmaps=0
                self.dmaps=[]
                self.detectedsrc=0
                self.detectedsrc=[]
                self.spurious=0
                self.spurious=[]
                self.notdetected=0
                self.notdetected=[]
                self.alldetection=0
                self.alldetection=[]
                self.finalflag=0
                self.finalflag=[]
                self.truedetection=0
                self.truedetection=[]
                self.truedetectioncoor=0
                self.truedetectioncoor=[]
                self.truedetectionflux=0
                self.truedetectionflux=[]
                self.stat=0
                self.stat=[]
                self.true=0
                self.true=[]

            self.finalfinalstat.append(self.finalstat)
            self.finalfinaltruedetectioncoor.append(self.finaltruedetectioncoor)
            self.finalfinaltruedetectionflux.append(self.finaltruedetectionflux)
            self.finalfinaltruedetection.append(self.finaltruedetection)
            self.finalfinalfinalflag.append(self.finalfinalflag)
            self.finalfinalbilan.append(self.finalbilan)
            self.finalfinaldetectedcoor.append(self.finaldetectedcoor)
            self.finalfinalgmaps.append(self.finalgmaps)
            self.finalfinaldmaps.append(self.finaldmaps)
            self.finalfinaldetectedsrc.append(self.finaldetectedsrc)
            self.finalfinalspurious.append(self.finalspurious)
            self.finalfinalnormspuriousdet.append(self.finalnormspuriousdet)
            self.finalfinalnormspuriouspix.append(self.finalnormspuriouspix)
            self.finalfinalnotdetected.append(self.finalnotdetected)
            self.finalfinalalldetection.append(self.finalalldetection)
            self.finalfinalnormspurious.append(self.finalnormspurious)
            self.finalfinalnormnotdetected.append(self.finalnormnotdetected)
            self.finalfinalefficiency.append(self.finalefficiency)  
            self.finalfinaltrue.append(self.finaltrue)


            self.finalfinalbilan.append(self.finalbilan)
   
            self.finalfinaldetectedsrc.append(self.finaldetectedsrc)
    
            self.finaltrue=0
            self.finaltrue=[]
            self.finalbbilan=0
            self.finalbbilan=[]
            self.finalfinalfinalflag=0
            self.finalfinalfinalflag=[]
            self.finalfinalfinalangdist=0
            self.finalfinalfinalangdist=[]
            self.finalfinalfinaldif=0                           
            self.finalfinalfinaldif=[]
            self.finalfinalmatchsrc=0
            self.finalfinalmatchsrc=[]
            self.finalfinalmatchsrcdist=0
            self.finalfinalmatchsrcdist=[] 
            self.finaldetectedcoor=0
            self.finaldetectedcoor=[]
            self.finaldetectedcoordeg=0
            self.finaldetectedcoordeg=[]
            self.finaltextarray=0
            self.finaltextarray=[]
            self.finaldmaps=0
            self.finaldmaps=[]
            self.finalgmaps=0
            self.finalgmaps=[]
            self.finaltxtfilesrcxy=0
            self.finaltxtfilesrcxy=[]
            self.finaldetectedsrc=0
            self.finaldetectedsrc=[]
            self.finalspurious=0
            self.finalspurious=[]
            self.finalnormspuriousdet=0
            self.finalnormspuriousdet=[]
            self.finalnormspuriouspix=0
            self.finalnormspuriouspix=[]
            self.finalnotdetected=0
            self.finalnotdetected=[]
            self.finalalldetection=0
            self.finalalldetection=[]
            self.finalnormspurious=0
            self.finalnormspurious=[]
            self.finalnormnotdetected=0
            self.finalnormnotdetected=[]
            self.finalefficiency=0
            self.finalefficiency=[]
            self.finalbilan=0
            self.finalbilan=[]
            self.finaltruedetection=0
            self.finaltruedetection=[]
            self.finalfinalflag=0
            self.finalfinalflag=[]
            self.finaltruedetectioncoor=0
            self.finaltruedetectioncoor=[]
            self.finaltruedetectionflux=0
            self.finaltruedetectionflux=[]
            self.finalstat=0
            self.finalstat=[]
          
                                
    def statistics(self):
        self.meanfinalfinalefficiencyeach=[]
        self.stdfinalfinalefficiencyeach=[]
        self.meanfinalfinalefficiency=[]
        self.stdfinalfinalefficiency=[]
        self.meanfinalfinalfinalefficiency=[]
        self.stdfinalfinalfinalefficiency=[]

        for noisex in range(0,np.shape(self.finalfinalefficiency)[0]):
            for i in range(0,np.shape(self.finalfinalefficiency)[1]):
                for j in range(0,np.shape(self.finalfinalefficiency)[2]):
                    self.meanfinalfinalefficiencyeach.append(np.mean(self.finalfinalefficiency[noisex][i][j]))
                    self.stdfinalfinalefficiencyeach.append(np.std(self.finalfinalefficiency[noisex][i][j]))
                self.meanfinalfinalefficiency.append(self.meanfinalfinalefficiencyeach)
                self.stdfinalfinalefficiency.append(self.stdfinalfinalefficiencyeach)
                self.meanfinalfinalefficiencyeach=0
                self.meanfinalfinalefficiencyeach=[]
                self.stdfinalfinalefficiencyeach=0
                self.stdfinalfinalefficiencyeach=[]

            self.meanfinalfinalfinalefficiency.append(self.meanfinalfinalefficiency)
            self.stdfinalfinalfinalefficiency.append(self.stdfinalfinalefficiency)
            self.meanfinalfinalefficiency=0
            self.meanfinalfinalefficiency=[]
            self.stdfinalfinalefficiency=0
            self.stdfinalfinalefficiency=[]
            ####    
        self.meanfinalfinalnormspuriouseach=[]
        self.stdfinalfinalnormspuriouseach=[]
        self.meanfinalfinalnormspurious=[]
        self.stdfinalfinalnormspurious=[]
        self.meanfinalfinalfinalnormspurious=[]
        self.stdfinalfinalfinalnormspurious=[]

        for noisex in range(0,np.shape(self.finalfinalnormspurious)[0]):
            for i in range(0,np.shape(self.finalfinalnormspurious)[1]):
                for j in range(0,np.shape(self.finalfinalnormspurious)[2]):
                    self.meanfinalfinalnormspuriouseach.append(np.mean(self.finalfinalnormspurious[noisex][i][j]))
                    self.stdfinalfinalnormspuriouseach.append(np.std(self.finalfinalnormspurious[noisex][i][j]))
        
                self.meanfinalfinalnormspurious.append(self.meanfinalfinalnormspuriouseach)
                self.stdfinalfinalnormspurious.append(self.stdfinalfinalnormspuriouseach)
                self.meanfinalfinalnormspuriouseach=0
                self.meanfinalfinalnormspuriouseach=[]
                self.stdfinalfinalnormspuriouseach=0
                self.stdfinalfinalnormspuriouseach=[]

            self.meanfinalfinalfinalnormspurious.append(self.meanfinalfinalnormspurious)
            self.stdfinalfinalfinalnormspurious.append(self.stdfinalfinalnormspurious)
            self.meanfinalfinalnormspurious=0
            self.meanfinalfinalnormspurious=[]
            self.stdfinalfinalnormspurious=0
            self.stdfinalfinalnormspurious=[]

        ####    
        self.meanfinalfinalnormspuriousdeteach=[]
        self.stdfinalfinalnormspuriousdeteach=[]
        self.meanfinalfinalnormspuriousdet=[]
        self.stdfinalfinalnormspuriousdet=[]
        self.meanfinalfinalfinalnormspuriousdet=[]
        self.stdfinalfinalfinalnormspuriousdet=[]

        for noisex in range(0,np.shape(self.finalfinalnormspurious)[0]):
            for i in range(0,np.shape(self.finalfinalnormspurious)[1]):
                for j in range(0,np.shape(self.finalfinalnormspurious)[2]):
                    self.meanfinalfinalnormspuriousdeteach.append(np.mean(self.finalfinalnormspuriousdet[noisex][i][j]))
                    self.stdfinalfinalnormspuriousdeteach.append(np.std(self.finalfinalnormspuriousdet[noisex][i][j]))
        
                self. meanfinalfinalnormspuriousdet.append(self.meanfinalfinalnormspuriousdeteach)
                self.stdfinalfinalnormspuriousdet.append(self.stdfinalfinalnormspuriousdeteach)
                self.meanfinalfinalnormspuriousdeteach=0
                self.meanfinalfinalnormspuriousdeteach=[]
                self.stdfinalfinalnormspuriousdeteach=0
                self.stdfinalfinalnormspuriousdeteach=[]

            self.meanfinalfinalfinalnormspuriousdet.append(self.meanfinalfinalnormspuriousdet)
            self.stdfinalfinalfinalnormspuriousdet.append(self.stdfinalfinalnormspuriousdet)
            self.meanfinalfinalnormspuriousdet=0
            self.meanfinalfinalnormspuriousdet=[]
            self.stdfinalfinalnormspuriousdet=0
            self.stdfinalfinalnormspuriousdet=[]

        ####    
        self.meanfinalfinalnormspuriouspixeach=[]
        self.stdfinalfinalnormspuriouspixeach=[]
        self.meanfinalfinalnormspuriouspix=[]
        self.stdfinalfinalnormspuriouspix=[]
        self.meanfinalfinalfinalnormspuriouspix=[]
        self.stdfinalfinalfinalnormspuriouspix=[]

        for noisex in range(0,np.shape(self.finalfinalnormspurious)[0]):
            for i in range(0,np.shape(self.finalfinalnormspurious)[1]):
                for j in range(0,np.shape(self.finalfinalnormspurious)[2]):
                    self.meanfinalfinalnormspuriouspixeach.append(np.mean(self.finalfinalnormspuriouspix[noisex][i][j]))
                    self.stdfinalfinalnormspuriouspixeach.append(np.std(self.finalfinalnormspuriouspix[noisex][i][j]))
        
                self.meanfinalfinalnormspuriouspix.append(self.meanfinalfinalnormspuriouspixeach)
                self.stdfinalfinalnormspuriouspix.append(self.stdfinalfinalnormspuriouspixeach)
                self.meanfinalfinalnormspuriouspixeach=0
                self.meanfinalfinalnormspuriouspixeach=[]
                self.stdfinalfinalnormspuriouspixeach=0
                self.stdfinalfinalnormspuriouspixeach=[]

            self.meanfinalfinalfinalnormspuriouspix.append(self.meanfinalfinalnormspuriouspix)
            self.stdfinalfinalfinalnormspuriouspix.append(self.stdfinalfinalnormspuriouspix)
            self.meanfinalfinalnormspuriouspix=0
            self.meanfinalfinalnormspuriouspix=[]
            self.stdfinalfinalnormspuriouspix=0
            self.stdfinalfinalnormspuriouspix=[]

        ##################
        self.meanfinalfinalstateach=[]
        self.stdfinalfinalstateach=[]
        self.meanfinalfinalstat=[]
        self.stdfinalfinalstat=[]
        self.meanfinalfinalfinalstat=[]
        self.stdfinalfinalfinalstat=[]

        for noisex in range(0,np.shape(self.finalfinalstat)[0]):
            for i in range(0,np.shape(self.finalfinalstat)[1]):
                for j in range(0,np.shape(self.finalfinalstat)[2]):
                    self.meanfinalfinalstateach.append(np.mean(self.finalfinalstat[noisex][i][j]))
                    self.stdfinalfinalstateach.append(np.std(self.finalfinalstat[noisex][i][j]))
        
                self.meanfinalfinalstat.append(self.meanfinalfinalstateach)
                self.stdfinalfinalstat.append(self.stdfinalfinalstateach)
                self.meanfinalfinalstateach=0
                self.meanfinalfinalstateach=[]
                self.stdfinalfinalstateach=0
                self.stdfinalfinalstateach=[]

            self.meanfinalfinalfinalstat.append(self.meanfinalfinalstat)
            self.stdfinalfinalfinalstat.append(self.stdfinalfinalstat)
            self.meanfinalfinalstat=0
            self.meanfinalfinalstat=[]
            self.stdfinalfinalstat=0
            self.stdfinalfinalstat=[]


        #creating  a set of flux to plot a scatter plot
        self.fluxeacheach=[]
        self.fluxeach=[]
        self.flux=[]
        self.finalflux=[]
        self.finalfinalflux=[]
        self.srcnum=10

        for sigx in range(0,len(self.signoise)):
            for nx in range(0,len(self.n)):
                print('for a given false detection rate n={}'.format(self.n[nx]))
                for i in range(0,len(self.jandec)):
                    for j in range(0,self.iteration):
                        for l in range(0,len(self.finalfinaltruedetectionflux[sigx][nx][i][j])):
                            self.fluxeacheach.append(self.jandec[i])
                self.fluxeach.append(self.fluxeacheach)
                self.fluxeacheach=0
                self.fluxeacheach=[]
            self.flux.append(self.fluxeach)
            self.fluxeach=0
            self.fluxeach=[]
                
            
    


        for sigx in range(0,np.shape(self.finalfinaltrue)[0]):
            for nx in range(0,np.shape(self.finalfinaltrue)[1]):
                for i in range(0,len(self.finalfinaltrue[sigx][nx])):
                    if len(self.finalfinaltrue[sigx][nx][i])>1:
                        self.finalfinaltrue[sigx][nx][i]=np.mean(self.finalfinaltrue[sigx][nx][i])


    def plot(self):
        plt.figure()
        for noisex in range(0,np.shape(self.meanfinalfinalfinalefficiency)[0]): 
            plt.subplot(310+(noisex+1))
            for i in range(0,np.shape(self.meanfinalfinalfinalefficiency)[1]):
                plt.errorbar(self.jandec,self.meanfinalfinalfinalefficiency[noisex][i],yerr=self.stdfinalfinalfinalefficiency[noisex][i],label='eff n={}'.format(self.n[i]))
                plt.xlabel('Sources flux [Jansky]')
                plt.ylabel('efficiency')
                plt.legend(loc='lower right')
                plt.title('efficiency as a function of sources flux for two false detection rate and noise={} Kelvin'.format(self.signoise[noisex]))
        plt.show()
            
        


        plt.figure()
        for noisex in range(0,np.shape(self.meanfinalfinalfinalefficiency)[0]): 
            plt.subplot(310+(noisex+1))
            for i in range(0,np.shape(self.meanfinalfinalfinalefficiency)[1]):
                plt.errorbar(self.jandec,self.meanfinalfinalfinalnormspurious[noisex][i],yerr=self.stdfinalfinalfinalnormspurious[noisex][i],label='spur n={}'.format(self.n[i]))
                plt.xlabel('Sources flux [Jansky]')
                plt.ylabel('spurious detections')
                #plt.plot(jandec,normnotdetected)
                plt.legend(loc='lower right')
                plt.title('spurious detections as a function of sources flux for two false detection rate and noise={} Kelvin'.format(self.signoise[noisex]))
        plt.show()

        # many n many noises ( fitting)
        nx2=1
        plt.figure()
        for noisex in range(0,np.shape(self.meanfinalfinalfinalefficiency)[0]): 
            for nx in range(0,np.shape(self.meanfinalfinalfinalefficiency)[1]):
                plt.errorbar(self.jandec,self.meanfinalfinalfinalefficiency[noisex][nx2],yerr=self.stdfinalfinalfinalefficiency[noisex][nx2],label='eff n={} noise={}'.format(self.n[nx2],self.signoise[noisex]))
                plt.xlabel('Sources flux')
                plt.ylabel('efficiency')
                plt.legend(loc='lower right')
                plt.title('efficiency in function of sources flux for many n  and many noises')
        plt.show()    


        # one n many noises ( fitting)
        plt.figure()
        for noisex in range(0,np.shape(self.meanfinalfinalfinalefficiency)[0]): 
            plt.errorbar(self.jandec,self.meanfinalfinalfinalefficiency[noisex][1],yerr=self.stdfinalfinalfinalefficiency[noisex][1],label='eff n={} noise={}'.format(self.n[1],self.signoise[noisex]))
            plt.xlabel('Sources flux')
            plt.ylabel('efficiency')
            plt.legend(loc='lower right')
           # plt.title('efficiency in function of sources flux for n={} and many noises'.format(self.n[]))
        plt.show()    

        
        plt.figure()
        plt.subplot(211)
        for i in range(0,np.shape(self.meanfinalfinalfinalefficiency)[1]):
            plt.errorbar(self.jandec,self.meanfinalfinalfinalefficiency[noisex][i],yerr=self.stdfinalfinalfinalefficiency[noisex][i],label='eff n={}'.format(self.n[i]))
            plt.errorbar(self.jandec,self.meanfinalfinalfinalnormspurious[noisex][i],yerr=self.stdfinalfinalfinalnormspurious[noisex][i],label='spur n={}'.format(self.n[i]))
            plt.xlabel('Sources flux [Jansky]')
            plt.ylabel('efficiency')
            #plt.plot(jandec,normnotdetected)
        plt.legend(loc='lower right')
        plt.title('efficiency as a function of sources flux for two false detection rate and noise={} Kelvin'.format(self.signoise[noisex]))
        plt.subplot(212)
        for noisex in range(0,np.shape(self.meanfinalfinalfinalefficiency)[0]):  
            plt.errorbar(self.jandec,self.meanfinalfinalfinalefficiency[noisex][1],yerr=self.stdfinalfinalfinalefficiency[noisex][1],label='eff n={} noise={}'.format(self.n[1],self.signoise[noisex]))
            #plt.errorbar(jandec,meanfinalfinalfinalnormspurious[noisex][1],yerr=stdfinalfinalfinalnormspurious[noisex][1],label='spur n={} noise={}'.format(n[1],signoise[noisex]))
            plt.xlabel('Sources flux')
            plt.ylabel('efficiency')
            #plt.plot(jandec,normnotdetected)
            plt.legend(loc='lower right')
            plt.title('efficiency in function of sources flux for n=1.2 and many noises'.format(self.signoise[noisex]))
                
        plt.show()

        ######################


        plt.figure()
        for noisex in range(0,np.shape(self.meanfinalfinalfinalefficiency)[0]): 
            plt.subplot(310+(noisex+1))
            ax=plt.subplot(310+(noisex+1))
            for i in range(0,np.shape(self.meanfinalfinalfinalefficiency)[1]):
                y = self.meanfinalfinalfinalefficiency[noisex][i]
                x = self.jandec
                ax.errorbar(self.jandec,self.meanfinalfinalfinalefficiency[noisex][i],yerr=self.stdfinalfinalfinalefficiency[noisex][i],label='n={}'.format(self.n[i]))
                ax2=ax.twiny()
                ax2.set_xticks(self.jandec/self.signoiseJden[noisex])
                ax2.set_xlabel('Flux density[sigma]')
                ax3=ax.twiny()
                ax3.spines['top'].set_position(('outward', 40))
                ax3.set_xticks(self.jandecmKamp)
                ax3.set_xlabel('amplitude [milliKelvin]')
                ax.set_xlabel('Flux density [Jansky]')
                ax.set_ylabel('efficiency')
                ax.legend(loc='lower right')
                ax.set_title('efficiency as a function of sources flux density and noise={} Kelvin'.format(self.signoise[noisex]))
            
        

        fig, ax = plt.subplots(constrained_layout=True)
        for i in range(0,np.shape(self.meanfinalfinalfinalefficiency)[1]):
            y = self.meanfinalfinalfinalefficiency[noisex][i]
            x = self.jandec
            ax.errorbar(self.jandec,self.meanfinalfinalfinalefficiency[noisex][i],yerr=self.stdfinalfinalfinalefficiency[noisex][i],label='n={}'.format(self.n[i]))
            ax2=ax.twiny()
            ax2.set_xticks(self.jandec/self.signoiseJden[noisex])
            ax2.set_xlabel('Flux density[sigma]')
            ax3=ax.twiny()
            ax3.spines['top'].set_position(('outward', 40))
            ax3.set_xticks(self.jandecmKamp)
            ax3.set_xlabel('amplitude [milliKelvin]')
            ax.set_xlabel('Flux density [Jansky]')
            ax.set_ylabel('efficiency')
            ax.legend(loc='lower right')
            ax.set_title('efficiency as a function of sources flux density and noise={} Kelvin'.format(self.signoise[noisex]))
        #plt.show()

        fig, ax = plt.subplots(constrained_layout=True)
        for i in range(0,np.shape(self.meanfinalfinalfinalefficiency)[1]):
            y = self.meanfinalfinalfinalefficiency[noisex][i]
            x = self.jandec
            ax.errorbar(self.jandec,self.meanfinalfinalfinalnormspurious[noisex][i],yerr=self.stdfinalfinalfinalnormspurious[noisex][i],label='n={}'.format(self.n[i]))
            ax2=ax.twiny()
            ax2.set_xticks(self.jandec/self.signoiseJden[noisex])
            ax2.set_xlabel('Flux density[sigma]')
            ax3=ax.twiny()
            ax3.spines['top'].set_position(('outward', 40))
            ax3.set_xticks(self.jandecmKamp)
            ax3.set_xlabel('amplitude [milliKelvin]')
            ax.set_xlabel('Flux density [Jansky]')
            ax.set_ylabel('[spurious detection / generated sources]')
            ax.legend(loc='lower right')
            ax.set_title('spurious detection as a function of sources flux density and noise={} Kelvin'.format(self.signoise[noisex]))
        #plt.show()

        fig, ax = plt.subplots(constrained_layout=True)
        for i in range(0,np.shape(self.meanfinalfinalfinalefficiency)[1]):
            y = self.meanfinalfinalfinalefficiency[noisex][i]
            x = self.jandec
            ax.errorbar(self.jandec,self.meanfinalfinalfinalnormspuriousdet[noisex][i],yerr=self.stdfinalfinalfinalnormspuriousdet[noisex][i],label='n={}'.format(self.n[i]))
            ax2=ax.twiny()
            ax2.set_xticks(self.jandec/self.signoiseJden[noisex])
            ax2.set_xlabel('Flux density[sigma]')
            ax3=ax.twiny()
            ax3.spines['top'].set_position(('outward', 40))
            ax3.set_xticks(self.jandecmKamp)
            ax3.set_xlabel('amplitude [milliKelvin]')
            #ax3.errorbar(xampK,meanfinalfinalfinalefficiency[noisex][i],yerr=stdfinalfinalfinalefficiency[noisex][i],label='eff n={}'.format(n[i]))
            ax.set_xlabel('Flux density [Jansky]')
            ax.set_ylabel('[spurious detection / detected sources]')
            ax.legend(loc='lower right')
            ax.set_title('spurious detection as a function of sources flux density and noise={} Kelvin'.format(self.signoise[noisex]))
        #plt.show()

        fig, ax = plt.subplots(constrained_layout=True)
        for i in range(0,np.shape(self.meanfinalfinalfinalefficiency)[1]):
            y = self.meanfinalfinalfinalefficiency[noisex][i]
            x = self.jandec
            ax.errorbar(self.jandec,self.meanfinalfinalfinalnormspuriouspix[noisex][i],yerr=self.stdfinalfinalfinalnormspuriouspix[noisex][i],label='n={}'.format(self.n[i]))
            ax2=ax.twiny()
            ax2.set_xticks(x/self.signoiseJden[noisex])
            ax2.set_xlabel('Flux density[sigma]')
            ax3=ax.twiny()
            ax3.spines['top'].set_position(('outward', 40))
            ax3.set_xticks(self.jandecmKamp)
            ax3.set_xlabel('amplitude [milliKelvin]')
            #ax3.errorbar(xampK,meanfinalfinalfinalefficiency[noisex][i],yerr=stdfinalfinalfinalefficiency[noisex][i],label='eff n={}'.format(n[i]))
            ax.set_xlabel('Flux density [Jansky]')
            ax.set_ylabel('spurious detection / map pixels  ')
            ax.legend(loc='lower right')
            ax.set_title('spurious detection,  noise={} Kelvin'.format(self.signoise[noisex]))
        plt.show()

        fig, ax = plt.subplots(constrained_layout=True)
        for i in range(0,np.shape(self.meanfinalfinalfinalefficiency)[1]):
            y=np.reshape(self.finalfinaltrue[noisex][i],(len(self.finalfinaltrue[noisex][i])))
            #x=np.arange(0,len(finalfinaltrue[noisex][i])*(np.max(jandec)/len(finalfinaltrue[noisex][i])),(np.max(jandec)/len(finalfinaltrue[noisex][i])))
            x=self.flux[noisex][i]
            plt.scatter(x,y,label='n={}'.format(self.n[i]))
            plt.plot(y,y)
            ax.set_xlabel('generated flux')
            ax.set_ylabel('detected flux')
            ax.legend(loc='lower right')
        plt.show()

        for i in range(0,np.shape(self.meanfinalfinalfinalefficiency)[1]):
            fig, ax = plt.subplots(constrained_layout=True)
            y=np.reshape(self.finalfinaltrue[noisex][i],(len(self.finalfinaltrue[noisex][i])))
            #x=np.arange(0,len(finalfinaltrue[noisex][i])*(np.max(jandec)/len(finalfinaltrue[noisex][i])),(np.max(jandec)/len(finalfinaltrue[noisex][i])))
            x=self.flux[noisex][i]
            plt.scatter(x,y,label='n={}'.format(self.n[i]))
            plt.plot(y,y)
            ax.set_xlabel('generated flux')
            ax.set_ylabel('detected flux')
            ax.legend(loc='lower right')
            plt.show()



        fig, ax = plt.subplots(constrained_layout=True)
        for i in range(0,np.shape(self.meanfinalfinalfinalefficiency)[1]):
            y = self.meanfinalfinalfinalefficiency[noisex][i]
            x = self.jandec
            #ax.scatter(jandec,meanfinalfinalfinalstat[noisex][i],label='n={}'.format(n[i]))
            ax.plot(x,x)
            ax.errorbar(self.jandec,self.meanfinalfinalfinalstat[noisex][i],yerr=self.stdfinalfinalfinalstat[noisex][i],ls='none',label='n={}'.format(self.n[i]),fmt='o')
            #ax.errorbar(jandec,meanfinalfinalfinalstat[noisex][i],yerr=stdfinalfinalfinalstat[noisex][i],label='n={}'.format(n[i]),fmt='o')
            ax2=ax.twiny()
            ax2.set_xticks(x/self.signoiseJden[noisex])
            ax2.set_xlabel('Flux density[sigma]')
            ax3=ax.twiny()
            ax3.spines['top'].set_position(('outward', 40))
            ax3.set_xticks(self.jandecmKamp)
            ax3.set_xlabel('amplitude [milliKelvin]')
            #ax3.errorbar(xampK,meanfinalfinalfinalefficiency[noisex][i],yerr=stdfinalfinalfinalefficiency[noisex][i],label='eff n={}'.format(n[i]))
            ax.set_xlabel('Flux density [Jansky]')
            ax.set_ylabel('Fluxes median')
            ax.legend(loc='lower right')
            ax.set_title('Fluxes fitting as a function of sources flux density and noise={} Kelvin'.format(self.signoise[noisex]))
            plt.show()


        plt.figure()
        for noisex in range(0,len(self.signoise)):
            ax = plt.subplot(310+(noisex+1))
            #ax = plt.subplot()
            plt.subplots_adjust(hspace=0.6)
            for i in range(0,np.shape(self.meanfinalfinalfinalefficiency)[1]):
                y = self.meanfinalfinalfinalefficiency[noisex][i]
                x = self.jandec
                ax.errorbar(self.jandec,self.meanfinalfinalfinalefficiency[noisex][i],yerr=self.stdfinalfinalfinalnormspurious[noisex][i],label='spur n={}'.format(self.n[i]))
            ax.set_xlabel('Flux [Jansky]')
            ax.set_ylabel('efficiency')
            ax.legend(loc='lower right')
            ax.set_title('efficiency as a function of sources flux for noise={} Kelvin'.format(self.signoise[noisex]))
        #secax = ax.secondary_xaxis('top', functions=(deg2rad, rad2deg))
        plt.show()



        plt.figure()
        for noisex in range(0,len(self.signoise)):
            ax = plt.subplot(310+(noisex+1))
            #ax = plt.subplot()
            plt.subplots_adjust(hspace=0.6)
            for i in range(0,np.shape(self.meanfinalfinalfinalefficiency)[1]):
                y = self.meanfinalfinalfinalefficiency[noisex][i]
                x = self.jandec
                ax.errorbar(self.jandec,self.meanfinalfinalfinalnormspurious[noisex][i],yerr=self.stdfinalfinalfinalnormspurious[noisex][i],label='spur n={}'.format(self.n[i]))
                ax2=ax.twiny()
                ax2.set_xticks(x/self.signoiseJden[noisex])
            ax.set_xlabel('Flux [Jansky]')
            ax.set_ylabel('spurious detections')
            ax.legend(loc='lower right')
            ax.set_title('normalized spurious detection number as a function of sources flux for  noise={} Kelvin'.format(self.signoise[noisex]))
    #secax = ax.secondary_xaxis('top', functions=(deg2rad, rad2deg))
            plt.show()

    def plot_all(self):

        # foreach noise a plot with many n
        for noisex in range(0,np.shape(self.meanfinalfinalfinalefficiency)[0]):
            fig, ax = plt.subplots(constrained_layout=True)
            for i in range(0,np.shape(self.meanfinalfinalfinalefficiency)[1]):
                ax.errorbar(self.jandec,self.meanfinalfinalfinalefficiency[noisex][i],yerr=self.stdfinalfinalfinalefficiency[noisex][i],label='n={}'.format(self.n[i]))
                ax2=ax.twiny()
                ax2.set_xticks(self.jandec/self.signoiseJden[noisex])
                ax2.set_xlabel('Flux density[sigma]')
                ax3=ax.twiny()
                ax3.spines['top'].set_position(('outward', 40))
                ax3.set_xticks(self.jandecmKamp)
                ax3.set_xlabel('amplitude [milliKelvin]')
                ax.set_xlabel('Flux density [Jansky]')
                ax.set_ylabel('efficiency')
                ax.legend(loc='lower right')
                ax.set_title('efficiency, noise={} Kelvin'.format(self.signoise[noisex]))
        #plt.show()


            fig, ax = plt.subplots(constrained_layout=True)
            for i in range(0,np.shape(self.meanfinalfinalfinalefficiency)[1]):
                ax.errorbar(self.jandec,self.meanfinalfinalfinalnormspurious[noisex][i],yerr=self.stdfinalfinalfinalnormspurious[noisex][i],label='n={}'.format(self.n[i]))
                ax2=ax.twiny()
                ax2.set_xticks(self.jandec/self.signoiseJden[noisex])
                ax2.set_xlabel('Flux density[sigma]')
                ax3=ax.twiny()
                ax3.spines['top'].set_position(('outward', 40))
                ax3.set_xticks(self.jandecmKamp)
                ax3.set_xlabel('amplitude [milliKelvin]')
                ax.set_xlabel('Flux density [Jansky]')
                ax.set_ylabel('[spurious detection / generated sources]')
                ax.legend(loc='lower right')
                ax.set_title('spurious detection, noise={} Kelvin'.format(self.signoise[noisex]))
        #plt.show()

            fig, ax = plt.subplots(constrained_layout=True)
            for i in range(0,np.shape(self.meanfinalfinalfinalefficiency)[1]):
                ax.errorbar(self.jandec,self.meanfinalfinalfinalnormspuriousdet[noisex][i],yerr=self.stdfinalfinalfinalnormspuriousdet[noisex][i],label='n={}'.format(self.n[i]))
                ax2=ax.twiny()
                ax2.set_xticks(self.jandec/self.signoiseJden[noisex])
                ax2.set_xlabel('Flux density[sigma]')
                ax3=ax.twiny()
                ax3.spines['top'].set_position(('outward', 40))
                ax3.set_xticks(self.jandecmKamp)
                ax3.set_xlabel('amplitude [milliKelvin]')
                ax.set_xlabel('Flux density [Jansky]')
                ax.set_ylabel('[spurious detection / detected sources]')
                ax.legend(loc='lower right')
                ax.set_title('spurious detection, noise={} Kelvin'.format(self.signoise[noisex]))
        #plt.show()

            fig, ax = plt.subplots(constrained_layout=True)
            for i in range(0,np.shape(self.meanfinalfinalfinalefficiency)[1]):
                ax.errorbar(self.jandec,self.meanfinalfinalfinalnormspuriouspix[noisex][i],yerr=self.stdfinalfinalfinalnormspuriouspix[noisex][i],label='n={}'.format(self.n[i]))
                ax2=ax.twiny()
                ax2.set_xticks(self.jandec/self.signoiseJden[noisex])
                ax2.set_xlabel('Flux density[sigma]')
                ax3=ax.twiny()
                ax3.spines['top'].set_position(('outward', 40))
                ax3.set_xticks(self.jandecmKamp)
                ax3.set_xlabel('amplitude [milliKelvin]')
                ax.set_xlabel('Flux density [Jansky]')
                ax.set_ylabel('spurious detection / map pixels  ')
                ax.legend(loc='lower right')
                ax.set_title('spurious detection, noise={} Kelvin'.format(self.signoise[noisex]))



            fig, ax = plt.subplots(constrained_layout=True)
            for i in range(0,np.shape(self.meanfinalfinalfinalefficiency)[1]):
                x=self.jandec
                ax.plot(x,x)
                ax.errorbar(self.jandec,self.meanfinalfinalfinalstat[noisex][i],yerr=self.stdfinalfinalfinalstat[noisex][i],ls='none',label='n={}'.format(self.n[i]),fmt='o')
                ax2=ax.twiny()
                ax2.set_xticks(self.jandec/self.signoiseJden[noisex])
                ax2.set_xlabel('Flux density[sigma]')
                ax3=ax.twiny()
                ax3.spines['top'].set_position(('outward', 40))
                ax3.set_xticks(self.jandecmKamp)
                ax3.set_xlabel('amplitude [milliKelvin]')
                ax.set_xlabel('Flux density [Jansky]')
                ax.set_ylabel('Fluxes median')
                ax.legend(loc='lower right')
                ax.set_title('Fluxes fitting as a function of sources flux density and noise={} Kelvin'.format(self.signoise[noisex]))
            

            fig, ax = plt.subplots(constrained_layout=True)
            for i in range(0,np.shape(self.meanfinalfinalfinalefficiency)[1]):
                y=np.reshape(self.finalfinaltrue[noisex][i],(len(self.finalfinaltrue[noisex][i])))
                x=self.flux[noisex][i]
                plt.scatter(x,y,label='n={}'.format(self.n[i]))
                plt.plot(y,y)
                ax.set_xlabel('generated flux')
                ax.set_ylabel('detected flux')
                ax.legend(loc='lower right')
                ax.set_title('detected flux vs generated flux, noise={} Kelvin'.format(self.signoise[noisex]))
       
           
    def plot_eachn(self):

        # one n many noises ( fitting)
        
        plt.figure()
        for noisex in range(0,np.shape(self.meanfinalfinalfinalefficiency)[0]): 
            plt.errorbar(self.jandec,self.meanfinalfinalfinalefficiency[noisex][self.nx2],yerr=self.stdfinalfinalfinalefficiency[noisex][self.nx2],label='noise={}'.format(self.signoise[noisex]))
            plt.xlabel('Sources flux [Jansky]')
            plt.ylabel('efficiency')
            plt.legend(loc='lower right')
            plt.title('efficiency, for n={} and many noises'.format(self.n[self.nx2]))
           
        plt.figure()
        for noisex in range(0,np.shape(self.meanfinalfinalfinalefficiency)[0]): 
            plt.errorbar(self.jandec,self.meanfinalfinalfinalnormspurious[noisex][self.nx2],yerr=self.stdfinalfinalfinalnormspurious[noisex][self.nx2],label='noise={}'.format(self.signoise[noisex]))
            plt.xlabel('Sources flux [Jansky]')
            plt.ylabel('[spurious detection / generated sources]')
            plt.legend(loc='lower right')
            plt.title('spurious detections, for n={} and many noises'.format(self.n[self.nx2]))

        plt.figure()
        for noisex in range(0,np.shape(self.meanfinalfinalfinalefficiency)[0]): 
            plt.errorbar(self.jandec,self.meanfinalfinalfinalnormspuriousdet[noisex][self.nx2],yerr=self.stdfinalfinalfinalnormspuriousdet[noisex][self.nx2],label='noise={}'.format(self.signoise[noisex]))
            plt.xlabel('Sources flux [Jansky]')
            plt.ylabel('[spurious detection / detected sources]')
            plt.legend(loc='lower right')
            plt.title('spurious detections, for n={} and many noises'.format(self.n[self.nx2]))

        plt.figure()
        for noisex in range(0,np.shape(self.meanfinalfinalfinalefficiency)[0]): 
            plt.errorbar(self.jandec,self.meanfinalfinalfinalnormspuriouspix[noisex][self.nx2],yerr=self.stdfinalfinalfinalnormspuriouspix[noisex][self.nx2],label='noise={}'.format(self.signoise[noisex]))
            plt.xlabel('Sources flux')
            plt.ylabel('spurious detection / map pixels')
            plt.legend(loc='lower right')
            plt.title('spurious detection, for n={} and many noises'.format(self.n[self.nx2]))



        fig, ax = plt.subplots(constrained_layout=True)
        for noisex in range(0,np.shape(self.meanfinalfinalfinalefficiency)[0]):
            x=self.jandec
            ax.plot(x,x)
            ax.errorbar(self.jandec,self.meanfinalfinalfinalstat[noisex][self.nx2],yerr=self.stdfinalfinalfinalstat[noisex][self.nx2],ls='none',label='noise={}'.format(self.signoise[noisex]),fmt='o')
            ax2=ax.twiny()
            ax2.set_xticks(self.jandec/self.signoiseJden[noisex])
            ax2.set_xlabel('Flux density[sigma]')
            ax3=ax.twiny()
            ax3.spines['top'].set_position(('outward', 40))
            ax3.set_xticks(self.jandecmKamp)
            ax3.set_xlabel('amplitude [milliKelvin]')
            ax.set_xlabel('Flux density [Jansky]')
            ax.set_ylabel('Fluxes median')
            ax.legend(loc='lower right')
            ax.set_title('Fluxes fitting as a function of sources flux density and n={} Kelvin'.format(self.n[self.nx2]))
            

           

        fig, ax = plt.subplots(constrained_layout=True)
        for noisex in range(0,np.shape(self.meanfinalfinalfinalefficiency)[0]): 
            y=np.reshape(self.finalfinaltrue[noisex][self.nx2],(len(self.finalfinaltrue[noisex][self.nx2])))
            x=self.flux[noisex][self.nx2]
            plt.scatter(x,y,label='noise={}'.format(self.signoise[noisex]))
            plt.plot(y,y)
            ax.set_xlabel('generated flux')
            ax.set_ylabel('detected flux')
            ax.legend(loc='lower right')
            ax.set_title('detected flux vs generated flux for many noises, n={} Kelvin'.format(self.n[self.nx2]))


            


