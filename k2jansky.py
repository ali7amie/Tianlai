import numpy as np
import scipy.constants


def k2jansky(use,freq,max_baseline,beam_surface,density_k):
    """ This function convert flux density in kelvin to flux density in jansky, if use=='freq' the used beam correspond to interferometer resolution lambda/D, then you must enter value for fréquency [Hz] and max_baseline [m] and omit beam_surface, if use=='surface' the beam_surface is entered manually and freq and max_baseline are omitted"""

    c=3*10**8
    lambdaa=c/freq
    if use=='freq':
        density_jansky=(2*scipy.constants.Boltzmann*density_k*(lambdaa/max_baseline)**2)/lambdaa**2
        density_jansky=density_jansky/10**(-26)
        return density_jansky

    if use=='surface':
        density_jansky=(2*scipy.constants.Boltzmann*density_k*beam_surface)/lambdaa**2
        density_jansky=density_jansky/10**(-26)
        return density_jansky
   
    #beam_surf=6*10**(-5)
    #beam_surf=(lambdaa/maxbaseline)**2#33(58 arcmin**2 in rad) boltzmann (J.k-1)
    #reso=12
    #beam_surf=(reso*(1/60)*(math.pi/180))**2
    #return S

#________________________________________________________________________________________________

def jansky2k(use,freq,max_baseline,beam_surface,density_jansky):
    """ This function convert flux density in jansky to flux density in kelvin, if use=='freq' the used beam correspond to interferometer resolution lambda/D, then you must enter value for fréquency [Hz] and max_baseline [m] and omit beam_surface, if use=='surface' the beam_surface is entered manually and freq and max_baseline are omitted"""

    c=3*10**8
    lambdaa=c/freq
    if use=='freq':
        density_k=density_jansky*10**(-26)*lambdaa**2/(2*scipy.constants.Boltzmann*(lambdaa/max_baseline)**2)
        return density_k

    if use=='surface':
        density_k=density_jansky*10**(-26)*lambdaa**2/(2*scipy.constants.Boltzmann*beam_surface)
        return density_k

        
