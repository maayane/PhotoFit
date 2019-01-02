

import astropy
from astropy import constants as const
import math
from . import distances_conversions
from . import extinction
import numpy as np
import pdb
import pylab

#def planck(wav, T):
#    a=2*6.626070040e-34*(3e8)**2
#    b=6.626070040e-34*(3e8)/(wav*T*1.38064852e-23)
#    #a = 2*const.h.value*const.c.value**2
#    #b = const.h.value*const.c.value/(wav*const.k_B.value*T) #convert into cgs
#    intensity = a/ ( (wav**5) * (np.exp(b) - 1.0) )
#    return intensity#*u.J/(u.s*u.m*u.m*u.m).cgs

#def planck_cgs(wav,T):
#    a=2*const.h.cgs.value*(const.c.cgs.value)**2
#    b=const.h.cgs.value*const.c.cgs.value/(wav*const.k_B.cgs.value*T)
#    intensity=a/( (wav**5) * (np.exp(b) - 1.0) )
#    return intensity

#def RayleighJeans(wav, T):
#    a = 2*c*k*T
#    intensity = a/wav**4
#    return intensity

#def Wien(wav, T):
#    a = 2*h*c**2
#    b=h*c/(wav*k*T)
#    intensity = (a/wav**5)*np.exp(-b)
#    return intensity

def black_body_flux_density(Temp,wavelength,type=None,verbose=False,distance_pc=None,Radius=None,Ebv=None,R_ext=None,redshift=0,plot=False):
    """Description: Given a temperature, calculates a black body flux density B_lambda.
    If a radius anda distance are given, calculate the apparent flux density (R/d)^2*B_lambda
    Input  :- Temperature [K]
            - numpy array of wavelengths [m], tipically np.linspace(1e-10,1e-6,num=1000)
            - type of formula:
                'P' Planck
                'RJ' Rayleigh-Jeans approximation
            - Radius (optionnal) in solar radius
            - distance (optionnal) in pc
            - Ebv: (optionnal, default is none) extinction to APPLY to the theoretical bb spectrum
            - redshift: (optionnal, default is none) z to apply to the theoretical spectrum with
    Output :array of numpy.arrays [spectrum_cgs,spectrum_Hz,spectrum_A,spectrum_mJy,spectrum_phot] CAREFULLL! confusing between spectrum_cgs and spectrum_A has caused so much arm in the past!
            - spectrum_cgs: wavelength [m], Emittance (flux density) in erg/sec/cm^2/cm(lambda)
            - spectrum_Hz: wavelength [m], Emittance in erg/sec/cm^2/Hz
            - spectrum_A: wavelength [m], Emittance in erg/sec/cm^2/Ang (lambda), 1e-8*Emittance (flux density) in erg/sec/cm^2/cm(lambda)
            - spectrum_mjy: wavelength [m], Emittance [mJy]
            - spectrum_phot: wavelength [m], number of photons [photons/sec/cm^2/Ang (lambda)]
    Tested : ?
         By : Maayane T. Soumagnac Nov 2016
        URL :
    Example:[E_cgs, E_Hz, E_A,Emjy, E_phot] = black_body_models.black_body_models(3000, wavelengths, 'P')
    Reliable:  """
    #if Ebv==0.:
    #    Ebv=None
    h_cgs=const.h.cgs.value
    c_cgs=const.c.cgs.value
    kB_cgs=const.k_B.cgs.value
    h_USI=const.h.value
    c_USI=const.c.value
    kB_USI=const.k_B.value
    wavelength_in_cm=wavelength*1e2 # wavelength in cgs
    wavelength_in_cm = wavelength_in_cm.astype(float)
    nu=c_cgs/wavelength_in_cm #frequency in s (because c is in cm/s and wavlength in cm)
    if (Radius!=None and distance_pc!=None):
        #print('pouet')
        #print 'the radius and distance were specified: I am calculating the apparent flux density'
        #pdb.set_trace()
        R_pc=distances_conversions.solar_radius_to_pc(Radius)
        coeff=(R_pc/distance_pc)**2
        #print coeff
        #pdb.set_trace()
    else:
        if verbose==True:
            print('the radius or distance or both were not specified')
        coeff=1.
        #print coeff
    #pdb.set_trace()
    if type.lower() in (None,'p'):
        if verbose==True:
            print('formula used for black body: Planck')
        b_cgs=h_cgs*c_cgs/(wavelength_in_cm*kB_cgs*Temp)
        #b_USI=h_USI*c_USI/(wavelength*kB_USI*Temp)

        #print 'b_USI is',b_USI
        #print 'exponent is', np.exp(b_cgs)
        #print np.shape(np.exp(b_cgs))
        #print np.isinf(np.exp(b_cgs))
        if verbose == True:
            print('b_cgs is', b_cgs)
            print('be aware that {0} elements in the exponent of the Planck formula lead to an infinite exponent'.format(np.shape(np.exp(b_cgs)[np.isinf(np.exp(b_cgs))==True])[0]))
            print('denom shape is',np.shape(h_cgs*c_cgs/(wavelength_in_cm*kB_cgs*Temp)))
        #np.exp(wavelength_in_cm)
        E_cgs=coeff*2*math.pi*h_cgs*c_cgs**2/(wavelength_in_cm**5 *(np.exp(h_cgs*c_cgs/(np.float64(wavelength_in_cm)*kB_cgs*Temp)) - 1.0))
        E_Hz=coeff*2*math.pi*h_cgs*nu**3/(c_cgs**2*(np.exp(h_cgs*nu/(kB_cgs*Temp))-1.0)) #this is the planck formula in Hz ()
        E_A=E_cgs*1e-8 # because cm-1 =(1e8 A)-1
        E_mjy=1e-26*E_Hz # because 1Jy=1e-26 J/(sec*m^2*Hz) and 1J=1e7erg
        E_phot=coeff*2*math.pi*nu**2/(c_cgs**2*(np.exp(h_cgs*nu/(kB_cgs*Temp))-1.0))
    elif type.lower() == 'rj':
        if verbose == True:
            print('formula used for black body: Rayleigh-Jeans')
        E_cgs=coeff*2*math.pi*c_cgs*kB_cgs*Temp/wavelength_in_cm**4
        E_Hz=coeff*2*math.pi*kB_cgs*Temp*(nu/c_cgs)**2
        E_A = E_cgs * 1e-8  # because cm-1 =(1e8 A)-1
        E_mjy = 1e-26 * E_Hz  # because 1Jy=1e-26 J/(sec*m^2*Hz) and 1J=1e7erg
        E_phot=None # I am not sure
    else:
        print('unknown formula')
        pdb.set_trace()
    #if redshift!=None:

        wavelength_fixed=wavelength*(redshift+1)
        #print(wavelength_fixed)
        #pdb.set_trace()
    #else:
    wavelength_fixed=wavelength*(redshift+1)
    E_A_fixed=E_A/(redshift+1)
    if Ebv==None:
        spectrum_cgs=np.array(list(zip(wavelength_fixed,E_cgs)))#not sure how z influences
        spectrum_Hz=np.array(list(zip(wavelength_fixed,E_Hz)))#not sure how z influences
        spectrum_A=np.array(list(zip(wavelength_fixed,E_A_fixed)))
        spectrum_mjy=np.array(list(zip(wavelength_fixed,E_mjy)))#not sure how z influences
        spectrum_phot=np.array(list(zip(wavelength_fixed,E_phot)))#not sure how z influences
    else:
        #print(wavelength)
        #print('***')
        #print(wavelength * 1e6)
        #pdb.set_trace()
        spectrum_cgs=np.array(list(zip(wavelength_fixed,extinction.apply_extinction_to_theoretical_flux(np.array(list(zip(wavelength*1e6,E_cgs))),Ebv,R=R_ext)[:,1])))# apply_extinction_to_theoretical_flux needs wavelengths in micropmeters
        spectrum_Hz=np.array(list(zip(wavelength_fixed,extinction.apply_extinction_to_theoretical_flux(np.array(list(zip(wavelength*1e6,E_Hz))),Ebv,R=R_ext)[:,1])))
        spectrum_A = np.array(list(zip(wavelength_fixed,extinction.apply_extinction_to_theoretical_flux(np.array(list(zip(wavelength * 1e6, E_A_fixed))),Ebv,R=R_ext)[:,1])))
        #spextrum_A_befor_E=np.array(list(zip(wavelength_fixed,E_A_fixed)))
        spectrum_mjy = np.array(list(zip(wavelength_fixed,extinction.apply_extinction_to_theoretical_flux(np.array(list(zip(wavelength * 1e6, E_mjy))),Ebv,R=R_ext)[:, 1])))
        spectrum_phot = np.array(list(zip(wavelength_fixed,extinction.apply_extinction_to_theoretical_flux(np.array(list(zip(wavelength * 1e6, E_phot))),Ebv,R=R_ext)[:, 1])))
    if plot==True:
        pylab.figure()
        pylab.plot(wavelength,E_A,label='sepctrum before applying z and E')
        pylab.plot(wavelength_fixed,E_A_fixed,label='sepctrum redshifted z={0}'.format(redshift))
        pylab.plot(spectrum_A[:,0],spectrum_A[:,1], label='sepctrum redshifted z={0} and extincted'.format(redshift))
        pylab.legend()
        pylab.show()
    #print('managed till here')
    #pdb.set_trace()
    return spectrum_cgs, spectrum_Hz, spectrum_A, spectrum_mjy, spectrum_phot


