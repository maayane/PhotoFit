
import numpy as np
from . import get_filter
from . import black_body_flux_density
from . import distances_conversions
import pylab
import pyphot
import os
import logging
import shutil
import pdb
from scipy.interpolate import interp2d
import matplotlib.pyplot as plt

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def calc_black_body_flux_filters(Temp,wavelengths,filters_directory=None,Filter_vector=None,P_vector=None,Radius=None,distance_pc=None,output_txt=True,output_plot=True,lib=None,show_plots=False,Ebv=None,R_ext=None,z=0,output_file=None):
    """This code calcultes the synthetic flux \bar{f}(P) of a black_body in a given filter P (filter family and filter name), given its temperature T, Radius and distance.
    Input  :- Temperature [K]
            - wavelengths [m]
            - Filter vector: N-long array of of arrays [[filter family, filtername],[filter family, filtername],[filter family, filtername],...etc]
            the family and names are the ones provided by the pyphot library #FIX THERE IS A PROBBLEM IF FILTER VECTOR IS ONLY 1D
            - alternatively, you can give a N-long array [[filter family, filtername,P1],[filter family, filtername,P2],..,[filter family, filtername.PN]], where P are pyphot.Filter objects. This can help to speed up the code, e.g. if it's used in a mcmc
            - Radius in cm default is None
            - distance: distance of the blackbody in pc, default is 10 pc
            - output_txt: if True (default), create an output file with txt files. In codes using this function
            many tims (e.g. fit_black_body_fluxfiters.py), you may want to set it to false.
            - output_plot: if True (default), create an output file with txt files. In codes using this function
            in loops (e.g. fit_black_body_fluxfiters.py), you may want to set it to false.
            - lib: library of filters in the same format as the pyphot library. If None, uses lib = pyphot.get_library().
            in codes using this function in loops (e.g. fit_black_body_fluxfiters.py), you may want to set lib = pyphot.get_library() outside of the loop.
            - CAREFULL! Ebv and z are exinctions and redshift to be APPLIED (not corrected) to a theoretical black boday before synthetic photometry is ran on it.
            i.e. the theoretical bb is extincted and redshifted by E abd z
    Output  :- N*4 array where each line is: [filter family, filtername,effective wavelength of P,\bar{f}(P)], f is in [erg/sec/cm^2/\AA ]
    Plots and output files: -plot of lambda_mean(P), \bar{f}(P)
	Tested : ?
	    By : Maayane T. Soumagnac Nov 2016
	   URL :
	Example: sun_flux=calc_black_body_flux_filters.calc_black_body_flux_filters(T,wavelengths,Filter_vector,distance=4.8481e-6,Radius=1.)
	Link to other functions: it is the same to do calc_anyspectrum_flux_filters(black_body_flux_density) and calc_black_body_flux_filters
	Reliable:
	Reference: my notes in observationnal astronomy, pargraph on synthetic photometry,
	            //anaconda/lib/python2.7/site-packages/pyphot/phot.py Filter class
	 TO DO: give an option for the speed: decrease the length of TempVec, give the spec units as options"""
    #print('I am running calc_blacl_body_flux_filter'
    if Filter_vector is None and P_vector is None:
        print('ERROR: you need to define either Filter_vector or P_vector')
        pdb.set_trace()
    #import ipdb; ipdb.set_trace()
    if output_txt==True:
        if os.path.exists('./outputs_from_calc_black_body_flux_filters_function'):
            logger.info('output_path/txt_files did exist, I am removing it and creating a new one')
            shutil.rmtree('./outputs_from_calc_black_body_flux_filters_function')
        else:
            logger.info('the output file file did not exist yet. I am creating it now')
        os.makedirs('./outputs_from_calc_black_body_flux_filters_function')
    #wavelengths = wavelengths*(z + 1)
    #print('wavelengths are',wavelengths)
    if Radius is not None:
        #print('Radius is not none'
        #print('Raius is {0}'.format(Radius)
        #pdb.set_trace()
        Radius=distances_conversions.cm_to_solar_radius(Radius)
    #if Ebv is None:
    #    black_body_spectrum=black_body_flux_density.black_body_flux_density(Temp,wavelengths,'P',distance_pc=distance_pc,Radius=Radius)[2]/(z+1)#,Ebv=Ebv)[2]# in erg/sec/cm^2/Ang
    #print(' black_body_spectrum is {0}'.format(black_body_spectrum)
    #else:
    #    black_body_spectrum =black_body_flux_density.black_body_flux_density(Temp, wavelengths, 'P', distance_pc=distance_pc, Radius=Radius,Ebv=Ebv)[2]/(z+1)  # in erg/sec/cm^2/Ang
    #black_body_spectrum = \
    #black_body_flux_density.black_body_flux_density(Temp, wavelengths, 'P', distance_pc=distance_pc, Radius=Radius,
    #                                                Ebv=Ebv,R_ext=R_ext,redshift=z)[2]   # in erg/sec/cm^2/Ang
    black_body_spectrum = \
    black_body_flux_density.black_body_flux_density_fast(Temp, wavelengths, 'P', distance_pc=distance_pc, Radius=Radius,
                                                    Ebv=Ebv,R_ext=R_ext,redshift=z)   # in erg/sec/cm^2/Ang
    #print('T is',Temp)
    #print('R is',Radius)
    #print('distance is',distance_pc)
    #print('EBV is',Ebv)
    #print(np.shape(wavelengths))
    #print(np.shape(black_body_spectrum))
    #print('black_body_spectrum[:,1] is',black_body_spectrum[:,1])
    #print 'bbflux density is',black_body_flux_density.black_body_flux_density(Temp, wavelengths, 'P', distance_pc=distance_pc, Radius=Radius,Ebv=Ebv)[2]
    #wavelengths_AA = wavelengths* 1e10#pyphot takes wavelengths in AA
    wavelengths_AA=black_body_spectrum[:,0]*1e10
    if lib is None:
        lib = pyphot.get_library()
        #print("Library contains: ", len(lib), " filters")
    else:
        #print('the library was defined outside'
        lib=lib
    wavelengths_AA = wavelengths_AA.astype(float)


    if Filter_vector is not None:
        #print('filters_directory:',filters_directory)
        [P, wav]=get_filter.make_filter_object(Filter_vector,filters_directory=filters_directory)
        #print(P,wav)
        #pdb.set_trace()
        output_array = np.empty([np.shape(Filter_vector)[0], 4], dtype=object)
        string = np.empty(np.shape(Filter_vector)[0], dtype=object)
        for i,j in enumerate(Filter_vector):
            '''
            #print j
            if j[0].lower()=='ptf_p48':
               if j[1].lower()=='r':
                   print('You gave the R filter of the PTF_P48 family')
                   Transmission=np.genfromtxt('/Users/maayanesoumagnac/Maayane_Astro_python_library/data/filters/P48_R_T.rtf',delimiter=None)
                   P=pyphot.Filter(Transmission[:,0], Transmission[:,1], name='PTF_P48_R', dtype='photon', unit='Angstrom')
            elif j[0].lower() == 'galex':
                if j[1].lower() == 'galex_nuv':
                    print('You gave the nuv filter of the galex family')
                    Transmission = np.genfromtxt(
                        '/Users/maayanesoumagnac/Maayane_Astro_python_library/data/filters//GALEX_NUV.dat',
                        delimiter=None)
                    P = pyphot.Filter(Transmission[:, 0], Transmission[:, 1], name='GALEX_NUV', dtype='photon',
                                      unit='Angstrom')
            else:
                f = lib.find(j[0].lower())#filter family
                P = lib[j[1]]#filter name
            # this calculates \bar{f}(P) through through filter P. wavelength must be given in the same units as wavelength units defined in filters
#            print('wavelengths in AA are',wavelengths_AA
#            print('the bb spectrum is', black_body_spectrum[:, 1]
            #pdb.set_trace()
            '''
            #P_vector = dict()
            #P_vector['filter_family'] = Filter_vector[:, 0]
            #P_vector['filter_name'] = Filter_vector[:, 1]
            #P_vector['filter_object'] = []
            fluxes = P['filter_object'][i].get_flux(black_body_spectrum[:, 0]*1e10, black_body_spectrum[:, 1])
            #print(fluxes)
            #pdb.set_trace()
            output_array[i,0]=Filter_vector[i,0]
            output_array[i,1]=Filter_vector[i,1]
            output_array[i, 2] = P['filter_object'][i].cl.item()#le eff_wl de Eran (auquel je veux pouvoir comparer) est le cl de pyphot
            #output_array[i,2]=P['filter_object'][i].leff.item()
            output_array[i,3]=fluxes
            #print('fluxes are',fluxes)
            #pdb.set_trace()
            string[i]=output_array[i,1]+'\n'+r'$(\lambda_c=$'+str(round(output_array[i,2],3))+'$)$'
            #if output_txt == True:

    else:
        #print('P vector is',P_vector)
        output_array = np.empty([np.shape(P_vector)[0], 4], dtype=object)
        string = np.empty(np.shape(P_vector)[0], dtype=object)
        for i, fil in enumerate(P_vector):
            #print P_vector[i,2]
            #print(wavelengths_AA)
            #print black_body_spectrum[:,1]
            #pdb.set_trace()
            #print(P_vector[i,2])
            fluxes = P_vector[i,2].get_flux(black_body_spectrum[:, 0]*1e10, black_body_spectrum[:, 1])
            #print('fluxes are',fluxes)
            #print fluxes
            output_array[i, 0] = P_vector[i, 0]
            output_array[i, 1] = P_vector[i, 1]
            output_array[i, 2] = P_vector[i, 2].cl.item()# le eff_wl de Eran (auquel je veux pouvoir comparer) est le cl de pyphot
            #output_array[i, 2] = P_vector[i, 2].leff.item()
            #pdb.set_trace()
            output_array[i, 3] = fluxes
            string[i] = output_array[i, 1] + '\n' + r'$(\lambda_c=$' + str(round(output_array[i, 2], 3)) + '$)$'
            #if output_txt == True:
            #    np.savetxt('./outputs_from_calc_black_body_flux_filters_function/output_array.txt', output_array,
            #               header=r'Filter family, Filter name, effective wavelength $(\AA)$, flux $\bar{f}(P)$ $[erg/sec/cm^2/\AA ]$',
            #               fmt="%s")
    if output_txt==True:
        np.savetxt('./outputs_from_calc_black_body_flux_filters_function/output_array.txt', output_array,header=r'Filter family, Filter name, effective wavelength $(\AA)$, flux $\bar{f}(P)$ $[erg/sec/cm^2/\AA ]$',fmt="%s")
    if output_plot==True:
        pylab.figure()
        pylab.plot(output_array[:,2], output_array[:,3], 'ro', label='synthetic flux through input filters')
        pylab.plot(black_body_spectrum[:,0]*1e10,black_body_spectrum[:,1],'b-', label='theoretical black body with input T ans R')
        #if Ebv != None:
        #    pylab.plot(output_array[:, 2], output_array[:,4], 'mo', label='without extinction')
            #pylab.plot(1 / wavelengths_allen, AlA_Alen, 'b', label='Allen')
        pylab.xlabel(r'Central wavelength of filters ($\AA$)')
        pylab.ylabel(r'flux $\bar{f}(P)$, $[erg/sec/cm^2/\AA ]$')
        pylab.title('result of calc_blac_body_flux_filters.py')
        labels = string[:]
        pylab.xticks(output_array[:,2],labels)#,rotation='vertical')
        pylab.legend(loc=4)
        pylab.savefig(output_file+'/spectrum_plot_from_calc_black_body_flux_filters.pdf', facecolor='w', edgecolor='w',
                  orientation='portrait', papertype=None, format='pdf', transparent=False, bbox_inches=None,
                  pad_inches=0.1)
    if show_plots==True:
        pylab.show()
    #print('output array is',output_array)
    return output_array #N*4 array where each line is: [filter family, filtername,effective wavelength of P,\bar{f}(P)], f is in [erg/sec/cm^2/\AA ]

def Integrate_flux_in_filter(black_body_spectrum,filters_directory,Filter_vector=None,P_vector=None):
    if Filter_vector is not None:
        [P, _ ]=get_filter.make_filter_object(Filter_vector,filters_directory=filters_directory)
        output_array = np.empty([np.shape(Filter_vector)[0], 4], dtype=object)
        string = np.empty(np.shape(Filter_vector)[0], dtype=object)
        for i in range(len(Filter_vector)):
            fluxes = P['filter_object'][i].get_flux(black_body_spectrum[:, 0]*1e10, black_body_spectrum[:, 1])

            output_array[i,0]=Filter_vector[i,0]
            output_array[i,1]=Filter_vector[i,1]
            output_array[i, 2] = P['filter_object'][i].cl.item()#le eff_wl de Eran (auquel je veux pouvoir comparer) est le cl de pyphot
            output_array[i,3]=fluxes
            string[i]=output_array[i,1]+'\n'+r'$(\lambda_c=$'+str(round(output_array[i,2],3))+'$)$'
    else:
        output_array = np.empty([np.shape(P_vector)[0], 4], dtype=object)
        string = np.empty(np.shape(P_vector)[0], dtype=object)
        for i in range(len(P_vector)):
            fluxes = P_vector[i,2].get_flux(black_body_spectrum[:, 0]*1e10, black_body_spectrum[:, 1])
            output_array[i, 0] = P_vector[i, 0]
            output_array[i, 1] = P_vector[i, 1]
            output_array[i, 2] = P_vector[i, 2].cl.item()# le eff_wl de Eran (auquel je veux pouvoir comparer) est le cl de pyphot
            output_array[i, 3] = fluxes
            string[i] = output_array[i, 1] + '\n' + r'$(\lambda_c=$' + str(round(output_array[i, 2], 3)) + '$)$'
    return output_array,string


def Plot_bb_fits(black_body_spectrum,string,output_array,output_file,show_plots):
    pylab.figure()
    pylab.plot(output_array[:,2], output_array[:,3], 'ro', label='synthetic flux through input filters')
    pylab.plot(black_body_spectrum[:,0]*1e10,black_body_spectrum[:,1],'b-', label='theoretical black body with input T ans R')
    #if Ebv != None:
    #    pylab.plot(output_array[:, 2], output_array[:,4], 'mo', label='without extinction')
        #pylab.plot(1 / wavelengths_allen, AlA_Alen, 'b', label='Allen')
    pylab.xlabel(r'Central wavelength of filters ($\AA$)')
    pylab.ylabel(r'flux $\bar{f}(P)$, $[erg/sec/cm^2/\AA ]$')
    pylab.title('result of calc_blac_body_flux_filters.py')
    labels = string[:]
    pylab.xticks(output_array[:,2],labels)#,rotation='vertical')
    pylab.legend(loc=4)
    pylab.savefig(output_file+'/spectrum_plot_from_calc_black_body_flux_filters.pdf', facecolor='w', edgecolor='w',
              orientation='portrait', papertype=None, format='pdf', transparent=False, bbox_inches=None,
              pad_inches=0.1)
              
    if show_plots==True:
        pylab.show()

def calc_black_body_flux_filters_clean(Temp,Radius,Ebv,distance_pc,Filter_vector=None,P_vector=None,wavelengths=np.arange(1e-7, 3e-6, 5e-9),output_txt=True,filters_directory=None,output_plot=True,lib=None,show_plots=False,R_ext=3.08,z=0,output_file=None,**kwargs):
    """This code calcultes a grid of synthetic fluxes \bar{f}(P) of a blackbody in a given filter P (filter family and filter name), given its temperature T, Radius and distance.
    Input  :- Temperature [K]
            - wavelengths [m]
            - Filter vector: N-long array of of arrays [[filter family, filtername],[filter family, filtername],[filter family, filtername],...etc]
            the family and names are the ones provided by the pyphot library #FIX THERE IS A PROBBLEM IF FILTER VECTOR IS ONLY 1D
            - alternatively, you can give a N-long array [[filter family, filtername,P1],[filter family, filtername,P2],..,[filter family, filtername.PN]], where P are pyphot.Filter objects. This can help to speed up the code, e.g. if it's used in a mcmc
            - Radius in cm default is None
            - distance: distance of the blackbody in pc, default is 10 pc
            - output_txt: if True (default), create an output file with txt files. In codes using this function
            many tims (e.g. fit_black_body_fluxfiters.py), you may want to set it to false.
            - output_plot: if True (default), create an output file with txt files. In codes using this function
            in loops (e.g. fit_black_body_fluxfiters.py), you may want to set it to false.
            - lib: library of filters in the same format as the pyphot library. If None, uses lib = pyphot.get_library().
            in codes using this function in loops (e.g. fit_black_body_fluxfiters.py), you may want to set lib = pyphot.get_library() outside of the loop.
            - CAREFULL! Ebv and z are exinctions and redshift to be APPLIED (not corrected) to a theoretical black boday before synthetic photometry is ran on it.
            i.e. the theoretical bb is extincted and redshifted by E abd z
    Output  :- N*4 array where each line is: [filter family, filtername,effective wavelength of P,\bar{f}(P)], f is in [erg/sec/cm^2/\AA ]
    Plots and output files: -plot of lambda_mean(P), \bar{f}(P)
	Tested : ?
	    By : Maayane T. Soumagnac Nov 2016
	   URL :
	Example: sun_flux=calc_black_body_flux_filters.calc_black_body_flux_filters(T,wavelengths,Filter_vector,distance=4.8481e-6,Radius=1.)
	Link to other functions: it is the same to do calc_anyspectrum_flux_filters(black_body_flux_density) and calc_black_body_flux_filters
	Reliable:
	Reference: my notes in observationnal astronomy, pargraph on synthetic photometry,
	            //anaconda/lib/python2.7/site-packages/pyphot/phot.py Filter class
	 TO DO: give an option for the speed: decrease the length of TempVec, give the spec units as options"""
    #print('I am running calc_blacl_body_flux_filter'
    if Filter_vector is None and P_vector is None:
        print('ERROR: you need to define either Filter_vector or P_vector')
        pdb.set_trace()
    if output_txt==True:
        if os.path.exists('./outputs_from_calc_black_body_flux_filters_function'):
            logger.info('output_path/txt_files did exist, I am removing it and creating a new one')
            shutil.rmtree('./outputs_from_calc_black_body_flux_filters_function')
        else:
            logger.info('the output file file did not exist yet. I am creating it now')
        os.makedirs('./outputs_from_calc_black_body_flux_filters_function')

    if Radius is not None:
        Radius=distances_conversions.cm_to_solar_radius(Radius)
    black_body_spectrum = \
    black_body_flux_density.black_body_flux_density_fast(Temp, wavelengths, 'P', distance_pc=distance_pc, Radius=Radius,
                                                    Ebv=Ebv,R_ext=R_ext,redshift=z)   # in erg/sec/cm^2/Ang
    if lib is None:
        lib = pyphot.get_library()
    else:
        lib=lib

    output_array,string=Integrate_flux_in_filter(black_body_spectrum,filters_directory=filters_directory,Filter_vector=Filter_vector,P_vector=P_vector)
            
    if output_txt==True:
        outpath='./outputs_from_calc_black_body_flux_filters_function/output_array.txt'
        head=r'Filter family, Filter name, effective wavelength $(\AA)$, flux $\bar{f}(P)$ $[erg/sec/cm^2/\AA ]$'
        np.savetxt(outpath, output_array,header=head,fmt="%s")
    if output_plot==True:
        Plot_bb_fits(black_body_spectrum,string,output_array,output_file,show_plots)
    return output_array #N*4 array where each line is: [filter family, filtername,effective wavelength of P,\bar{f}(P)], f is in [erg/sec/cm^2/\AA ]




def generate_grid(T_lims,R_lims,Ebv,z,grid_path,filters_directory,distance_pc,lib,R_ext,Filter_vector=None,P_vector=None, res=100):
    logT_min=np.log10(np.min(T_lims))
    logT_max=np.log10(np.max(T_lims))
    logR_min=np.log10(np.min(R_lims))
    logR_max=np.log10(np.max(R_lims))

    T_vec=np.linspace(logT_min,logT_max,res)
    R_vec=np.linspace(logR_min,logR_max,res)
    grid={}
    f_grid={}
    for j in range(len(P_vector)):
        filt=P_vector[j,0]+'+'+P_vector[j,1]
        grid[filt]=np.zeros((res,res))

    
    wl=np.arange(1e-7, 3e-6, 5e-9)
    Tv, Rv = np.meshgrid(10**T_vec, 10**R_vec)
    for n in range(res):
        for m in range(res):
            TT=Tv[n,m]
            RR=Rv[n,m]
            out = calc_black_body_flux_filters(TT,wl,Filter_vector=Filter_vector, P_vector=P_vector,
                                                                            Radius=RR, distance_pc=distance_pc,
                                                                            output_plot=False,
                                                                            show_plots=False, output_txt=False,
                                                                            lib=lib, z=z,Ebv=Ebv,R_ext=R_ext)

            for j in range(len(P_vector)):
                filt=P_vector[j,0]+'+'+P_vector[j,1]
                grid[filt][n,m]=out[j, 3]	


    #head='T_lim={0},R_lims={1},Ebv={2}'.format(T_lims,R_lims,Ebv)
    #np.savetxt(grid_path,grid,header=head,fmt="%s")
    return Tv, Rv, grid


