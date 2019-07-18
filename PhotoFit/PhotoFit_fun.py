
"""*******************************************************
PhotoFit: calculates the evolution in time of the effective radius
and temperature, by interpolating to a common epoch, and fitting a
blackbody (after applying E and z corrections and synthetic photometry)
to each SED.
******************************************************"""
print(__doc__)


import numpy as np
import pylab
import pdb
import os
import shutil
import logging
import pyphot
from pathlib import Path
import matplotlib.pyplot as plt
import math
from . import interpolate_errors
from . import read_data_from_file
from . import fit_black_body_flux_filters_mcmc
from . import fit_black_body_flux_filters
from . import get_filter
from . import black_body_flux_density
from . import distances_conversions
from . import energy_conversions
from . import fitter_general

#print('Hi')
#pdb.set_trace()

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

__all__=['calculate_T_and_R_in_time','plot_T_and_R_in_time','plot_L_in_time','plot_SEDs']

#mcmc=params.mcmc

#output_file_mcmc=params.output_file_mcmc
#output_file_linear=params.output_file_linear
#output_file_interpolation=params.output_file_interpolation

#filters_directroy=params.filters_directory

#if mcmc==True:
#    output=output_file_mcmc
#else:
#    output=output_file_linear
#if os.path.exists(output):
#    print('the output directory, '+output + ' exists already')
#else:
#    os.mkdir(output)

########## Informations on the object ##########

#redshift = params.z
#distance_modulus=params.distance_modulus
#distance_pc=distances_conversions.DM_to_pc(distance_modulus)   # sn distance in pc
#explosion_date=params.explosion_date
#EBV=params.EBV

########################################### Read and process the photometric data ############################################

# make sure you don't have cases with repeated x,y, and different yerr
#lower_limit_on_flux=params.lower_limit_on_flux
def calculate_T_and_R_in_time(data_file=None,dates_file=None,already_run_interp_errors_from_param=None,already_run_mcmc=None,
                              already_run_matrix=None,show_underlying_plots=True,verbose=False,redshift=None,
                              distance_modulus=None,explosion_date=None,EBV=None,output=None,filters_directory=None,
                              mcmc=False,output_file_interpolation=None,lower_limit_on_flux=None,num_iterations=None,
                              num_steps=None,nwalkers=None,csm=False,already_run_fit=None,excluded_bands=[]):

    #print('num_steps is',num_steps)
    #print('nwalkers is',nwalkers)
    #pdb.set_trace()
    if data_file is None:
        print('you need to procide a data_file')
        exit()
    if dates_file is None:
        print('you need to procide a dates file')
        exit()
    if already_run_interp_errors_from_param is None:
        print('you need to procide a already_run_interp_errors_from_param array')
        exit()
    if already_run_mcmc is None:
        print('you need to procide a already_run_mcmc param')
        exit()
    if already_run_matrix is None:
        print('you need to procide a already_run_matrix')
        exit()

    distance_pc = distances_conversions.DM_to_pc(distance_modulus)

    data_dico = \
    read_data_from_file.read_data_into_numpy_array(data_file, header=True, delimiter=',', no_repeat_rows=True)[2]


    data_dicts=dict()
    for i,j in enumerate(data_dico['filter']):
        data_dicts[j]=dict()
        data_dicts[j]['jd']=data_dico['jd'][(data_dico['filter']==j) & (data_dico['flux']>lower_limit_on_flux)]
        #data_dicts[j]['mag'] = data_dico['mag'][(data_dico['filter']==j) & (data_dico['flux']>lower_limit_on_flux)]
        data_dicts[j]['flux'] = data_dico['flux'][(data_dico['filter']==j) & (data_dico['flux']>lower_limit_on_flux)]
        #data_dicts[j]['magerr'] = data_dico['magerr'][(data_dico['filter']==j) & (data_dico['flux']>lower_limit_on_flux)]
        data_dicts[j]['fluxerr'] = data_dico['fluxerr'][(data_dico['filter']==j) & (data_dico['flux']>lower_limit_on_flux)]
        #data_dicts[j]['absmag'] = data_dico['absmag'][(data_dico['filter']==j) & (data_dico['flux']>lower_limit_on_flux)]
        #data_dicts[j]['absmagerr'] = data_dico['absmagerr'][(data_dico['filter']==j) & (data_dico['flux']>lower_limit_on_flux)]
        #data_dicts[j]['filter'] = data_dico['filter'][(data_dico['filter']==j) & (data_dico['flux']>lower_limit_on_flux)]
        #data_dicts[j]['instr'] = data_dico['instr'][(data_dico['filter']==j) & (data_dico['flux']>lower_limit_on_flux)]

    print('the bands of your multiple-bands data are:',data_dicts.keys())
    print('the bands to be excluded are:',excluded_bands)

    ########################################### Make the spectrum to fit for each epoch ############################################

    ##### definition of interpolation dates

    interpolation_dates=np.genfromtxt(dates_file)#param.dates_file

    print('**********************************')
    print('******* INTERPOLATION STEP *******')
    print('**********************************')


    print('I am interpolating on dates', interpolation_dates)
    print('If I subtract the explosion dates you provided, this comes to interpolating on t_exp + ', interpolation_dates - explosion_date)
    condition=dict()
    JD_basis_interpolation=dict()
    already_run_interp_errors=dict()
    jd_flux_fluxerr=dict()
    interp_and_errors_array=dict()
    interp=dict()

    for k in [x for x in data_dicts.keys() if x not in excluded_bands]:
        already_run_interp_errors[k] = already_run_interp_errors_from_param[k]#params.already_run_interp_errors[k]

    ### Interpolate all data on the interpolation dates dates ####

    #print('data_dicts.keys() are',data_dicts.keys())
    #print('without the excluded bands:',[x for x in data_dicts.keys() if x not in excluded_bands])
    #pdb.set_trace()

    if os.path.exists(output_file_interpolation) == False:
        os.mkdir(output_file_interpolation)
    #[x for x in a if x not in [2, 3, 7]]

    for k in [x for x in data_dicts.keys() if x not in excluded_bands]:
        #if verbose==True:
        #print('I am looking at the {0} band of the data'.format(k))
        if already_run_interp_errors[k]==True:
            print('******')
            print('Band {0}: You indicated that the interpolation was already done for this band, so I am not doing it again.'.format(k))
        else:
            print('******')
            print('Band {0}: interpolating now'.format(
                k))
        #print('the days of {0} are {1}'.format(k,data_dicts[k]['jd']-explosion_date))
        condition[k]=np.logical_and(interpolation_dates>=np.min(data_dicts[k]['jd']),interpolation_dates<=np.max(data_dicts[k]['jd']))
        JD_basis_interpolation[k] = interpolation_dates[condition[k]]
        if verbose == True:
            print('The interpolation dates relevant to this band are {1}'.format(k,JD_basis_interpolation[k]-explosion_date))
        #print('and interpolating on {0}'.format(interpolation_dates-explosion_date))
        #print('len(JD_basis_interpolation[k]) is',len(JD_basis_interpolation[k]))
        if len(JD_basis_interpolation[k])>0:


            a=np.array(list(zip(data_dicts[k]['jd'],data_dicts[k]['flux'],data_dicts[k]['fluxerr'])))
            jd_flux_fluxerr[k]=a[a[:, 0].argsort()]
            output_path=output_file_interpolation+'/errors_interpolation_results_{0}'.format(k)
            if already_run_interp_errors[k]==False:
                if os.path.exists(output_path)==True:
                    shutil.rmtree(output_path)
                os.mkdir(output_path)
            interp_and_errors_array[k]=interpolate_errors.interpolate_errors\
            (jd_flux_fluxerr[k],JD_basis_interpolation[k],output_path=output_file_interpolation+'/errors_interpolation_results_{0}'.format(k),
             already_run=already_run_interp_errors[k],show_plot=False,title='{0}'.format(k),verbose=verbose,band_name='{0}'.format(k))#array with days_nuv, interp1d P48, interp P48 with another method, lower limit 1sigma, upper limit 1 sigma
            interp[k]=dict()
            #print('interp_and_errors_array[k] is',interp_and_errors_array[k])
            #print('np.shape(interp_and_errors_array[k]) is',np.shape(interp_and_errors_array[k]))
            #print('k is',k)

            for i, j in enumerate(JD_basis_interpolation[k]):
                if interp_and_errors_array[k].ndim==1:
                    interp[k][str(round(j - explosion_date, 8))] = interp_and_errors_array[k][:]
                else:
                    interp[k][str(round(j-explosion_date,8))]=interp_and_errors_array[k][i,:]
            #print('interp_and_errors_array[k] is',interp_and_errors_array[k])
            #print('JD_basis_interpolation[k] is',JD_basis_interpolation[k]-explosion_date)


    Spectra=[]
    if verbose == True:
        print('I am creating the Spectra')
    for i,j in enumerate(interpolation_dates):
        #print('********')
        #print('i is',i)
        #print('********')
        epoch=dict()
        epoch['time']=round(interpolation_dates[i]-explosion_date,8)#jd_flux_fluxerr_UVW2[2,0]#
        #print('epoch[time] is',str(epoch['time']))
        for k in [x for x in data_dicts.keys() if x not in excluded_bands]:
            if len(JD_basis_interpolation[k]) > 0:
                if str(epoch['time']) in interp[k].keys():
                    epoch[k]=interp[k][str(epoch['time'])][1] #flux
                    epoch[k+'_err']=[interp[k][str(epoch['time'])][3],interp[k][str(epoch['time'])][4]] # lower and higher limit of 2 sigma area. the error is half of this
        Spectra.append(epoch)

    #for i,j in enumerate(Spectra):
    #    print('At day {0}, the Spectrum is {1}'.format(j['time'],Spectra[i]))
    #    pdb.set_trace()

    print('***************************************************')
    print('******* IMPORTANT CHECK BEFORE YOU PROCEED ********')
    print('To check that your interpolation is correct: go to the directory where you store the interpolation results. For each band, PhotoFit produces a plot called "Plot_w_interpolated_errors.png". Check, for all bands, that '
          'the green and blue points are superimposed on this plot (or almost superimposed). If you do not see any blue, you are good :-) If this is not the case, try changing your dates file so that none of the dates exactely coincides with one of'
          'your observation dates. This should solve the problem.')
    print('***************************************************')


    ########################################### Filters definition ###########################################

    # Filter_vector for all the filters
    Filter_vector = np.empty([21, 2], dtype=object)
    Filter_vector[0] = [str('Swift'), str('UVW1')]
    Filter_vector[1] = [str('Swift'), str('UVW2')]
    Filter_vector[2] = [str('Swift'), str('UVM2')]
    Filter_vector[3] = [str('Swift'), str('u_swift')]
    Filter_vector[4] = [str('Swift'), str('b_swift')]
    Filter_vector[5] = [str('Swift'), str('v_swift')]
    Filter_vector[6] = [str('johnson'), str('u_johnson')]
    Filter_vector[7] = [str('johnson'), str('v_johnson')]
    Filter_vector[8] = [str('johnson'), str('b_johnson')]
    Filter_vector[9] = [str('ztf_p48'), str('r_p48')]
    Filter_vector[10] = [str('ztf_p48'), str('g_p48')]
    Filter_vector[11] = [str('ztf_p48'), str('i_p48')]
    Filter_vector[12] = [str('SDSS'), str('g_sdss')]
    Filter_vector[13] = [str('SDSS'), str('i_sdss')]
    Filter_vector[14] = [str('SDSS'), str('z_sdss')]
    Filter_vector[15] = [str('SDSS'), str('r_sdss')]
    Filter_vector[16] = [str('SDSS'), str('u_sdss')]
    Filter_vector[17] = [str('cousin'), str('i_cousin')]
    Filter_vector[18] = [str('cousin'), str('r_cousin')]
    Filter_vector[19] = [str('2MASS'), str('j_2mass')]
    Filter_vector[20] = [str('2MASS'), str('h_2mass')]

    lib = pyphot.get_library()
    #print("Library contains: ", len(lib), " filters")

    [P,wavelengths_filter]=get_filter.make_filter_object(Filter_vector,filters_directory=filters_directory)
    if verbose==True:
        print('wavelengths are',wavelengths_filter)

    ################################## plot the x points spectra with error bars ###########################################
    print('**********************************')
    print('******* PLOTTING THE SEDs ********')
    print('**********************************')
    if show_underlying_plots==True:
        print('You have set the parameter "show_underlying_plot" to True, so I will show you the SEDs before we proceed')
    else:
        print('You have set the parameter "show_underlying_plot" to False, so I will not show you the SEDs and store them instead in {0}/day_%/SED_date_%.png'.format(output))
    # TO DO: put in parameters file
    color=dict()
    color['UVW1']='blue'
    color['UVW2']='black'
    color['UVM2']='purple'
    color['u_swift']='orange'
    color['b_swift']='cyan'
    color['v_swift']='magenta'
    color['r_sdss']='red'
    color['i_sdss']='black'
    color['g_sdss']='green'
    color['z_sdss'] = 'grey'
    color['u_sdss'] = 'orange'
    color['r_p48']='red'
    color['g_p48']='green'
    color['i_p48']='black'
    color['u_johnson']='orange'
    color['b_johnson']='cyan'
    color['v_johnson']='magenta'
    color['r_cousin']='red'
    color['i_cousin']='black'
    color['h_2mass']='darkcyan'
    color['j_2mass']='sienna'

    symbol=dict()
    symbol['UVW1']='d'
    symbol['UVW2']='d'
    symbol['UVM2']='d'
    symbol['u_swift']='d'
    symbol['b_swift']='d'
    symbol['v_swift']='d'
    symbol['r_cousin']='v'
    symbol['i_cousin']='v'
    symbol['r_sdss']='*'
    symbol['i_sdss']='*'
    symbol['g_sdss']='*'
    symbol['z_sdss'] = '*'
    symbol['u_sdss'] = '*'
    symbol['r_p48']='o'
    symbol['g_p48']='o'
    symbol['i_p48'] = 'o'
    symbol['j_2mass']='x'
    symbol['h_2mass']='x'
    symbol['u_johnson']='s'
    symbol['b_johnson']='s'
    symbol['v_johnson']='s'

    for i,j in enumerate(Spectra):#j is an epoch, epoch['UVW1']=flux, epoch['UVW1_err']=fluxerr
        if verbose==True:
            print('the keys of this epoch are',j.keys())
        pylab.figure()
        for k in [x for x in data_dicts.keys() if x not in excluded_bands]:
            if k in j.keys():
                if verbose==True:
                    print('the filter name {0} is therefore in epoch'.format(k))
                    print('epoch[{0}] is {1}'.format(k,j[k]))
                pylab.plot(wavelengths_filter[k], j[k], symbol[k], color=color[k],label=k)
                pylab.vlines(wavelengths_filter[k], j[k + '_err'][0], j[k + '_err'][1], color=color[k])
        pylab.legend()
        pylab.title('flux spectrum on JD={0}'.format(j['time']))
        pylab.grid()
        pylab.ylabel('flux $F\; [erg/s/cm^2/\AA]$', fontsize=20)
        pylab.xlabel(r'wavelength [$\AA$]', fontsize=20)
        pylab.tight_layout()
        if os.path.exists(output + '/day_' + str(round(j['time'], 3))):
            if verbose==True:
                print(output + '/day_' + str(round(j['time'], 3)) + 'exists')
        else:
            os.mkdir(output + '/day_' + str(round(j['time'], 3)))
        pylab.savefig(output + '/day_'+str(round(j['time'],3))+'/SED_date_'+str(round(j['time'],3))+'.png',
                      facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format='png',
                      transparent=False,
                      bbox_inches=None, pad_inches=0.5)


    if show_underlying_plots==True:
        pylab.show()
            #pdb.set_trace()

    ################################################# Fit black body using mcmc/linear fit ###############################

    print('**********************************')
    print('************ FIT STEP ************')
    print('**********************************')

    already_run_mcmc=already_run_mcmc#params.already_run_mcmc
    if already_run_fit is None:
        already_run_fit=[False]*len(Spectra)
    else:
        already_run_fit=already_run_fit
    already_run_matrix=already_run_matrix#params.already_run_matrix

    if mcmc==True:
        fast=[False]*len(Spectra)
        already_run=already_run_mcmc
    else:
        already_run=already_run_matrix

    if already_run==False:
        Best=np.zeros((np.shape(Spectra)[0],8))
        #print(Best)
        #pdb.set_trace()
        for i,j in enumerate(Spectra): #goes through all epochs
            Spectrum = []
            if 'UVW2' in j.keys():
                Spectrum.append(np.array(['Swift', 'UVW2',j['UVW2'],max(j['UVW2_err'][1]-j['UVW2'],j['UVW2']-j['UVW2_err'][0])],dtype=object))#the error is half od the length upper_sigma-lower_sigma
            if 'UVW1' in j.keys():
                Spectrum.append(np.array(['Swift', 'UVW1',j['UVW1'],max(j['UVW1_err'][1]-j['UVW1'],j['UVW1']-j['UVW1_err'][0]) ],dtype=object))
            if 'UVM2' in j.keys():
                Spectrum.append(np.array(['Swift', 'UVM2',j['UVM2'],max(j['UVM2_err'][1]-j['UVM2'],j['UVM2']-j['UVM2_err'][0])],dtype=object))
            if 'u_swift' in j.keys():
                Spectrum.append(np.array(['Swift', 'u_swift',j['u_swift'],max(j['u_swift_err'][1]-j['u_swift'],j['u_swift']-j['u_swift_err'][0]) ],dtype=object))
            if 'v_swift' in j.keys():
                Spectrum.append(np.array(['Swift', 'v_swift',j['v_swift'],max(j['v_swift_err'][1]-j['v_swift'],j['v_swift']-j['v_swift_err'][0]) ],dtype=object))
            if 'b_swift' in j.keys():
                Spectrum.append(np.array(['Swift', 'b_swift',j['b_swift'],max(j['b_swift_err'][1]-j['b_swift'],j['b_swift']-j['b_swift_err'][0]) ],dtype=object))
            if 'r_p48' in j.keys():
                Spectrum.append(np.array(['ztf_p48', 'r_p48',j['r_p48'],max(j['r_p48_err'][1]-j['r_p48'],j['r_p48']-j['r_p48_err'][0])],dtype=object))
            if 'g_p48' in j.keys():
                Spectrum.append(np.array(['ztf_p48', 'g_p48',j['g_p48'],max(j['g_p48_err'][1]-j['g_p48'],j['g_p48']-j['g_p48_err'][0])],dtype=object))
            if 'i_p48' in j.keys():
                Spectrum.append(np.array(['ztf_p48', 'i_p48',j['i_p48'],max(j['i_p48_err'][1]-j['i_p48'],j['i_p48']-j['i_p48_err'][0])],dtype=object))
            if 'r_sdss' in j.keys():
                Spectrum.append(np.array(['sdss', 'r_sdss',j['r_sdss'],max(j['r_sdss_err'][1]-j['r_sdss'],j['r_sdss']-j['r_sdss_err'][0])],dtype=object))
            if 'g_sdss' in j.keys():
                Spectrum.append(np.array(['sdss', 'g_sdss',j['g_sdss'],max(j['g_sdss_err'][1]-j['g_sdss'],j['g_sdss']-j['g_sdss_err'][0])],dtype=object))
            if 'i_sdss' in j.keys():
                Spectrum.append(np.array(['sdss', 'i_sdss',j['i_sdss'],max(j['i_sdss_err'][1]-j['i_sdss'],j['i_sdss']-j['i_sdss_err'][0])],dtype=object))
            if 'z_sdss' in j.keys():
                Spectrum.append(np.array(['sdss', 'z_sdss',j['z_sdss'],max(j['z_sdss_err'][1]-j['z_sdss'],j['z_sdss']-j['z_sdss_err'][0])],dtype=object))
            if 'u_sdss' in j.keys():
                Spectrum.append(np.array(['sdss', 'u_sdss', j['u_sdss'],
                                          max(j['u_sdss_err'][1] - j['u_sdss'], j['u_sdss'] - j['u_sdss_err'][0])],
                                         dtype=object))

            if 'u_johnson' in j.keys():
                Spectrum.append(np.array(['johnson', 'u_johnson',j['u_johnson'],max(j['u_johnson_err'][1]-j['u_johnson'],j['u_johnson']-j['u_johnson_err'][0])],dtype=object))
            if 'b_johnson' in j.keys():
                Spectrum.append(np.array(['johnson', 'b_johnson', j['b_johnson'], max(j['b_johnson_err'][1]-j['b_johnson'],j['b_johnson']-j['b_johnson_err'][0])],dtype=object))
            if 'v_johnson' in j.keys():
                Spectrum.append(np.array(['johnson', 'v_johnson', j['v_johnson'], max(j['v_johnson_err'][1]-j['v_johnson'],j['v_johnson']-j['v_johnson_err'][0])],dtype=object))
            if 'i_cousin' in j.keys():
                Spectrum.append(np.array(
                    ['cousin', 'i_cousin', j['i_cousin'],  max(j['i_cousin_err'][1]-j['i_cousin'],j['i_cousin']-j['i_cousin_err'][0])],dtype=object))
            if 'r_cousin' in j.keys():
                Spectrum.append(np.array(
                    ['cousin', 'r_cousin', j['r_cousin'],  max(j['r_cousin_err'][1]-j['r_cousin'],j['r_cousin']-j['r_cousin_err'][0])],dtype=object))
            if 'h_2mass' in j.keys():
                Spectrum.append(np.array(['2MASS', 'h_2mass', j['h_2mass'],  max(j['h_2mass_err'][1]-j['h_2mass'],j['h_2mass']-j['h_2mass_err'][0])],dtype=object))
            if 'j_2mass' in j.keys():
                Spectrum.append(np.array(['2MASS', 'j_2mass', j['j_2mass'],  max(j['j_2mass_err'][1]-j['j_2mass'],j['j_2mass']-j['j_2mass_err'][0])],dtype=object))

            Spectrum_right_format=np.array(Spectrum)
            if csm==False:
                if (j['time']< 1):
                    hitemp = 3e5
                    hirad=1e15
                elif(j['time'] < 2):
                    hitemp = 5e4
                    hirad = 2e15
                else:
                    hitemp = 2e4
                    hirad = 2e15
            else:
                if (j['time']< 1):
                    hitemp = 3e5
                    hirad=2e15#7e15
                elif(j['time'] < 3):
                    hitemp = 5e4
                    hirad =2e15#7e15
                else:
                    hitemp = 4e4
                    hirad = 2e15#8e15
            #print('hitemp is',hitemp)
            #pdb.set_trace()

            num_iterations=num_iterations
            #math.log10(3000), math.log10(hitemp)

            if mcmc==True:

                print('******* method chosen: MCMC *******')

                if already_run_fit[i]==False:
                    if fast[i]==False:
                        [best_temp, best_radius, best_luminosity,best_coeff,chi_square,chi_square_dof]=fit_black_body_flux_filters_mcmc.fit_black_body_flux_filters_mcmc\
                        (Spectrum_right_format,triangle_plot_title=r'$JD-t_{ref}=$'+str(round(j['time'],2)),nwalkers=nwalkers,num_steps=num_steps,num_winners=20,
                         already_run_mcmc=already_run_mcmc,already_run_calc_all_chis=False,Temp_prior=np.array([3000,hitemp]),
                         Radius_prior=np.array([1e14,hirad]),initial_conditions=np.array([15000,3e14]),distance_pc=distance_pc,Ebv=EBV,ndof=None,show_plot=False,output_mcmc=output+'/day_'+str(round(j['time'],3)),show_mag_AB=True,z=redshift,
                         path_to_txt_file=None,fast=fast[i],dilution_factor=10,filters_directory=filters_directory)
                    else:
                        [best_temp, best_radius, best_luminosity, best_coeff] = fit_black_body_flux_filters_mcmc.fit_black_body_flux_filters_mcmc \
                            (Spectrum_right_format,
                             triangle_plot_title=r'$JD-t_{ref}=$' + str(round(j['time'] - explosion_date, 2)), nwalkers=nwalkers,
                             num_steps=num_steps, num_winners=20, already_run_mcmc=already_run_mcmc,
                             already_run_calc_all_chis=False, Temp_prior=np.array([7000, 35000]),
                             Radius_prior=np.array([1e12, 1e15]), initial_conditions=np.array([15000, 8e13]),
                             distance_pc=distance_pc,
                             Ebv=EBV, ndof=None, show_plot=False,
                             output_mcmc=output+'/day_' + str(round(j['time'] - explosion_date, 3)),
                             show_mag_AB=True, z=redshift, path_to_txt_file=None, fast=fast[i], dilution_factor=10)

                else:
                    if fast[i]==False:
                        best_temp=np.genfromtxt(output+'/day_'+str(round(j['time'],3))+'/best_fit_results.txt')[0]
                        best_radius=np.genfromtxt(output+'/day_'+str(round(j['time'],3))+'/best_fit_results.txt')[1]
                        best_luminosity=np.genfromtxt(output+'/day_'+str(round(j['time'],3))+'/best_fit_results.txt')[2]
                    else:
                        best_temp=np.genfromtxt(output+'/day_'+str(round(j['time'],3))+'/best_fit_results_from_median.txt')[0]
                        best_radius=np.genfromtxt(output+'/day_'+str(round(j['time'],3))+'/best_fit_results_from_median.txt')[1]
                        best_luminosity=np.genfromtxt(output+'/day_'+str(round(j['time'],3))+'/best_fit_results_from_median.txt')[2]        #print('I got here 0')

                Best[i, 0] = j['time']
                Best[i, 1] = best_temp
                Best[i, 2] = np.genfromtxt(output+'/day_' + str(round(j['time'], 3)) + '/1sigma.txt',skip_header=1)[0, 0]
                Best[i, 3] =np.genfromtxt(output+'/day_' + str(round(j['time'], 3)) + '/1sigma.txt',skip_header=1)[0, 1]
                Best[i, 4] = best_radius
                Best[i, 5] = np.genfromtxt(output+'/day_' + str(round(j['time'], 3)) + '/1sigma.txt',skip_header=1)[1, 0]
                Best[i, 6] = np.genfromtxt(output+'/day_' + str(round(j['time'], 3)) + '/1sigma.txt',skip_header=1)[1, 1]
                Best[i, 7] = best_luminosity

            else:
                print('******* method chosen: LINEAR FIT *******')
                #print(Spectrum_right_format)
                #pdb.set_trace()
                if already_run_fit[i]==False:
                    [Xi_array, best_temp, index_min, best_coeff, best_radius,best_luminosity]=fit_black_body_flux_filters.fit_black_body_flux_filters(
                        Spectrum_right_format,TempVec=np.logspace(math.log10(3000),math.log10(hitemp),num_iterations),num_temp_iterations=None,
                        distance=distance_pc,uncertainties=Spectrum_right_format[:,3],Ebv=EBV,z=redshift,output_file_path=output+'/day_'+str(round(j['time'],3))
                        ,path_to_txt_file=output+'/day_'+str(round(j['time'],3))+'/best_fit_result.txt',filters_directory=filters_directory)
                else:
                    best_temp=np.genfromtxt(output+'/day_'+str(round(j['time'],3))+'/best_fit_result.txt')[1]
                    best_radius=np.genfromtxt(output+'/day_'+str(round(j['time'],3))+'/best_fit_result.txt')[4]
                    best_luminosity=np.genfromtxt(output+'/day_'+str(round(j['time'],3))+'/best_fit_result.txt')[5]

                Best[i,0]=j['time']
                Best[i,1]=best_temp
                Best[i,2]='0'#np.genfromtxt(output_file_linear+'/day_'+str(round(j['time']-explosion_date,3))+'/1sigma.txt',skip_header=1)[0,0]
                Best[i,3]='0'#np.genfromtxt(output_file_linear+'/day_'+str(round(j['time']-explosion_date,3))+'/1sigma.txt',skip_header=1)[0,1]
                Best[i,4]=best_radius
                Best[i,5]='0'#np.genfromtxt(output_file_linear+'/day_'+str(round(j['time']-explosion_date,3))+'/1sigma.txt',skip_header=1)[1,0]
                Best[i,6] ='0'#np.genfromtxt(output_file_linear+'/day_'+str(round(j['time']-explosion_date,3))+'/1sigma.txt',skip_header=1)[1,1]
                Best[i,7]=best_luminosity
            #pylab.show()


        np.savetxt(output + '/Results.txt', Best,
                       header='JD, best T, lower sigma_T, upper sigma_T,best R, lower sigma_R, upper sigma_R, best L')


    else:
        print('******* You have ran the fit already *******')

        Best = np.genfromtxt(output+'/Results.txt', skip_header=1)
    if show_underlying_plots==True:
        pylab.show()
    else:
        pylab.close('all')
    return Best #JD, best T, lower sigma_T, upper sigma_T,best R, lower sigma_R, upper sigma_R, best L

def plot_T_and_R_in_time(Best,data_compare=None,compare=False,label_comparision=None,output=None):
    print('**********************************')
    print('***** PLOT CALCULATED T AND R ****')
    print('**********************************')
    #if data_compare is None:
    #    print('you need to provide a "data_compare" parameter')
    #    exit()
    #Best_comare=np.genfromtxt('/Users/maayanesoumagnac/PostDoc/projects/2018fif/comparision_with_others/results_fit_sed_mat_PTF13dqy_test/Results.txt')

    if compare==True:
        if os.path.isfile(data_compare)==False:
            print('I do not find any file at the location of the data to compare,',data_compare)
        Best_compare = np.genfromtxt(data_compare)

    #####################################   Plot best T and best R with error bars, compare with other methods ###############
    pylab.figure()
    pylab.plot(Best[:, 0], Best[:, 1], 'ro', label='Best fit')
    if (compare==True) & (label_comparision is not None):
        pylab.plot(Best_compare[:, 0], Best_compare[:, 1], 'bo', label=label_comparision)
        pylab.legend()
    elif (compare==True) & (label_comparision is None):
        pylab.plot(Best_compare[:, 0], Best_compare[:, 1], 'bo',label='comparision file')
        pylab.legend()
    for i, j in enumerate(Best):
        pylab.vlines(Best[i, 0] , Best[i, 2], Best[i, 3], color='red')
    pylab.xlabel(r'time $[days]$')
    pylab.ylabel(r'$T_{BB}\; [K]$')
    pylab.grid()
    #pylab.legend()
    #pylab.xlim(np.min(Best[i, 0]) - 1, np.max(Best[i, 0]) + 1)
    #pylab.xlim(0, np.max(Best[i, 0]) + 1)
    pylab.savefig(output + '/T_bb_evo.png',
                  facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format='png', transparent=False,
                  bbox_inches=None, pad_inches=0.5)
    #pylab.show()

    pylab.figure()
    pylab.plot(Best[:, 0], Best[:, 4], 'go', label='Best fit')
    for i, j in enumerate(Best):
        pylab.vlines(Best[i, 0] , Best[i, 5], Best[i, 6], color='green')
    if (compare==True) & (label_comparision is not None):
        pylab.plot(Best_compare[:, 0], Best_compare[:, 2], 'bo', label=label_comparision)
        pylab.legend()
    elif (compare==True) & (label_comparision is None):
        pylab.plot(Best_compare[:, 0], Best_compare[:, 2], 'bo',label='comparision file')
        pylab.legend()
    pylab.xlabel(r'time $[days]$')
    pylab.ylabel(r'$r_{BB}\; [cm]$')
    pylab.grid()
    #pylab.legend()
    #pylab.xlim(np.min(Best[i, 0])-1,np.max(Best[i, 0])+1)
    pylab.savefig(output+'/r_bb_evo.png',
                  facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format='png', transparent=False,
                  bbox_inches=None, pad_inches=0.5)

    pylab.show()
    '''
    pylab.figure()
    pylab.plot(Best[:,0],Best[:,7],'ro',label='best fit from mat')
    #for i,j in enumerate(Best):
    #    pylab.vlines(Best[i,0],Best[i,2],Best[i,3],color='red')
    #pylab.grid()
    pylab.xlabel(r'time $[days]$')
    pylab.ylabel(r'$L_{BB}$')
    pylab.grid()
    ax=pylab.gca()
    ax.set_xscale("log")
    ax.set_yscale("log")
    #pylab.legend()
    #pylab.savefig(output+'/L_bb_evo_log.png',
    #              facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format='png', transparent=False,
    #              bbox_inches=None, pad_inches=0.5)
    pylab.show()
    '''

def plot_L_in_time(Best,dates_file=None,data_file=None,lower_limit_on_flux=None,error_lum_ran=False,explosion_date=None,output=None,mcmc=False,output_file_interpolation=None,ylim=None,excluded_bands=[],verbose=False):
    ################ L #######
    print('**********************************')
    print('******* PLOT CALCULATED L ********')
    print('**********************************')

    sigma_Boltzmann = 5.670367e-8

    data_dico = \
        read_data_from_file.read_data_into_numpy_array(data_file, header=True, delimiter=',', no_repeat_rows=True)[2]

    if mcmc==True:
        data_dicts = dict()
        for i, j in enumerate(data_dico['filter']):
            data_dicts[j] = dict()
            data_dicts[j]['jd'] = data_dico['jd'][
                (data_dico['filter'] == j) & (data_dico['flux'] > lower_limit_on_flux)]
            data_dicts[j]['mag'] = data_dico['mag'][
                (data_dico['filter'] == j) & (data_dico['flux'] > lower_limit_on_flux)]
            data_dicts[j]['flux'] = data_dico['flux'][
                (data_dico['filter'] == j) & (data_dico['flux'] > lower_limit_on_flux)]
            data_dicts[j]['magerr'] = data_dico['magerr'][
                (data_dico['filter'] == j) & (data_dico['flux'] > lower_limit_on_flux)]
            data_dicts[j]['fluxerr'] = data_dico['fluxerr'][
                (data_dico['filter'] == j) & (data_dico['flux'] > lower_limit_on_flux)]
            data_dicts[j]['absmag'] = data_dico['absmag'][
                (data_dico['filter'] == j) & (data_dico['flux'] > lower_limit_on_flux)]
            data_dicts[j]['absmagerr'] = data_dico['absmagerr'][
                (data_dico['filter'] == j) & (data_dico['flux'] > lower_limit_on_flux)]
            data_dicts[j]['filter'] = data_dico['filter'][
                (data_dico['filter'] == j) & (data_dico['flux'] > lower_limit_on_flux)]
            data_dicts[j]['instr'] = data_dico['instr'][
                (data_dico['filter'] == j) & (data_dico['flux'] > lower_limit_on_flux)]

        #print(data_dicts.keys())

        ########################################### Make the spectrum to fit for each epoch ############################################

        ##### definition of interpolation dates
        if dates_file is None:
            print('you need to provide a dates_file')
            exit()
        interpolation_dates = np.genfromtxt(dates_file)
        condition = dict()
        JD_basis_interpolation = dict()
        already_run_interp_errors = dict()
        jd_flux_fluxerr = dict()
        interp_and_errors_array = dict()
        interp = dict()

        already_run_interp_errors['UVW1'] = True
        already_run_interp_errors['UVW2'] = True
        already_run_interp_errors['UVM2'] = True
        already_run_interp_errors['u_swift'] = True
        already_run_interp_errors['b_swift'] = True
        already_run_interp_errors['v_swift'] = True
        already_run_interp_errors['r_sdss'] = True
        already_run_interp_errors['g_sdss'] = True
        already_run_interp_errors['i_sdss'] = True
        already_run_interp_errors['r_p48'] = True
        already_run_interp_errors['g_p48'] = True
        already_run_interp_errors['z_sdss'] = True

        ### Interpolate all data on the interpolation dates dates ####
        for k in [x for x in data_dicts.keys() if x not in excluded_bands]:
            if verbose==True:
                print('I am looking at the {0} band of the data'.format(k))
                print('the days of {0} are {1}'.format(k, data_dicts[k]['jd'] - explosion_date))
            condition[k] = np.logical_and(interpolation_dates >= np.min(data_dicts[k]['jd']),
                                          interpolation_dates <= np.max(data_dicts[k]['jd']))
            JD_basis_interpolation[k] = interpolation_dates[condition[k]]
            if verbose == True:
                print('JD basis interp. of {0} over the interpolation days are {1}'.format(k, JD_basis_interpolation[
                k] - explosion_date))
                print('and interpolating on {0}'.format(interpolation_dates - explosion_date))
            if len(JD_basis_interpolation[k]) > 0:
                a = np.array(list(zip(data_dicts[k]['jd'], data_dicts[k]['flux'], data_dicts[k]['fluxerr'])))
                jd_flux_fluxerr[k] = a[a[:, 0].argsort()]
                output_path = output_file_interpolation + '/errors_interpolation_results_{0}'.format(k)
                if os.path.exists(output_path):
                    if verbose==True:
                        print(output_path + 'exists')
                else:
                    os.mkdir(output_path)
                #print('already_run_interp_errors', already_run_interp_errors)
                interp_and_errors_array[k] = interpolate_errors.interpolate_errors \
                    (jd_flux_fluxerr[k], JD_basis_interpolation[k],
                     output_path=output_file_interpolation + '/errors_interpolation_results_{0}'.format(k),
                     already_run=already_run_interp_errors[k], show_plot=False,
                     title='{0} on the interpolation dates'.format(
                         k))  # array with days_nuv, interp1d P48, interp P48 with another method, lower limit 1sigma, upper limit 1 sigma
                interp[k] = dict()
                for i, j in enumerate(JD_basis_interpolation[k]):
                    interp[k][str(round(j - explosion_date, 8))] = interp_and_errors_array[k][i, :]
        pylab.close('all')
        Spectra = []

        for i, j in enumerate(interpolation_dates):
            if verbose == True:
                print('********')
                print('i is', i)
                print('********')
            epoch = dict()
            epoch['time'] = round(interpolation_dates[i] - explosion_date, 8)  # jd_flux_fluxerr_UVW2[2,0]#
            # print('epoch[time] is',str(epoch['time']))
            for k in [x for x in data_dicts.keys() if x not in excluded_bands]:
                if len(JD_basis_interpolation[k]) > 0:
                    if str(epoch['time']) in interp[k].keys():
                        epoch[k] = interp[k][str(epoch['time'])][1]  # flux
                        epoch[k + '_err'] = [interp[k][str(epoch['time'])][3], interp[k][str(epoch['time'])][4]]  # fluxerr
            Spectra.append(epoch)

        path_to_errors=Path(output+'/results_errors_lumi')
        if path_to_errors.is_dir() == False:
            os.mkdir(output+'/results_errors_lumi')
        ####################### Calc errors on L ################

        errors_luminosity_early=np.zeros((np.shape(Spectra)[0],3))

        errors_lum_ran=error_lum_ran
        if errors_lum_ran==False:
            for i, j in enumerate(Spectra):
                print('epoch #{0}'.format(i+1))
                #print('file is output+/day_' + str(round(j['time'], 3)) + '/flatchain.txt')
                output_txt_file_path=output+'/results_errors_lumi/day_' + str(round(j['time'], 3))
                output_png_file_path=output_txt_file_path
                my_path=Path(output_txt_file_path)
                if my_path.is_dir()==False:
                    os.mkdir(output_txt_file_path)
                temperatures_early=np.genfromtxt(output+'/day_' + str(round(j['time'], 3))+'/flatchain.txt')[:,0]
                radii_early=np.genfromtxt(output+'/day_' + str(round(j['time'], 3))+'/flatchain.txt')[:,1]
                luminosities_early=energy_conversions.convert_energy(4*math.pi*np.power(radii_early*1e-2,2)*sigma_Boltzmann*np.power(temperatures_early,4),'J','erg')
                #histos de luminosities
                histos=fitter_general.plot_1D_marginalized_distribution(luminosities_early, bests=None, output_pdf_file_path=output_png_file_path, output_txt_file_path=output_txt_file_path,
                                                  parameters_labels=None, number_bins=10000)
        for i, j in enumerate(Spectra):
            errors_luminosity_early[i,0]=round(j['time'],3)
            errors_luminosity_early[i,1]=np.genfromtxt(output+'/results_errors_lumi/day_'+str(round(j['time'],3))+'/1sigma.txt',comments='#')[0]
            errors_luminosity_early[i,2]=np.genfromtxt(output+'/results_errors_lumi/day_'+str(round(j['time'],3))+'/1sigma.txt',comments='#')[1]

        results=np.genfromtxt(output+'/Results.txt')
        Lum=Best[:,7]

        pylab.figure()
        pylab.plot(Best[:,0],Lum,'bo',label='Photometry')
        pylab.xlabel(r'time [days]',fontsize=25)
        pylab.ylabel(r'$\rm{L_{BB}}$ [erg/s]',fontsize=25)
        #pylab.legend(prop={'size': 12})
        ax=pylab.gca()
        #ax.set_yscale("log")
        #pylab.grid()
        pylab.xticks(fontsize=16)
        pylab.yticks(fontsize=16)
        for i,j in enumerate(Spectra):
            pylab.vlines(j['time'],errors_luminosity_early[i,1],errors_luminosity_early[i,2],color='blue')
        pylab.grid(True, which="both")
        if ylim is not None:
            pylab.ylim(ylim[0],ylim[1])
        pylab.tight_layout()
        pylab.savefig(output+'/L_bb_evo.png',
                      facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format='png', transparent=False,
                      bbox_inches=None, pad_inches=0.5)
    else:
        print('Because you did not use mcmc, I am plotting L without errors. To use mcmc change the "mcmc" parameter in params.py and re-run calculate_T_and_R_in_time')
        Lum = Best[:, 7]
        pylab.figure()
        pylab.plot(Best[:, 0], Lum, 'bo', label='Photometry')
        pylab.xlabel(r'time [days]', fontsize=25)
        pylab.ylabel(r'$\rm{L_{BB}}$ [erg/s]', fontsize=25)
        pylab.legend(prop={'size': 12})
        ax = pylab.gca()
        # ax.set_yscale("log")
        # pylab.grid()
        pylab.xticks(fontsize=16)
        pylab.yticks(fontsize=16)
        pylab.grid(True, which="both")
        pylab.tight_layout()
        pylab.savefig(output + '/L_bb_evo.png',
                      facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format='png',
                      transparent=False,
                      bbox_inches=None, pad_inches=0.5)


    plt.show()

def plot_SEDs(Best,already_interpolated=False,data_file=None,lower_limit_on_flux=None,dates_file=None,already_run_interp_errors_from_params=None, number_of_plot=9,redshift=None,distance_modulus=None,explosion_date=None,output=None,
              filters_directory=None,output_file_interpolation=None,EBV=None,ylim=None,excluded_bands=[],verbose=False):

    distance_pc = distances_conversions.DM_to_pc(distance_modulus)

    print('**********************************')
    print('***** MOSAIC PLOT OF SEDs ********')
    print('**********************************')

    #if ylim is not None:
        #print('ylim is not none!')
        #print('ylim is', ylim)
        #pdb.set_trace()

    if dates_file is None:
        print('you need to provide a dates file')
        exit()
    if data_file is None:
        print('you need to provide a data file')
        exit()
    if already_run_interp_errors_from_params is None:
        print('you need to provide a already_run_interp_errors_from_params array')
        exit()

    data_dico = \
        read_data_from_file.read_data_into_numpy_array(data_file, header=True, delimiter=',', no_repeat_rows=True)[2]

    color = dict()
    color['UVW1'] = 'blue'
    color['UVW2'] = 'black'
    color['UVM2'] = 'purple'
    color['u_swift'] = 'orange'
    color['b_swift'] = 'cyan'
    color['v_swift'] = 'magenta'
    color['r_sdss'] = 'red'
    color['i_sdss'] = 'black'
    color['g_sdss'] = 'green'
    color['z_sdss'] = 'grey'
    color['u_sdss'] = 'orange'
    color['r_p48'] = 'red'
    color['g_p48'] = 'green'
    color['i_p48'] = 'black'
    color['u_johnson'] = 'orange'
    color['b_johnson'] = 'cyan'
    color['v_johnson'] = 'magenta'
    color['r_cousin'] = 'red'
    color['i_cousin'] = 'black'
    color['h_2mass'] = 'darkcyan'
    color['j_2mass'] = 'sienna'

    symbol = dict()
    symbol['UVW1'] = 'd'
    symbol['UVW2'] = 'd'
    symbol['UVM2'] = 'd'
    symbol['u_swift'] = 'd'
    symbol['b_swift'] = 'd'
    symbol['v_swift'] = 'd'
    symbol['r_cousin'] = 'v'
    symbol['i_cousin'] = 'v'
    symbol['r_sdss'] = '*'
    symbol['i_sdss'] = '*'
    symbol['g_sdss'] = '*'
    symbol['z_sdss'] = '*'
    symbol['u_sdss'] = '*'
    symbol['r_p48'] = 'o'
    symbol['g_p48'] = 'o'
    symbol['i_p48'] = 'o'
    symbol['j_2mass'] = 'x'
    symbol['h_2mass'] = 'x'
    symbol['u_johnson'] = 's'
    symbol['b_johnson'] = 's'
    symbol['v_johnson'] = 's'
    
    name_bands = dict()
    name_bands['UVW1'] = 'UW1'
    name_bands['UVW2'] = 'UW2'
    name_bands['UVM2'] = 'UM2'
    name_bands['u_swift'] = 'u (UVOT)'
    name_bands['b_swift'] = 'B (UVOT)'
    name_bands['v_swift'] = 'V (UVOT)'
    name_bands['r_cousin'] = 'R'
    name_bands['i_cousin'] = 'I'
    name_bands['r_sdss'] = "r'"
    name_bands['i_sdss'] = "i'"
    name_bands['g_sdss'] = "g'"
    name_bands['z_sdss'] = "z'"
    name_bands['u_sdss'] = "u'"
    name_bands['r_p48'] = 'r (P48)'
    name_bands['g_p48'] = 'g (P48)'
    name_bands['i_p48'] = 'i (P48)'
    name_bands['j_2mass'] = 'J'
    name_bands['h_2mass'] = 'H'
    name_bands['u_johnson'] = 'U'
    name_bands['b_johnson'] = 'B'
    name_bands['v_johnson'] = 'V'
    
    order_bands=dict()
    order_bands['UVW1'] = 3
    order_bands['UVW2'] = 1
    order_bands['UVM2'] = 2
    order_bands['u_swift'] = 4
    order_bands['b_swift'] = 5
    order_bands['v_swift'] = 8
    order_bands['r_cousin'] = 11
    order_bands['i_cousin'] = 14
    order_bands['r_sdss'] = 9
    order_bands['i_sdss'] = 12
    order_bands['g_sdss'] = 6
    order_bands['z_sdss'] = 15
    order_bands['u_sdss'] = 4.5
    order_bands['r_p48'] = 10
    order_bands['g_p48'] = 7
    order_bands['i_p48'] = 13
    order_bands['j_2mass'] = 15
    order_bands['h_2mass'] = 16
    order_bands['u_johnson'] = 4.7
    order_bands['b_johnson'] = 5.5
    order_bands['v_johnson'] = 8.5

    data_dicts = dict()
    for i, j in enumerate(data_dico['filter']):
        data_dicts[j] = dict()
        data_dicts[j]['jd'] = data_dico['jd'][
            (data_dico['filter'] == j) & (data_dico['flux'] > lower_limit_on_flux)]
        #data_dicts[j]['mag'] = data_dico['mag'][
        #    (data_dico['filter'] == j) & (data_dico['flux'] > lower_limit_on_flux)]
        data_dicts[j]['flux'] = data_dico['flux'][
            (data_dico['filter'] == j) & (data_dico['flux'] > lower_limit_on_flux)]
        #data_dicts[j]['magerr'] = data_dico['magerr'][
        #    (data_dico['filter'] == j) & (data_dico['flux'] > lower_limit_on_flux)]
        data_dicts[j]['fluxerr'] = data_dico['fluxerr'][
            (data_dico['filter'] == j) & (data_dico['flux'] > lower_limit_on_flux)]
        #data_dicts[j]['absmag'] = data_dico['absmag'][
        #    (data_dico['filter'] == j) & (data_dico['flux'] > lower_limit_on_flux)]
        #data_dicts[j]['absmagerr'] = data_dico['absmagerr'][
        #    (data_dico['filter'] == j) & (data_dico['flux'] > lower_limit_on_flux)]
        #data_dicts[j]['filter'] = data_dico['filter'][
        #    (data_dico['filter'] == j) & (data_dico['flux'] > lower_limit_on_flux)]
        #data_dicts[j]['instr'] = data_dico['instr'][
        #    (data_dico['filter'] == j) & (data_dico['flux'] > lower_limit_on_flux)]
    if verbose==True:
        print(data_dicts.keys())

    ########################################### Make the spectrum to fit for each epoch ############################################

    ##### definition of interpolation dates

    interpolation_dates = np.genfromtxt(dates_file)
    condition = dict()
    JD_basis_interpolation = dict()
    already_run_interp_errors = dict()
    jd_flux_fluxerr = dict()
    interp_and_errors_array = dict()
    interp = dict()

    for k in [x for x in data_dicts.keys() if x not in excluded_bands]:
        already_run_interp_errors[k] = already_run_interp_errors_from_params[k]

    ### Interpolate all data on the interpolation dates dates ####
    for k in [x for x in data_dicts.keys() if x not in excluded_bands]:

        #print('I am looking at the {0} band of the data'.format(k))
        #print('the days of {0} are {1}'.format(k, data_dicts[k]['jd'] - explosion_date))
        condition[k] = np.logical_and(interpolation_dates >= np.min(data_dicts[k]['jd']),
                                      interpolation_dates <= np.max(data_dicts[k]['jd']))
        JD_basis_interpolation[k] = interpolation_dates[condition[k]]
        #print('JD basis interp. of {0} over the interpolation days are {1}'.format(k, JD_basis_interpolation[
        #    k] - explosion_date))
        #print('and interpolating on {0}'.format(interpolation_dates - explosion_date))
        if already_interpolated==False:
            if len(JD_basis_interpolation[k]) > 0:
                a = np.array(list(zip(data_dicts[k]['jd'], data_dicts[k]['flux'], data_dicts[k]['fluxerr'])))
                jd_flux_fluxerr[k] = a[a[:, 0].argsort()]
                output_path = output_file_interpolation + '/errors_interpolation_results_{0}'.format(k)
                if os.path.exists==False:
                    os.mkdir(output_path)
                interp_and_errors_array[k] = interpolate_errors.interpolate_errors \
                    (jd_flux_fluxerr[k], JD_basis_interpolation[k],
                     output_path=output_file_interpolation + '/errors_interpolation_results_{0}'.format(k),
                     already_run=already_run_interp_errors[k], show_plot=False,
                     title='{0} on the interpolation dates'.format(
                         k))  # array with days_nuv, interp1d P48, interp P48 with another method, lower limit 1sigma, upper limit 1 sigma
                interp[k] = dict()
                for i, j in enumerate(JD_basis_interpolation[k]):
                    interp[k][str(round(j - explosion_date, 8))] = interp_and_errors_array[k][i, :]


        else:
            if len(JD_basis_interpolation[k]) > 0:
                a = np.array(list(zip(data_dicts[k]['jd'], data_dicts[k]['flux'], data_dicts[k]['fluxerr'])))
                jd_flux_fluxerr[k] = a[a[:, 0].argsort()]
                output_path = output_file_interpolation + '/errors_interpolation_results_{0}'.format(k)
                if os.path.exists == False:
                    os.mkdir(output_path)
                interp_and_errors_array[k] = interpolate_errors.interpolate_errors \
                    (jd_flux_fluxerr[k], JD_basis_interpolation[k],
                     output_path=output_file_interpolation + '/errors_interpolation_results_{0}'.format(k),
                     already_run=already_interpolated, show_plot=False,
                     title='{0} on the interpolation dates'.format(
                         k))  # array with days_nuv, interp1d P48, interp P48 with another method, lower limit 1sigma, upper limit 1 sigma
                interp[k] = dict()
                for i, j in enumerate(JD_basis_interpolation[k]):
                    if interp_and_errors_array[k].ndim == 1:
                        interp[k][str(round(j - explosion_date, 8))] = interp_and_errors_array[k][:]
                    else:
                        interp[k][str(round(j - explosion_date, 8))] = interp_and_errors_array[k][i, :]
                pylab.close('all')
                    #interp[k][str(round(j - explosion_date, 8))] = interp_and_errors_array[k][i, :]


    Spectra = []

    for i, j in enumerate(interpolation_dates):
        #print('********')
        #print('i is', i)
        #print('********')
        epoch = dict()
        epoch['time'] = round(interpolation_dates[i] - explosion_date, 8)  # jd_flux_fluxerr_UVW2[2,0]#
        # print('epoch[time] is',str(epoch['time']))
        for k in [x for x in data_dicts.keys() if x not in excluded_bands]:
            if len(JD_basis_interpolation[k]) > 0:
                if str(epoch['time']) in interp[k].keys():
                    epoch[k] = interp[k][str(epoch['time'])][1]  # flux
                    epoch[k + '_err'] = [interp[k][str(epoch['time'])][3], interp[k][str(epoch['time'])][4]]  # fluxerr
        Spectra.append(epoch)

    ########################################### Filters definition ###########################################

    # Filter_vector for all the filters
    Filter_vector = np.empty([21, 2], dtype=object)
    Filter_vector[0] = [str('Swift'), str('UVW1')]
    Filter_vector[1] = [str('Swift'), str('UVW2')]
    Filter_vector[2] = [str('Swift'), str('UVM2')]
    Filter_vector[3] = [str('Swift'), str('u_swift')]
    Filter_vector[4] = [str('Swift'), str('b_swift')]
    Filter_vector[5] = [str('Swift'), str('v_swift')]
    Filter_vector[6] = [str('johnson'), str('u_johnson')]
    Filter_vector[7] = [str('johnson'), str('v_johnson')]
    Filter_vector[8] = [str('johnson'), str('b_johnson')]
    Filter_vector[9] = [str('ztf_p48'), str('r_p48')]
    Filter_vector[10] = [str('ztf_p48'), str('g_p48')]
    Filter_vector[11] = [str('ztf_p48'), str('i_p48')]
    Filter_vector[12] = [str('SDSS'), str('g_sdss')]
    Filter_vector[13] = [str('SDSS'), str('i_sdss')]
    Filter_vector[14] = [str('SDSS'), str('z_sdss')]
    Filter_vector[15] = [str('SDSS'), str('r_sdss')]
    Filter_vector[16] = [str('SDSS'), str('u_sdss')]
    Filter_vector[17] = [str('cousin'), str('i_cousin')]
    Filter_vector[18] = [str('cousin'), str('r_cousin')]
    Filter_vector[19] = [str('2MASS'), str('j_2mass')]
    Filter_vector[20] = [str('2MASS'), str('h_2mass')]

    lib = pyphot.get_library()
    #print("Library contains: ", len(lib), " filters")

    [P,wavelengths_filter]=get_filter.make_filter_object(Filter_vector,filters_directory=filters_directory)
    #print('wavelengths are',wavelengths_filter)

    #n=len([name for name in os.listdir(output)])
    #print('output is',output)

    if number_of_plot==9:
        number=0
        number=len(open(dates_file).readlines())
        #for name in os.listdir(output):
        #    #print('os.path.splitext(name)[1] is ',os.path.splitext(name)[1])
        #    if os.path.splitext(name)[1] not in ['.png','.txt']:
        #        if name[0]=='d':#directories starting with "day"
        #            #print('it is neither .png or .txt')
        #            number=number+1
                #print(name)
            #else:
                #print('it is png or txt')

        print('there are {0} epochs to plot'.format(number))

        indexes=range(number)
        #for i in indexes:
        #    print(i)
        #print(indexes)
        if number>=9:
            a = number // 9
            print('I will take every {0} file from this directory'.format(a))
            indexes = range(number)[0::a]
            fig, axes2d = plt.subplots(nrows=3, ncols=3,
                                       sharex=True, sharey=True,
                                       figsize=(10, 10))
            Spectra2D = np.empty((3, 3), dtype=object)
            for i,j in enumerate( Spectra2D[0:3, 0]):
                Spectra2D[i, 0] = Spectra[indexes[i]]
            for i, j in enumerate(Spectra2D[0:3, 1]):
                Spectra2D[i, 1] = Spectra[indexes[i+3]]
            for i, j in enumerate(Spectra2D[0:3, 2]):
                Spectra2D[i, 2] = Spectra[indexes[i+6]]


            Result_T2D = np.empty((3, 3), dtype=object)
            #for i,j in enumerate(Spectra2)
            Result_T2D[0:3, 0] = Best[indexes[0:3], 1]
            Result_T2D[0:3, 1] = Best[indexes[3:6], 1]
            Result_T2D[0:3, 2] = Best[indexes[6:9], 1]


            Result_R2D = np.empty((3, 3), dtype=object)
            Result_R2D[0:3, 0] = Best[indexes[0:3], 4]
            Result_R2D[0:3, 1] = Best[indexes[3:6], 4]
            Result_R2D[0:3, 2] = Best[indexes[6:9], 4]


            bands_in_legend=[]
            order_in_legend=[]
            line=np.empty((3,3))



            #pylab.show()

            '''
            g[:, 1] = calc_black_body_flux_filters.calc_black_body_flux_filters(self.T, np.arange(1e-7, 3e-6, 1e-9),
                                                                                Filter_vector=None, P_vector=P_vector,
                                                                                Radius=self.r, distance_pc=distance_pc,
                                                                                output_plot=False,
                                                                                show_plots=False, output_txt=False,
                                                                                lib=lib, z=z, Ebv=Ebv, R_ext=R_ext)[:,
                      3]

            '''

            for i, row in enumerate(Spectra2D):
                for k, j in enumerate(row):
                    if verbose==True:
                        print('i is',i)
                        print('k is',k)
                        print('Result_T2D[i, k] is',Result_T2D[i, k])
                        print('Result_R2D[i, k] is',Result_R2D[i, k])
                    #pdb.set_trace()
                    # print(i)
                    #print(j['time'])
                    #print(j.keys())

                    best_fit_full = black_body_flux_density.black_body_flux_density(Result_T2D[i, k],
                                                                                    np.arange(1e-7, 3e-6, 1e-9),
                                                                                    'P',
                                                                                    distance_pc=distance_pc,
                                                                                    Radius=distances_conversions.cm_to_solar_radius(
                                                                                        Result_R2D[i, k]),
                                                                                    redshift=redshift, Ebv=EBV)[2]
                    axes2d[k, i].plot(best_fit_full[:, 0] * 1e10, best_fit_full[:, 1], 'k-')

                    #axes2d[k, i].plot(np.arange(1e-7, 3e-6, 1e-9) * 1e10,
                    #                  black_body_flux_density.black_body_flux_density(Result_T2D[i, k], np.arange(1e-7, 3e-6, 1e-9),
                    #                                                                  'P',
                    #                                                                  distance_pc=distance_pc,
                    #                                                                  Radius=distances_conversions.cm_to_solar_radius(
                    #                                                                      Result_R2D[i, k]), redshift=redshift, Ebv=EBV)[
                    #                      2][:, 1], 'k-')
                    if 'UVW1' in j.keys():
                        axes2d[k, i].plot([wavelengths_filter['UVW1']], [j['UVW1']], color=color['UVW1'], marker=symbol['UVW1'],label='UW1')
                        axes2d[k, i].vlines(wavelengths_filter['UVW1'], j['UVW1_err'][0], j['UVW1_err'][1])
                        if 'UVW1' not in bands_in_legend:
                            bands_in_legend.append('UVW1')
                            order_in_legend.append('UVW1')
                        # axes2d[k,i].errorbar([wavelengths_filter['UVW1']], [j['UVW1']],yerr=[j['UVW1_err']], color='purple',fmt='o', label='UVW1')
                    if 'UVW2' in j.keys():
                        axes2d[k, i].plot([wavelengths_filter['UVW2']], [j['UVW2']], color=color['UVW2'], marker=symbol['UVW2'],label='UW2')  # , label='UVW2')
                        axes2d[k, i].vlines(wavelengths_filter['UVW2'], j['UVW2_err'][0], j['UVW2_err'][1])
                        if 'UVW2' not in bands_in_legend:
                            bands_in_legend.append('UVW2')
                            order_in_legend.append('UVW2')
                        # axes2d[k,i].errorbar([wavelengths_filter['UVW2']], [j['UVW2']],yerr=[j['UVW2_err']], color='k',fmt='o', label='UVW2')
                    if 'UVM2' in j.keys():
                        axes2d[k, i].plot([wavelengths_filter['UVM2']], [j['UVM2']], color=color['UVM2'], marker=symbol['UVM2'],label='UM2')  # , label='UVM2')
                        axes2d[k, i].vlines(wavelengths_filter['UVM2'], j['UVM2_err'][0], j['UVM2_err'][1])
                        if 'UVM2' not in bands_in_legend:
                            bands_in_legend.append('UVM2')
                            order_in_legend.append('UVM2')
                    if 'r_p48' in j.keys():
                        axes2d[k, i].plot([wavelengths_filter['r_p48']], [j['r_p48']], color=color['r_p48'],marker=symbol['r_p48'],label='P48 r')  # ,label='P48 r')
                        axes2d[k, i].vlines(wavelengths_filter['r_p48'], j['r_p48_err'][0], j['r_p48_err'][1])
                        if 'r_p48' not in bands_in_legend:
                            bands_in_legend.append('r_p48')
                            order_in_legend.append('r_p48')
                    if 'g_p48' in j.keys():
                        axes2d[k, i].plot([wavelengths_filter['g_p48']], [j['g_p48']], color=color['g_p48'],marker=symbol['g_p48'],label='P48 g')  # ,label='P48 g')
                        axes2d[k, i].vlines(wavelengths_filter['g_p48'], j['g_p48_err'][0], j['g_p48_err'][1])
                        if 'g_p48' not in bands_in_legend:
                            bands_in_legend.append('g_p48')
                            order_in_legend.append('g_p48')
                    if 'g_sdss' in j.keys():
                        axes2d[k, i].plot([wavelengths_filter['g_sdss']], [j['g_sdss']], color=color['g_sdss'],marker=symbol['g_sdss'],label='P60 g')  # ,label='P60 g')
                        axes2d[k, i].vlines(wavelengths_filter['g_sdss'], j['g_sdss_err'][0], j['g_sdss_err'][1])
                        if 'g_sdss' not in bands_in_legend:
                            bands_in_legend.append('g_sdss')
                            order_in_legend.append('g_sdss')
                    if 'i_sdss' in j.keys():
                        axes2d[k, i].plot([wavelengths_filter['i_sdss']], [j['i_sdss']], color=color['i_sdss'],marker=symbol['i_sdss'],label='P60 i')  # ,label='P60 i')
                        axes2d[k, i].vlines(wavelengths_filter['i_sdss'], j['i_sdss_err'][0], j['i_sdss_err'][1])
                        if 'i_sdss' not in bands_in_legend:
                            bands_in_legend.append('i_sdss')
                            order_in_legend.append('i_sdss')
                    if 'r_sdss' in j.keys():
                        axes2d[k, i].plot([wavelengths_filter['r_sdss']], [j['r_sdss']], color=color['r_sdss'],marker=symbol['r_sdss'],label='P60 r')  # ,label='P60 r')
                        axes2d[k, i].vlines(wavelengths_filter['r_sdss'], j['r_sdss_err'][0], j['r_sdss_err'][1])
                        if 'r_sdss' not in bands_in_legend:
                            bands_in_legend.append('r_sdss')
                            order_in_legend.append('r_sdss')
                    if 'z_sdss' in j.keys():
                        axes2d[k, i].plot([wavelengths_filter['z_sdss']], [j['z_sdss']], color=color['z_sdss'],
                                          marker=symbol['z_sdss'], label='P60 z')  # ,label='P60 r')
                        axes2d[k, i].vlines(wavelengths_filter['z_sdss'], j['z_sdss_err'][0], j['z_sdss_err'][1])
                        if 'z_sdss' not in bands_in_legend:
                            bands_in_legend.append('z_sdss')
                            order_in_legend.append('z_sdss')
                    if 'u_sdss' in j.keys():
                        axes2d[k, i].plot([wavelengths_filter['u_sdss']], [j['u_sdss']], color=color['u_sdss'],
                                          marker=symbol['u_sdss'], label='P60 u')  # ,label='P60 r')
                        axes2d[k, i].vlines(wavelengths_filter['u_sdss'], j['u_sdss_err'][0], j['u_sdss_err'][1])
                        if 'u_sdss' not in bands_in_legend:
                            bands_in_legend.append('u_sdss')
                            order_in_legend.append('u_sdss')
                    if 'u_swift' in j.keys():
                        axes2d[k, i].plot([wavelengths_filter['u_swift']], [j['u_swift']], color=color['u_swift'],marker=symbol['u_swift'],label='u')  # ,label='u')
                        axes2d[k, i].vlines(wavelengths_filter['u_swift'], j['u_swift_err'][0], j['u_swift_err'][1])
                        if 'u_swift' not in bands_in_legend:
                            bands_in_legend.append('u_swift')
                            order_in_legend.append('u_swift')
                    if 'v_swift' in j.keys():
                        axes2d[k, i].plot([wavelengths_filter['v_swift']], [j['v_swift']], color=color['v_swift'],marker=symbol['v_swift'],label='V')  # ,label='v_swift')
                        axes2d[k, i].vlines(wavelengths_filter['v_swift'], j['v_swift_err'][0], j['v_swift_err'][1])
                        if 'v_swift' not in bands_in_legend:
                            bands_in_legend.append('v_swift')
                            order_in_legend.append('v_swift')
                    if 'b_swift' in j.keys():
                        axes2d[k, i].plot([wavelengths_filter['b_swift']], [j['b_swift']], color=color['b_swift'],marker=symbol['b_swift'],label='B')  # ,label='b_swift')
                        axes2d[k, i].vlines(wavelengths_filter['b_swift'], j['b_swift_err'][0], j['b_swift_err'][1])
                        if 'b_swift' not in bands_in_legend:
                            bands_in_legend.append('b_swift')
                            order_in_legend.append('b_swift')
                    if 'b_johnson' in j.keys():
                        axes2d[k, i].plot([wavelengths_filter['b_johnson']], [j['b_johnson']],color=color['b_johnson'],marker=symbol['b_johnson'],label='B')  # ,label='b_swift')
                        axes2d[k, i].vlines(wavelengths_filter['b_johnson'], j['b_johnson_err'][0], j['b_johnson_err'][1])
                        if 'b_johnson' not in bands_in_legend:
                            bands_in_legend.append('b_johnson')
                            order_in_legend.append('b_johnson')
                    if 'v_johnson' in j.keys():
                        axes2d[k, i].plot([wavelengths_filter['v_johnson']], [j['v_johnson']],color=color['v_johnson'],marker=symbol['v_johnson'],label='V')  # ,label='b_swift')
                        axes2d[k, i].vlines(wavelengths_filter['v_johnson'], j['v_johnson_err'][0], j['v_johnson_err'][1])
                        if 'v_johnson' not in bands_in_legend:
                            bands_in_legend.append('v_johnson')
                            order_in_legend.append('v_johnson')
                    if 'u_johnson' in j.keys():
                        axes2d[k, i].plot([wavelengths_filter['u_johnson']], [j['u_johnson']],color=color['u_johnson'],marker=symbol['u_johnson'],label='U')  # ,label='b_swift')
                        axes2d[k, i].vlines(wavelengths_filter['u_johnson'], j['u_johnson_err'][0], j['u_johnson_err'][1])
                        if 'u_johnson' not in bands_in_legend:
                            bands_in_legend.append('u_johnson')
                            order_in_legend.append('u_johnson')
                    if 'i_cousin' in j.keys():
                        axes2d[k, i].plot([wavelengths_filter['i_cousin']], [j['i_cousin']], color=color['i_cousin'],
                                          marker=symbol['i_cousin'], label='P60 r')  # ,label='P60 r')
                        axes2d[k, i].vlines(wavelengths_filter['i_cousin'], j['i_cousin_err'][0], j['i_cousin_err'][1])
                        if 'i_cousin' not in bands_in_legend:
                            bands_in_legend.append('i_cousin')
                            order_in_legend.append('i_cousin')
                    if 'r_cousin' in j.keys():
                        axes2d[k, i].plot([wavelengths_filter['r_cousin']], [j['r_cousin']], color=color['r_cousin'],
                                          marker=symbol['r_cousin'], label='P60 r')  # ,label='P60 r')
                        axes2d[k, i].vlines(wavelengths_filter['r_cousin'], j['r_cousin_err'][0], j['r_cousin_err'][1])
                        if 'r_cousin' not in bands_in_legend:
                            bands_in_legend.append('r_cousin')
                            order_in_legend.append('r_cousin')
                    if 'j_2mass' in j.keys():
                        axes2d[k, i].plot([wavelengths_filter['j_2mass']], [j['j_2mass']], color=color['j_2mass'],
                                          marker=symbol['j_2mass'], label='P60 r')  # ,label='P60 r')
                        axes2d[k, i].vlines(wavelengths_filter['j_2mass'], j['j_2mass_err'][0], j['j_2mass_err'][1])
                        if 'j_2mass' not in bands_in_legend:
                            bands_in_legend.append('j_2mass')
                            order_in_legend.append('j_2mass')
                    if 'h_2mass' in j.keys():
                        axes2d[k, i].plot([wavelengths_filter['h_2mass']], [j['h_2mass']], color=color['h_2mass'],
                                          marker=symbol['h_2mass'], label='P60 r')  # ,label='P60 r')
                        axes2d[k, i].vlines(wavelengths_filter['h_2mass'], j['h_2mass_err'][0], j['h_2mass_err'][1])
                        if 'h_2mass' not in bands_in_legend:
                            bands_in_legend.append('h_2mass')
                            order_in_legend.append('h_2mass')


                    axes2d[k, i].set_xscale("log")#, nonposx='clip')
                    axes2d[k, i].set_yscale("log")#, nonposx='clip')


                    #axes2d[k, i].legend(loc='top right')
                    axes2d[k, i].grid()
                    axes2d[k, i].set_title(r'JD-t$_0$={0}'.format(round(j['time'],2)))
                    #if (ylim1 is not None)&(ylim2 is not None) & (ylim3 is not None):
                    #    if k==0:
                    #       axes2d[k, i].set_ylim(ylim1[0], ylim1[1])
                    #    if k==1:
                    #        axes2d[k, i].set_ylim(ylim2[0], ylim2[1])
                    #    if k==2:
                    #        axes2d[k, i].set_ylim(ylim3[0], ylim3[1])
                    if ylim is not None:
                        axes2d[k, i].set_ylim(ylim[0], ylim[1])
                    #axes2d[k,i].set_ylim(1e-17,1e-15)
                    # axes2d[k,i].set_ylabel('flux $F\; [erg/s/cm^2/\AA]$', fontsize=20)
                    # axes2d[k,i].set_xlabel(r'wavelength [$\AA$]', fontsize=20)
                    # axes2d[k,i].savefig('results_fit_sed_mat/day_'+str(round(j['time'],3))+'/spectrum.png',
                    #              facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format='png', transparent=False,
                    #              bbox_inches=None, pad_inches=0.5)
                    # axes2d[k,i].show()
            #ax = pylab.gca()
            #box = ax.get_position()
            #ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
            #axes2d[1, 2 ].legend(loc='lower right',bbox_to_anchor=(1.2, 0.5))
            #art=[]

            from matplotlib.lines import Line2D
            #handles, labels = plt.gca().get_legend_handles_labels()
            #print(handles)
            #fig.legend(handles, labels, loc='lower left')
            #print(bands_in_legend)
            #print(order_in_legend)
            #pdb.set_trace()
            labels=[name_bands[i] for i in bands_in_legend]
            lines = [Line2D([0], [0], color=color[i],marker=symbol[i] ) for i in bands_in_legend]
            fig.legend(lines,labels,loc='lower left')
            #axes2d[1, 1].legend(loc='lower center', bbox_to_anchor = (0.5, -0.3),ncol=8)
            #art.append(lgd)
            axes2d[2,1].set_xlabel(r'wavelength [$\AA$]', fontsize=20)
            axes2d[1,0].set_ylabel(r'flux $F\; \rm [erg/s/cm^2/\AA]$', fontsize=20)
            #plt.tight_layout()
            fig.subplots_adjust(left=0.2)

            pylab.savefig(output + '/2D_SEDs_9.png',
                              facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format='png', transparent=False,
                              bbox_inches=None, pad_inches=1.5)
            plt.show()
        else:
            #print('number is',number)
            fig, axes2d = plt.subplots(nrows=number, ncols=1,
                                       sharex=True, sharey=True,
                                       figsize=(10, 10))
            #print('the shape of axes2d is '.format(axes2d))
            #pdb.set_trace()
            Spectra2D = np.empty((number, 1), dtype=object)
            #print(np.shape(Spectra2D))
            #print(np.shape(Spectra))
            for i, j in enumerate(Spectra2D):
                Spectra2D[i] = Spectra[i]

            # Result_T2D = np.empty((number, 1), dtype=object)
            # for i,j in enumerate(Spectra2)
            Result_T2D = Best[:, 1]
            Result_R2D = Best[:, 4]
            bands_in_legend = []
            #line = np.empty((3, 3))
            for i, row in enumerate(Spectra2D):
                #print('ble')
                #pdb.set_trace()
                for k, j in enumerate(row):
                    #print(i)
                    # print(j['time'])
                    # print(j.keys())
                    #axes2d[i].plot(np.arange(1e-7, 3e-6, 1e-9) * 1e10,
                    #                  black_body_flux_density.black_body_flux_density(Result_T2D[i],
                    #                                                                  np.arange(1e-7, 3e-6, 1e-9),
                    #                                                                  'P',
                    #                                                                  distance_pc=distance_pc,
                    #                                                                  Radius=distances_conversions.cm_to_solar_radius(
                    #                                                                      Result_R2D[i]),
                    #                                                                  redshift=redshift, Ebv=EBV)[
                    #                      2][:, 1], 'k-')
                    best_fit_full = black_body_flux_density.black_body_flux_density(Result_T2D[i],
                                                                                    np.arange(1e-7, 3e-6, 1e-9),
                                                                                    'P',
                                                                                    distance_pc=distance_pc,
                                                                                    Radius=distances_conversions.cm_to_solar_radius(
                                                                                        Result_R2D[i]),
                                                                                    redshift=redshift, Ebv=EBV)[2]
                    axes2d[i].plot(best_fit_full[:, 0] * 1e10, best_fit_full[:, 1], 'k-')

                    if 'UVW1' in j.keys():
                        axes2d[i].plot([wavelengths_filter['UVW1']], [j['UVW1']], color=color['UVW1'],
                                          marker=symbol['UVW1'], label='UW1')
                        axes2d[i].vlines(wavelengths_filter['UVW1'], j['UVW1_err'][0], j['UVW1_err'][1])
                        if 'UVW1' not in bands_in_legend:
                            bands_in_legend.append('UVW1')
                            # axes2d[k,i].errorbar([wavelengths_filter['UVW1']], [j['UVW1']],yerr=[j['UVW1_err']], color='purple',fmt='o', label='UVW1')
                    if 'UVW2' in j.keys():
                        axes2d[i].plot([wavelengths_filter['UVW2']], [j['UVW2']], color=color['UVW2'],
                                          marker=symbol['UVW2'], label='UW2')  # , label='UVW2')
                        axes2d[i].vlines(wavelengths_filter['UVW2'], j['UVW2_err'][0], j['UVW2_err'][1])
                        if 'UVW2' not in bands_in_legend:
                            bands_in_legend.append('UVW2')
                            # axes2d[k,i].errorbar([wavelengths_filter['UVW2']], [j['UVW2']],yerr=[j['UVW2_err']], color='k',fmt='o', label='UVW2')
                    if 'UVM2' in j.keys():
                        axes2d[i].plot([wavelengths_filter['UVM2']], [j['UVM2']], color=color['UVM2'],
                                          marker=symbol['UVM2'], label='UM2')  # , label='UVM2')
                        axes2d[i].vlines(wavelengths_filter['UVM2'], j['UVM2_err'][0], j['UVM2_err'][1])
                        if 'UVM2' not in bands_in_legend:
                            bands_in_legend.append('UVM2')
                    if 'r_p48' in j.keys():
                        axes2d[i].plot([wavelengths_filter['r_p48']], [j['r_p48']], color=color['r_p48'],
                                          marker=symbol['r_p48'], label='P48 r')  # ,label='P48 r')
                        axes2d[i].vlines(wavelengths_filter['r_p48'], j['r_p48_err'][0], j['r_p48_err'][1])
                        if 'r_p48' not in bands_in_legend:
                            bands_in_legend.append('r_p48')
                    if 'g_p48' in j.keys():
                        axes2d[i].plot([wavelengths_filter['g_p48']], [j['g_p48']], color=color['g_p48'],
                                          marker=symbol['g_p48'], label='P48 g')  # ,label='P48 g')
                        axes2d[i].vlines(wavelengths_filter['g_p48'], j['g_p48_err'][0], j['g_p48_err'][1])
                        if 'g_p48' not in bands_in_legend:
                            bands_in_legend.append('g_p48')
                    if 'i_p48' in j.keys():
                        axes2d[i].plot([wavelengths_filter['i_p48']], [j['i_p48']], color=color['i_p48'],
                                          marker=symbol['i_p48'], label='P48 i')  # ,label='P48 g')
                        axes2d[i].vlines(wavelengths_filter['i_p48'], j['i_p48_err'][0], j['i_p48_err'][1])
                        if 'i_p48' not in bands_in_legend:
                            bands_in_legend.append('i_p48')
                    if 'g_sdss' in j.keys():
                        axes2d[i].plot([wavelengths_filter['g_sdss']], [j['g_sdss']], color=color['g_sdss'],
                                          marker=symbol['g_sdss'], label='P60 g')  # ,label='P60 g')
                        axes2d[i].vlines(wavelengths_filter['g_sdss'], j['g_sdss_err'][0], j['g_sdss_err'][1])
                        if 'g_sdss' not in bands_in_legend:
                            bands_in_legend.append('g_sdss')
                    if 'i_sdss' in j.keys():
                        axes2d[i].plot([wavelengths_filter['i_sdss']], [j['i_sdss']], color=color['i_sdss'],
                                          marker=symbol['i_sdss'], label='P60 i')  # ,label='P60 i')
                        axes2d[i].vlines(wavelengths_filter['i_sdss'], j['i_sdss_err'][0], j['i_sdss_err'][1])
                        if 'i_sdss' not in bands_in_legend:
                            bands_in_legend.append('i_sdss')
                    if 'r_sdss' in j.keys():
                        axes2d[i].plot([wavelengths_filter['r_sdss']], [j['r_sdss']], color=color['r_sdss'],
                                          marker=symbol['r_sdss'], label='P60 r')  # ,label='P60 r')
                        axes2d[i].vlines(wavelengths_filter['r_sdss'], j['r_sdss_err'][0], j['r_sdss_err'][1])
                        if 'r_sdss' not in bands_in_legend:
                            bands_in_legend.append('r_sdss')
                    if 'z_sdss' in j.keys():
                        axes2d[i].plot([wavelengths_filter['z_sdss']], [j['z_sdss']], color=color['z_sdss'],
                                          marker=symbol['z_sdss'], label='P60 z')  # ,label='P60 r')
                        axes2d[i].vlines(wavelengths_filter['z_sdss'], j['z_sdss_err'][0], j['z_sdss_err'][1])
                        if 'z_sdss' not in bands_in_legend:
                            bands_in_legend.append('z_sdss')
                    if 'u_sdss' in j.keys():
                        axes2d[i].plot([wavelengths_filter['u_sdss']], [j['u_sdss']], color=color['u_sdss'],
                                          marker=symbol['u_sdss'], label='P60 u')  # ,label='P60 r')
                        axes2d[i].vlines(wavelengths_filter['u_sdss'], j['u_sdss_err'][0], j['u_sdss_err'][1])
                        if 'u_sdss' not in bands_in_legend:
                            bands_in_legend.append('u_sdss')
                    if 'u_swift' in j.keys():
                        axes2d[i].plot([wavelengths_filter['u_swift']], [j['u_swift']], color=color['u_swift'],
                                          marker=symbol['u_swift'], label='u')  # ,label='u')
                        axes2d[i].vlines(wavelengths_filter['u_swift'], j['u_swift_err'][0], j['u_swift_err'][1])
                        if 'u_swift' not in bands_in_legend:
                            bands_in_legend.append('u_swift')
                    if 'v_swift' in j.keys():
                        axes2d[i].plot([wavelengths_filter['v_swift']], [j['v_swift']], color=color['v_swift'],
                                          marker=symbol['v_swift'], label='V')  # ,label='v_swift')
                        axes2d[i].vlines(wavelengths_filter['v_swift'], j['v_swift_err'][0], j['v_swift_err'][1])
                        if 'v_swift' not in bands_in_legend:
                            bands_in_legend.append('v_swift')
                    if 'b_swift' in j.keys():
                        axes2d[i].plot([wavelengths_filter['b_swift']], [j['b_swift']], color=color['b_swift'],
                                          marker=symbol['b_swift'], label='B')  # ,label='b_swift')
                        axes2d[i].vlines(wavelengths_filter['b_swift'], j['b_swift_err'][0], j['b_swift_err'][1])
                        if 'b_swift' not in bands_in_legend:
                            bands_in_legend.append('b_swift')
                    if 'b_johnson' in j.keys():
                        axes2d[i].plot([wavelengths_filter['b_johnson']], [j['b_johnson']], color=color['b_johnson'],
                                          marker=symbol['b_johnson'], label='B')  # ,label='b_swift')
                        axes2d[i].vlines(wavelengths_filter['b_johnson'], j['b_johnson_err'][0],
                                            j['b_johnson_err'][1])
                        if 'b_johnson' not in bands_in_legend:
                            bands_in_legend.append('b_johnson')
                    if 'v_johnson' in j.keys():
                        axes2d[i].plot([wavelengths_filter['v_johnson']], [j['v_johnson']], color=color['v_johnson'],
                                          marker=symbol['v_johnson'], label='V')  # ,label='b_swift')
                        axes2d[i].vlines(wavelengths_filter['v_johnson'], j['v_johnson_err'][0],
                                            j['v_johnson_err'][1])
                        if 'v_johnson' not in bands_in_legend:
                            bands_in_legend.append('v_johnson')
                    if 'u_johnson' in j.keys():
                        axes2d[i].plot([wavelengths_filter['u_johnson']], [j['u_johnson']], color=color['u_johnson'],
                                          marker=symbol['u_johnson'], label='U')  # ,label='b_swift')
                        axes2d[i].vlines(wavelengths_filter['u_johnson'], j['u_johnson_err'][0],
                                            j['u_johnson_err'][1])
                        if 'u_johnson' not in bands_in_legend:
                            bands_in_legend.append('u_johnson')
                    if 'i_cousin' in j.keys():
                        axes2d[i].plot([wavelengths_filter['i_cousin']], [j['i_cousin']], color=color['i_cousin'],
                                          marker=symbol['i_cousin'], label='P60 r')  # ,label='P60 r')
                        axes2d[i].vlines(wavelengths_filter['i_cousin'], j['i_cousin_err'][0], j['i_cousin_err'][1])
                        if 'i_cousin' not in bands_in_legend:
                            bands_in_legend.append('i_cousin')
                    if 'r_cousin' in j.keys():
                        axes2d[i].plot([wavelengths_filter['r_cousin']], [j['r_cousin']], color=color['r_cousin'],
                                          marker=symbol['r_cousin'], label='P60 r')  # ,label='P60 r')
                        axes2d[i].vlines(wavelengths_filter['r_cousin'], j['r_cousin_err'][0], j['r_cousin_err'][1])
                        if 'r_cousin' not in bands_in_legend:
                            bands_in_legend.append('r_cousin')
                    if 'j_2mass' in j.keys():
                        axes2d[i].plot([wavelengths_filter['j_2mass']], [j['j_2mass']], color=color['j_2mass'],
                                          marker=symbol['j_2mass'], label='P60 r')  # ,label='P60 r')
                        axes2d[i].vlines(wavelengths_filter['j_2mass'], j['j_2mass_err'][0], j['j_2mass_err'][1])
                        if 'j_2mass' not in bands_in_legend:
                            bands_in_legend.append('j_2mass')
                    if 'h_2mass' in j.keys():
                        axes2d[i].plot([wavelengths_filter['h_2mass']], [j['h_2mass']], color=color['h_2mass'],
                                          marker=symbol['h_2mass'], label='P60 r')  # ,label='P60 r')
                        axes2d[i].vlines(wavelengths_filter['h_2mass'], j['h_2mass_err'][0], j['h_2mass_err'][1])
                        if 'h_2mass' not in bands_in_legend:
                            bands_in_legend.append('h_2mass')


                    axes2d[i].set_xscale("log")#, nonposx='clip')
                    axes2d[i].set_yscale("log")#, nonposx='clip')


                    # axes2d[k, i].legend(loc='top right')
                    axes2d[i].grid()
                    axes2d[i].set_title(r'JD-t$_0$={0}'.format(round(j['time'], 2)))
                    if ylim is not None:
                        print('ylim is not none!')
                        print('ylim is',ylim)
                        axes2d[i].set_ylim(ylim[0],ylim[1])
                    # axes2d[k,i].set_ylabel('flux $F\; [erg/s/cm^2/\AA]$', fontsize=20)
                    # axes2d[k,i].set_xlabel(r'wavelength [$\AA$]', fontsize=20)
                    # axes2d[k,i].savefig('results_fit_sed_mat/day_'+str(round(j['time'],3))+'/spectrum.png',
                    #              facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format='png', transparent=False,
                    #              bbox_inches=None, pad_inches=0.5)
                    # axes2d[k,i].show()
            # ax = pylab.gca()
            # box = ax.get_position()
            # ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
            # axes2d[1, 2 ].legend(loc='lower right',bbox_to_anchor=(1.2, 0.5))
            # art=[]

            from matplotlib.lines import Line2D
            # handles, labels = plt.gca().get_legend_handles_labels()
            # print(handles)
            # fig.legend(handles, labels, loc='lower left')
            labels = [name_bands[i] for i in bands_in_legend]
            lines = [Line2D([0], [0], color=color[i], marker=symbol[i]) for i in bands_in_legend]
            fig.legend(lines, labels, loc='lower left')
            # axes2d[1, 1].legend(loc='lower center', bbox_to_anchor = (0.5, -0.3),ncol=8)
            # art.append(lgd)
            axes2d[number-1].set_xlabel(r'wavelength [$\AA$]', fontsize=20)
            axes2d[1].set_ylabel(r'flux $F\; \rm [erg/s/cm^2/\AA]$', fontsize=20)
            # plt.tight_layout()
            fig.subplots_adjust(left=0.2)
            pylab.savefig(output + '/2D_SEDs_9.png',facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format='png',
                          transparent=False,
                          bbox_inches=None, pad_inches=1.5)
            plt.show()

    elif number_of_plot==16:
        number=0
        for name in os.listdir(output):
            #print('os.path.splitext(name)[1] is ',os.path.splitext(name)[1])
            if os.path.splitext(name)[1] not in ['.png','.txt']:
                if name[0]=='d':#directories starting with "day"
                    #print('it is neither .png or .txt')
                    number=number+1
                #print(name)
            #else:
                #print('it is png or txt')

        print('there are {0} directories in the output file'.format(number))

        indexes=range(number)
        #for i in indexes:
        #    print(i)
        #print(indexes)
        if number>=16:
            a = number // 16
            print('I will take every {0} file from this directory'.format(a))
            indexes = range(number)[0::a]
            fig, axes2d = plt.subplots(nrows=4, ncols=4,
                                       sharex=True, sharey=True,
                                       figsize=(10, 10))
            Spectra2D = np.empty((4, 4), dtype=object)
            for i,j in enumerate( Spectra2D[0:4, 0]):
                Spectra2D[i, 0] = Spectra[indexes[i]]
            for i, j in enumerate(Spectra2D[0:4, 1]):
                Spectra2D[i, 1] = Spectra[indexes[i+4]]
            for i, j in enumerate(Spectra2D[0:4, 2]):
                Spectra2D[i, 2] = Spectra[indexes[i+8]]
            for i, j in enumerate(Spectra2D[0:4, 3]):
                Spectra2D[i, 3] = Spectra[indexes[i+12]]


            Result_T2D = np.empty((4, 4), dtype=object)
            #for i,j in enumerate(Spectra2)
            Result_T2D[0:4, 0] = Best[indexes[0:4], 1]
            Result_T2D[0:4, 1] = Best[indexes[4:8], 1]
            Result_T2D[0:4, 2] = Best[indexes[8:12], 1]
            Result_T2D[0:4, 3] = Best[indexes[12:16], 1]

            Result_R2D = np.empty((4, 4), dtype=object)
            Result_R2D[0:4, 0] = Best[indexes[0:4], 4]
            Result_R2D[0:4, 1] = Best[indexes[4:8], 4]
            Result_R2D[0:4, 2] = Best[indexes[8:12], 4]
            Result_R2D[0:4, 3] = Best[indexes[12:16], 4]
            bands_in_legend=[]
            order_in_legend=[]
            line=np.empty((4,4))

            for i, row in enumerate(Spectra2D):
                for k, j in enumerate(row):
                    # print(i)
                    #print(j['time'])
                    #print(j.keys())
                    best_fit_full = black_body_flux_density.black_body_flux_density(Result_T2D[i, k],
                                                                                    np.arange(1e-7, 3e-6, 1e-9),
                                                                                    'P',
                                                                                    distance_pc=distance_pc,
                                                                                    Radius=distances_conversions.cm_to_solar_radius(
                                                                                        Result_R2D[i, k]),
                                                                                    redshift=redshift, Ebv=EBV)[2]
                    axes2d[k, i].plot(best_fit_full[:, 0] * 1e10, best_fit_full[:, 1], 'k-')
                    #axes2d[k, i].plot(np.arange(1e-7, 3e-6, 1e-9) * 1e10,
                    #                  black_body_flux_density.black_body_flux_density(Result_T2D[i, k], np.arange(1e-7, 3e-6, 1e-9),
                    #                                                                  'P',
                    #                                                                  distance_pc=distance_pc,
                    #                                                                  Radius=distances_conversions.cm_to_solar_radius(
                    #                                                                      Result_R2D[i, k]), redshift=redshift, Ebv=EBV)[
                    #                      2][:, 1], 'k-')
                    if 'UVW1' in j.keys():
                        axes2d[k, i].plot([wavelengths_filter['UVW1']], [j['UVW1']], color=color['UVW1'], marker=symbol['UVW1'],label='UW1')
                        axes2d[k, i].vlines(wavelengths_filter['UVW1'], j['UVW1_err'][0], j['UVW1_err'][1])
                        if 'UVW1' not in bands_in_legend:
                            bands_in_legend.append('UVW1')
                            order_in_legend.append(order_bands['UVW1'])
                        # axes2d[k,i].errorbar([wavelengths_filter['UVW1']], [j['UVW1']],yerr=[j['UVW1_err']], color='purple',fmt='o', label='UVW1')
                    if 'UVW2' in j.keys():
                        axes2d[k, i].plot([wavelengths_filter['UVW2']], [j['UVW2']], color=color['UVW2'], marker=symbol['UVW2'],label='UW2')  # , label='UVW2')
                        axes2d[k, i].vlines(wavelengths_filter['UVW2'], j['UVW2_err'][0], j['UVW2_err'][1])
                        if 'UVW2' not in bands_in_legend:
                            bands_in_legend.append('UVW2')
                            order_in_legend.append(order_bands['UVW2'])
                        # axes2d[k,i].errorbar([wavelengths_filter['UVW2']], [j['UVW2']],yerr=[j['UVW2_err']], color='k',fmt='o', label='UVW2')
                    if 'UVM2' in j.keys():
                        axes2d[k, i].plot([wavelengths_filter['UVM2']], [j['UVM2']], color=color['UVM2'], marker=symbol['UVM2'],label='UM2')  # , label='UVM2')
                        axes2d[k, i].vlines(wavelengths_filter['UVM2'], j['UVM2_err'][0], j['UVM2_err'][1])
                        if 'UVM2' not in bands_in_legend:
                            bands_in_legend.append('UVM2')
                            order_in_legend.append(order_bands['UVM2'])
                    if 'r_p48' in j.keys():
                        axes2d[k, i].plot([wavelengths_filter['r_p48']], [j['r_p48']], color=color['r_p48'],marker=symbol['r_p48'],label='P48 r')  # ,label='P48 r')
                        axes2d[k, i].vlines(wavelengths_filter['r_p48'], j['r_p48_err'][0], j['r_p48_err'][1])
                        if 'r_p48' not in bands_in_legend:
                            bands_in_legend.append('r_p48')
                            order_in_legend.append(order_bands['r_p48'])
                    if 'g_p48' in j.keys():
                        axes2d[k, i].plot([wavelengths_filter['g_p48']], [j['g_p48']], color=color['g_p48'],marker=symbol['g_p48'],label='P48 g')  # ,label='P48 g')
                        axes2d[k, i].vlines(wavelengths_filter['g_p48'], j['g_p48_err'][0], j['g_p48_err'][1])
                        if 'g_p48' not in bands_in_legend:
                            bands_in_legend.append('g_p48')
                            order_in_legend.append(order_bands['g_p48'])
                    if 'i_p48' in j.keys():
                        axes2d[k, i].plot([wavelengths_filter['i_p48']], [j['i_p48']], color=color['i_p48'],marker=symbol['i_p48'],label='P48 i')  # ,label='P48 g')
                        axes2d[k, i].vlines(wavelengths_filter['i_p48'], j['i_p48_err'][0], j['i_p48_err'][1])
                        if 'i_p48' not in bands_in_legend:
                            bands_in_legend.append('i_p48')
                            order_in_legend.append(order_bands['i_p48'])
                    if 'g_sdss' in j.keys():
                        axes2d[k, i].plot([wavelengths_filter['g_sdss']], [j['g_sdss']], color=color['g_sdss'],marker=symbol['g_sdss'],label='P60 g')  # ,label='P60 g')
                        axes2d[k, i].vlines(wavelengths_filter['g_sdss'], j['g_sdss_err'][0], j['g_sdss_err'][1])
                        if 'g_sdss' not in bands_in_legend:
                            bands_in_legend.append('g_sdss')
                            order_in_legend.append(order_bands['g_sdss'])
                    if 'i_sdss' in j.keys():
                        axes2d[k, i].plot([wavelengths_filter['i_sdss']], [j['i_sdss']], color=color['i_sdss'],marker=symbol['i_sdss'],label='P60 i')  # ,label='P60 i')
                        axes2d[k, i].vlines(wavelengths_filter['i_sdss'], j['i_sdss_err'][0], j['i_sdss_err'][1])
                        if 'i_sdss' not in bands_in_legend:
                            bands_in_legend.append('i_sdss')
                            order_in_legend.append(order_bands['i_sdss'])
                    if 'r_sdss' in j.keys():
                        axes2d[k, i].plot([wavelengths_filter['r_sdss']], [j['r_sdss']], color=color['r_sdss'],marker=symbol['r_sdss'],label='P60 r')  # ,label='P60 r')
                        axes2d[k, i].vlines(wavelengths_filter['r_sdss'], j['r_sdss_err'][0], j['r_sdss_err'][1])
                        if 'r_sdss' not in bands_in_legend:
                            bands_in_legend.append('r_sdss')
                            order_in_legend.append(order_bands['r_sdss'])
                    if 'z_sdss' in j.keys():
                        axes2d[k, i].plot([wavelengths_filter['z_sdss']], [j['z_sdss']], color=color['z_sdss'],
                                          marker=symbol['z_sdss'], label='P60 z')  # ,label='P60 r')
                        axes2d[k, i].vlines(wavelengths_filter['z_sdss'], j['z_sdss_err'][0], j['z_sdss_err'][1])
                        if 'z_sdss' not in bands_in_legend:
                            bands_in_legend.append('z_sdss')
                            order_in_legend.append(order_bands['z_sdss'])
                    if 'u_sdss' in j.keys():
                        axes2d[k, i].plot([wavelengths_filter['u_sdss']], [j['u_sdss']], color=color['u_sdss'],
                                          marker=symbol['u_sdss'], label='sdss u')  # ,label='P60 r')
                        axes2d[k, i].vlines(wavelengths_filter['u_sdss'], j['u_sdss_err'][0], j['u_sdss_err'][1])
                        if 'u_sdss' not in bands_in_legend:
                            bands_in_legend.append('u_sdss')
                            order_in_legend.append(order_bands['u_sdss'])
                    if 'u_swift' in j.keys():
                        axes2d[k, i].plot([wavelengths_filter['u_swift']], [j['u_swift']], color=color['u_swift'],marker=symbol['u_swift'],label='u')  # ,label='u')
                        axes2d[k, i].vlines(wavelengths_filter['u_swift'], j['u_swift_err'][0], j['u_swift_err'][1])
                        if 'u_swift' not in bands_in_legend:
                            bands_in_legend.append('u_swift')
                            order_in_legend.append(order_bands['u_swift'])
                    if 'v_swift' in j.keys():
                        axes2d[k, i].plot([wavelengths_filter['v_swift']], [j['v_swift']], color=color['v_swift'],marker=symbol['v_swift'],label='V')  # ,label='v_swift')
                        axes2d[k, i].vlines(wavelengths_filter['v_swift'], j['v_swift_err'][0], j['v_swift_err'][1])
                        if 'v_swift' not in bands_in_legend:
                            bands_in_legend.append('v_swift')
                            order_in_legend.append(order_bands['v_swift'])
                    if 'b_swift' in j.keys():
                        axes2d[k, i].plot([wavelengths_filter['b_swift']], [j['b_swift']], color=color['b_swift'],marker=symbol['b_swift'],label='B')  # ,label='b_swift')
                        axes2d[k, i].vlines(wavelengths_filter['b_swift'], j['b_swift_err'][0], j['b_swift_err'][1])
                        if 'b_swift' not in bands_in_legend:
                            bands_in_legend.append('b_swift')
                            order_in_legend.append(order_bands['b_swift'])
                    if 'b_johnson' in j.keys():
                        axes2d[k, i].plot([wavelengths_filter['b_johnson']], [j['b_johnson']],color=color['b_johnson'],marker=symbol['b_johnson'],label='B')  # ,label='b_swift')
                        axes2d[k, i].vlines(wavelengths_filter['b_johnson'], j['b_johnson_err'][0], j['b_johnson_err'][1])
                        if 'b_johnson' not in bands_in_legend:
                            bands_in_legend.append('b_johnson')
                            order_in_legend.append(order_bands['b_johnson'])
                    if 'v_johnson' in j.keys():
                        axes2d[k, i].plot([wavelengths_filter['v_johnson']], [j['v_johnson']],color=color['v_johnson'],marker=symbol['v_johnson'],label='V')  # ,label='b_swift')
                        axes2d[k, i].vlines(wavelengths_filter['v_johnson'], j['v_johnson_err'][0], j['v_johnson_err'][1])
                        if 'v_johnson' not in bands_in_legend:
                            bands_in_legend.append('v_johnson')
                            order_in_legend.append(order_bands['v_johnson'])
                    if 'u_johnson' in j.keys():
                        axes2d[k, i].plot([wavelengths_filter['u_johnson']], [j['u_johnson']],color=color['u_johnson'],marker=symbol['u_johnson'],label='U')  # ,label='b_swift')
                        axes2d[k, i].vlines(wavelengths_filter['u_johnson'], j['u_johnson_err'][0], j['u_johnson_err'][1])
                        if 'u_johnson' not in bands_in_legend:
                            bands_in_legend.append('u_johnson')
                            order_in_legend.append(order_bands['u_johnson'])
                    if 'i_cousin' in j.keys():
                        axes2d[k, i].plot([wavelengths_filter['i_cousin']], [j['i_cousin']], color=color['i_cousin'],
                                          marker=symbol['i_cousin'], label='P60 r')  # ,label='P60 r')
                        axes2d[k, i].vlines(wavelengths_filter['i_cousin'], j['i_cousin_err'][0], j['i_cousin_err'][1])
                        if 'i_cousin' not in bands_in_legend:
                            bands_in_legend.append('i_cousin')
                            order_in_legend.append(order_bands['i_cousin'])
                    if 'r_cousin' in j.keys():
                        axes2d[k, i].plot([wavelengths_filter['r_cousin']], [j['r_cousin']], color=color['r_cousin'],
                                          marker=symbol['r_cousin'], label='P60 r')  # ,label='P60 r')
                        axes2d[k, i].vlines(wavelengths_filter['r_cousin'], j['r_cousin_err'][0], j['r_cousin_err'][1])
                        if 'r_cousin' not in bands_in_legend:
                            bands_in_legend.append('r_cousin')
                            order_in_legend.append(order_bands['r_cousin'])
                    if 'j_2mass' in j.keys():
                        axes2d[k, i].plot([wavelengths_filter['j_2mass']], [j['j_2mass']], color=color['j_2mass'],
                                          marker=symbol['j_2mass'], label='P60 r')  # ,label='P60 r')
                        axes2d[k, i].vlines(wavelengths_filter['j_2mass'], j['j_2mass_err'][0], j['j_2mass_err'][1])
                        if 'j_2mass' not in bands_in_legend:
                            bands_in_legend.append('j_2mass')
                            order_in_legend.append(order_bands['j_2mass'])
                    if 'h_2mass' in j.keys():
                        axes2d[k, i].plot([wavelengths_filter['h_2mass']], [j['h_2mass']], color=color['h_2mass'],
                                          marker=symbol['h_2mass'], label='P60 r')  # ,label='P60 r')
                        axes2d[k, i].vlines(wavelengths_filter['h_2mass'], j['h_2mass_err'][0], j['h_2mass_err'][1])
                        if 'h_2mass' not in bands_in_legend:
                            bands_in_legend.append('h_2mass')
                            order_in_legend.append(order_bands['h_2mass'])


                    axes2d[k, i].set_xscale("log")#, nonposx='clip')
                    axes2d[k, i].set_yscale("log")#, nonposx='clip')


                    #axes2d[k, i].legend(loc='top right')
                    axes2d[k, i].grid()
                    axes2d[k, i].set_title(r'JD-t$_0$={0}'.format(round(j['time'],2)))
                    #if (ylim1 is not None)&(ylim2 is not None) & (ylim3 is not None):
                    #    if k==0:
                    #       axes2d[k, i].set_ylim(ylim1[0], ylim1[1])
                    #    if k==1:
                    #        axes2d[k, i].set_ylim(ylim2[0], ylim2[1])
                    #    if k==2:
                    #        axes2d[k, i].set_ylim(ylim3[0], ylim3[1])
                    if ylim is not None:
                        axes2d[k, i].set_ylim(ylim[0], ylim[1])
                    #axes2d[k,i].set_ylim(1e-17,1e-15)
                    # axes2d[k,i].set_ylabel('flux $F\; [erg/s/cm^2/\AA]$', fontsize=20)
                    # axes2d[k,i].set_xlabel(r'wavelength [$\AA$]', fontsize=20)
                    # axes2d[k,i].savefig('results_fit_sed_mat/day_'+str(round(j['time'],3))+'/spectrum.png',
                    #              facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format='png', transparent=False,
                    #              bbox_inches=None, pad_inches=0.5)
                    # axes2d[k,i].show()
            #ax = pylab.gca()
            #box = ax.get_position()
            #ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
            #axes2d[1, 2 ].legend(loc='lower right',bbox_to_anchor=(1.2, 0.5))
            #art=[]

            from matplotlib.lines import Line2D
            #handles, labels = plt.gca().get_legend_handles_labels()
            #print(handles)
            #fig.legend(handles, labels, loc='lower left')
            #print(bands_in_legend)
            #print(order_in_legend)
            #pdb.set_trace()

            labels=[name_bands[i] for i in np.array(bands_in_legend)[np.argsort(order_in_legend)]]
            lines = [Line2D([0], [0], color=color[i],marker=symbol[i] ) for i in np.array(bands_in_legend)[np.argsort(order_in_legend)]]
            fig.legend(lines,labels,loc='lower left')
            #axes2d[1, 1].legend(loc='lower center', bbox_to_anchor = (0.5, -0.3),ncol=8)
            #art.append(lgd)
            axes2d[3,1].set_xlabel(r'wavelength [$\rm \AA$]', fontsize=20)
            axes2d[1,0].set_ylabel(r'flux $\rm F\; [erg/s/cm^2/\AA]$', fontsize=20)
            #plt.tight_layout()
            fig.subplots_adjust(left=0.2)

            pylab.savefig(output + '/2D_SEDs_16.png',
                              facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format='png', transparent=False,
                              bbox_inches=None, pad_inches=1.5)
            plt.show()
        else:
            #print('number is',number)
            fig, axes2d = plt.subplots(nrows=number, ncols=1,
                                       sharex=True, sharey=True,
                                       figsize=(10, 10))
            #print('the shape of axes2d is '.format(axes2d))
            #pdb.set_trace()
            Spectra2D = np.empty((number, 1), dtype=object)
            #print(np.shape(Spectra2D))
            #print(np.shape(Spectra))
            for i, j in enumerate(Spectra2D):
                Spectra2D[i] = Spectra[i]

            # Result_T2D = np.empty((number, 1), dtype=object)
            # for i,j in enumerate(Spectra2)
            Result_T2D = Best[:, 1]
            Result_R2D = Best[:, 4]
            bands_in_legend = []
            #line = np.empty((3, 3))
            for i, row in enumerate(Spectra2D):
                #print('ble')
                #pdb.set_trace()
                for k, j in enumerate(row):
                    #print(i)
                    # print(j['time'])
                    # print(j.keys())
                    best_fit_full = black_body_flux_density.black_body_flux_density(Result_T2D[i],
                                                                                    np.arange(1e-7, 3e-6, 1e-9),
                                                                                    'P',
                                                                                    distance_pc=distance_pc,
                                                                                    Radius=distances_conversions.cm_to_solar_radius(
                                                                                        Result_R2D[i]),
                                                                                    redshift=redshift, Ebv=EBV)[2]
                    axes2d[i].plot(best_fit_full[:, 0] * 1e10, best_fit_full[:, 1], 'k-')

                    #axes2d[i].plot(np.arange(1e-7, 3e-6, 1e-9) * 1e10,
                    #                  black_body_flux_density.black_body_flux_density(Result_T2D[i],
                     #                                                                 np.arange(1e-7, 3e-6, 1e-9),
                     #                                                                 'P',
                     #                                                                 distance_pc=distance_pc,
                     #                                                                 Radius=distances_conversions.cm_to_solar_radius(
                     #                                                                     Result_R2D[i]),
                     #                                                                 redshift=redshift, Ebv=EBV)[
                     #                     2][:, 1], 'k-')
                    if 'UVW1' in j.keys():
                        axes2d[i].plot([wavelengths_filter['UVW1']], [j['UVW1']], color=color['UVW1'],
                                          marker=symbol['UVW1'], label='UW1')
                        axes2d[i].vlines(wavelengths_filter['UVW1'], j['UVW1_err'][0], j['UVW1_err'][1])
                        if 'UVW1' not in bands_in_legend:
                            bands_in_legend.append('UVW1')
                            # axes2d[k,i].errorbar([wavelengths_filter['UVW1']], [j['UVW1']],yerr=[j['UVW1_err']], color='purple',fmt='o', label='UVW1')
                    if 'UVW2' in j.keys():
                        axes2d[i].plot([wavelengths_filter['UVW2']], [j['UVW2']], color=color['UVW2'],
                                          marker=symbol['UVW2'], label='UW2')  # , label='UVW2')
                        axes2d[i].vlines(wavelengths_filter['UVW2'], j['UVW2_err'][0], j['UVW2_err'][1])
                        if 'UVW2' not in bands_in_legend:
                            bands_in_legend.append('UVW2')
                            # axes2d[k,i].errorbar([wavelengths_filter['UVW2']], [j['UVW2']],yerr=[j['UVW2_err']], color='k',fmt='o', label='UVW2')
                    if 'UVM2' in j.keys():
                        axes2d[i].plot([wavelengths_filter['UVM2']], [j['UVM2']], color=color['UVM2'],
                                          marker=symbol['UVM2'], label='UM2')  # , label='UVM2')
                        axes2d[i].vlines(wavelengths_filter['UVM2'], j['UVM2_err'][0], j['UVM2_err'][1])
                        if 'UVM2' not in bands_in_legend:
                            bands_in_legend.append('UVM2')
                    if 'r_p48' in j.keys():
                        axes2d[i].plot([wavelengths_filter['r_p48']], [j['r_p48']], color=color['r_p48'],
                                          marker=symbol['r_p48'], label='P48 r')  # ,label='P48 r')
                        axes2d[i].vlines(wavelengths_filter['r_p48'], j['r_p48_err'][0], j['r_p48_err'][1])
                        if 'r_p48' not in bands_in_legend:
                            bands_in_legend.append('r_p48')
                    if 'g_p48' in j.keys():
                        axes2d[i].plot([wavelengths_filter['g_p48']], [j['g_p48']], color=color['g_p48'],
                                          marker=symbol['g_p48'], label='P48 g')  # ,label='P48 g')
                        axes2d[i].vlines(wavelengths_filter['g_p48'], j['g_p48_err'][0], j['g_p48_err'][1])
                        if 'g_p48' not in bands_in_legend:
                            bands_in_legend.append('g_p48')
                    if 'i_p48' in j.keys():
                        axes2d[i].plot([wavelengths_filter['i_p48']], [j['i_p48']], color=color['i_p48'],
                                          marker=symbol['i_p48'], label='P48 i')  # ,label='P48 g')
                        axes2d[i].vlines(wavelengths_filter['i_p48'], j['i_p48_err'][0], j['i_p48_err'][1])
                        if 'i_p48' not in bands_in_legend:
                            bands_in_legend.append('i_p48')
                    if 'g_sdss' in j.keys():
                        axes2d[i].plot([wavelengths_filter['g_sdss']], [j['g_sdss']], color=color['g_sdss'],
                                          marker=symbol['g_sdss'], label='P60 g')  # ,label='P60 g')
                        axes2d[i].vlines(wavelengths_filter['g_sdss'], j['g_sdss_err'][0], j['g_sdss_err'][1])
                        if 'g_sdss' not in bands_in_legend:
                            bands_in_legend.append('g_sdss')
                    if 'i_sdss' in j.keys():
                        axes2d[i].plot([wavelengths_filter['i_sdss']], [j['i_sdss']], color=color['i_sdss'],
                                          marker=symbol['i_sdss'], label='P60 i')  # ,label='P60 i')
                        axes2d[i].vlines(wavelengths_filter['i_sdss'], j['i_sdss_err'][0], j['i_sdss_err'][1])
                        if 'i_sdss' not in bands_in_legend:
                            bands_in_legend.append('i_sdss')
                    if 'r_sdss' in j.keys():
                        axes2d[i].plot([wavelengths_filter['r_sdss']], [j['r_sdss']], color=color['r_sdss'],
                                          marker=symbol['r_sdss'], label='P60 r')  # ,label='P60 r')
                        axes2d[i].vlines(wavelengths_filter['r_sdss'], j['r_sdss_err'][0], j['r_sdss_err'][1])
                        if 'r_sdss' not in bands_in_legend:
                            bands_in_legend.append('r_sdss')
                    if 'z_sdss' in j.keys():
                        axes2d[i].plot([wavelengths_filter['z_sdss']], [j['z_sdss']], color=color['z_sdss'],
                                          marker=symbol['z_sdss'], label='P60 z')  # ,label='P60 r')
                        axes2d[i].vlines(wavelengths_filter['z_sdss'], j['z_sdss_err'][0], j['z_sdss_err'][1])
                        if 'z_sdss' not in bands_in_legend:
                            bands_in_legend.append('z_sdss')
                    if 'u_sdss' in j.keys():
                        axes2d[i].plot([wavelengths_filter['u_sdss']], [j['u_sdss']], color=color['u_sdss'],
                                          marker=symbol['u_sdss'], label='P60 u')  # ,label='P60 r')
                        axes2d[i].vlines(wavelengths_filter['u_sdss'], j['u_sdss_err'][0], j['u_sdss_err'][1])
                        if 'u_sdss' not in bands_in_legend:
                            bands_in_legend.append('u_sdss')
                    if 'u_swift' in j.keys():
                        axes2d[i].plot([wavelengths_filter['u_swift']], [j['u_swift']], color=color['u_swift'],
                                          marker=symbol['u_swift'], label='u')  # ,label='u')
                        axes2d[i].vlines(wavelengths_filter['u_swift'], j['u_swift_err'][0], j['u_swift_err'][1])
                        if 'u_swift' not in bands_in_legend:
                            bands_in_legend.append('u_swift')
                    if 'v_swift' in j.keys():
                        axes2d[i].plot([wavelengths_filter['v_swift']], [j['v_swift']], color=color['v_swift'],
                                          marker=symbol['v_swift'], label='V')  # ,label='v_swift')
                        axes2d[i].vlines(wavelengths_filter['v_swift'], j['v_swift_err'][0], j['v_swift_err'][1])
                        if 'v_swift' not in bands_in_legend:
                            bands_in_legend.append('v_swift')
                    if 'b_swift' in j.keys():
                        axes2d[i].plot([wavelengths_filter['b_swift']], [j['b_swift']], color=color['b_swift'],
                                          marker=symbol['b_swift'], label='B')  # ,label='b_swift')
                        axes2d[i].vlines(wavelengths_filter['b_swift'], j['b_swift_err'][0], j['b_swift_err'][1])
                        if 'b_swift' not in bands_in_legend:
                            bands_in_legend.append('b_swift')
                    if 'b_johnson' in j.keys():
                        axes2d[i].plot([wavelengths_filter['b_johnson']], [j['b_johnson']], color=color['b_johnson'],
                                          marker=symbol['b_johnson'], label='B')  # ,label='b_swift')
                        axes2d[i].vlines(wavelengths_filter['b_johnson'], j['b_johnson_err'][0],
                                            j['b_johnson_err'][1])
                        if 'b_johnson' not in bands_in_legend:
                            bands_in_legend.append('b_johnson')
                    if 'v_johnson' in j.keys():
                        axes2d[i].plot([wavelengths_filter['v_johnson']], [j['v_johnson']], color=color['v_johnson'],
                                          marker=symbol['v_johnson'], label='V')  # ,label='b_swift')
                        axes2d[i].vlines(wavelengths_filter['v_johnson'], j['v_johnson_err'][0],
                                            j['v_johnson_err'][1])
                        if 'v_johnson' not in bands_in_legend:
                            bands_in_legend.append('v_johnson')
                    if 'u_johnson' in j.keys():
                        axes2d[i].plot([wavelengths_filter['u_johnson']], [j['u_johnson']], color=color['u_johnson'],
                                          marker=symbol['u_johnson'], label='U')  # ,label='b_swift')
                        axes2d[i].vlines(wavelengths_filter['u_johnson'], j['u_johnson_err'][0],
                                            j['u_johnson_err'][1])
                        if 'u_johnson' not in bands_in_legend:
                            bands_in_legend.append('u_johnson')
                    if 'i_cousin' in j.keys():
                        axes2d[i].plot([wavelengths_filter['i_cousin']], [j['i_cousin']], color=color['i_cousin'],
                                          marker=symbol['i_cousin'], label='P60 r')  # ,label='P60 r')
                        axes2d[i].vlines(wavelengths_filter['i_cousin'], j['i_cousin_err'][0], j['i_cousin_err'][1])
                        if 'i_cousin' not in bands_in_legend:
                            bands_in_legend.append('i_cousin')
                    if 'r_cousin' in j.keys():
                        axes2d[i].plot([wavelengths_filter['r_cousin']], [j['r_cousin']], color=color['r_cousin'],
                                          marker=symbol['r_cousin'], label='P60 r')  # ,label='P60 r')
                        axes2d[i].vlines(wavelengths_filter['r_cousin'], j['r_cousin_err'][0], j['r_cousin_err'][1])
                        if 'r_cousin' not in bands_in_legend:
                            bands_in_legend.append('r_cousin')
                    if 'j_2mass' in j.keys():
                        axes2d[i].plot([wavelengths_filter['j_2mass']], [j['j_2mass']], color=color['j_2mass'],
                                          marker=symbol['j_2mass'], label='P60 r')  # ,label='P60 r')
                        axes2d[i].vlines(wavelengths_filter['j_2mass'], j['j_2mass_err'][0], j['j_2mass_err'][1])
                        if 'j_2mass' not in bands_in_legend:
                            bands_in_legend.append('j_2mass')
                    if 'h_2mass' in j.keys():
                        axes2d[i].plot([wavelengths_filter['h_2mass']], [j['h_2mass']], color=color['h_2mass'],
                                          marker=symbol['h_2mass'], label='P60 r')  # ,label='P60 r')
                        axes2d[i].vlines(wavelengths_filter['h_2mass'], j['h_2mass_err'][0], j['h_2mass_err'][1])
                        if 'h_2mass' not in bands_in_legend:
                            bands_in_legend.append('h_2mass')


                    axes2d[i].set_xscale("log")#, nonposx='clip')
                    axes2d[i].set_yscale("log")#, nonposx='clip')


                    # axes2d[k, i].legend(loc='top right')
                    axes2d[i].grid()
                    axes2d[i].set_title(r'JD-t$_0$={0}'.format(round(j['time'], 2)))
                    if ylim is not None:
                        print('ylim is not none!')
                        print('ylim is',ylim)
                        axes2d[i].set_ylim(ylim[0],ylim[1])
                    # axes2d[k,i].set_ylabel('flux $F\; [erg/s/cm^2/\AA]$', fontsize=20)
                    # axes2d[k,i].set_xlabel(r'wavelength [$\AA$]', fontsize=20)
                    # axes2d[k,i].savefig('results_fit_sed_mat/day_'+str(round(j['time'],3))+'/spectrum.png',
                    #              facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format='png', transparent=False,
                    #              bbox_inches=None, pad_inches=0.5)
                    # axes2d[k,i].show()
            # ax = pylab.gca()
            # box = ax.get_position()
            # ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
            # axes2d[1, 2 ].legend(loc='lower right',bbox_to_anchor=(1.2, 0.5))
            # art=[]

            from matplotlib.lines import Line2D
            # handles, labels = plt.gca().get_legend_handles_labels()
            # print(handles)
            # fig.legend(handles, labels, loc='lower left')
            labels = [name_bands[i] for i in bands_in_legend]
            lines = [Line2D([0], [0], color=color[i], marker=symbol[i]) for i in bands_in_legend]
            fig.legend(lines, labels, loc='lower left')
            # axes2d[1, 1].legend(loc='lower center', bbox_to_anchor = (0.5, -0.3),ncol=8)
            # art.append(lgd)
            axes2d[number-1].set_xlabel(r'wavelength [$\AA$]', fontsize=20)
            axes2d[1].set_ylabel('flux $F\; \rm [erg/s/cm^2/\AA]$', fontsize=20)
            # plt.tight_layout()
            fig.subplots_adjust(left=0.2)
            pylab.savefig(output + '/2D_SEDs_9.png',
                          facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format='png',
                          transparent=False,
                          bbox_inches=None, pad_inches=1.5)
            plt.show()
    else:
        print('Sorry, the number of plots can be either 9 or 16')



