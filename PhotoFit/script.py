
"""*******************************************************
PhotoFit: calculates the evolution in time of the effective radius
and temperature, by interpolating to a common epoch, and fitting a
blackbody (after applying E and z corrections and synthetic photometry)
to each SED.
******************************************************"""
print(__doc__)


import numpy as np
import pylab
import os
import logging
import pyphot
from pathlib import Path
import pdb
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
from . import params
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

__all__=['calculate_T_and_R_in_time']

sigma_Boltzmann = 5.670367e-8

mcmc=params.mcmc

output_file_mcmc=params.output_file_mcmc
output_file_linear=params.output_file_linear
output_file_interpolation=params.output_file_interpolation

filters_directroy=params.filters_directory

if mcmc==True:
    output=output_file_mcmc
else:
    output=output_file_linear

########## Informations on the object ##########

redshift = params.z
distance_modulus=params.distance_modulus
distance_pc=distances_conversions.DM_to_pc(distance_modulus)   # sn distance in pc
explosion_date=params.explosion_date
EBV=params.EBV

########################################### Read and process the photometric data ############################################


data_dico=read_data_from_file.read_data_into_numpy_array(params.data_file,header=True,delimiter=',',no_repeat_rows=True)[2]
# make sure you don't have cases with repeated x,y, and different yerr

lower_limit_on_flux=params.lower_limit_on_flux

data_dicts=dict()
for i,j in enumerate(data_dico['filter']):
    data_dicts[j]=dict()
    data_dicts[j]['jd']=data_dico['jd'][(data_dico['filter']==j) & (data_dico['flux']>lower_limit_on_flux)]
    data_dicts[j]['mag'] = data_dico['mag'][(data_dico['filter']==j) & (data_dico['flux']>lower_limit_on_flux)]
    data_dicts[j]['flux'] = data_dico['flux'][(data_dico['filter']==j) & (data_dico['flux']>lower_limit_on_flux)]
    data_dicts[j]['magerr'] = data_dico['magerr'][(data_dico['filter']==j) & (data_dico['flux']>lower_limit_on_flux)]
    data_dicts[j]['fluxerr'] = data_dico['fluxerr'][(data_dico['filter']==j) & (data_dico['flux']>lower_limit_on_flux)]
    data_dicts[j]['absmag'] = data_dico['absmag'][(data_dico['filter']==j) & (data_dico['flux']>lower_limit_on_flux)]
    data_dicts[j]['absmagerr'] = data_dico['absmagerr'][(data_dico['filter']==j) & (data_dico['flux']>lower_limit_on_flux)]
    data_dicts[j]['filter'] = data_dico['filter'][(data_dico['filter']==j) & (data_dico['flux']>lower_limit_on_flux)]
    data_dicts[j]['instr'] = data_dico['instr'][(data_dico['filter']==j) & (data_dico['flux']>lower_limit_on_flux)]

print(data_dicts.keys())

########################################### Make the spectrum to fit for each epoch ############################################

##### definition of interpolation dates

interpolation_dates=np.genfromtxt(params.dates_file)

print('I am interpolating on dates',interpolation_dates-explosion_date)

condition=dict()
JD_basis_interpolation=dict()
already_run_interp_errors=dict()
jd_flux_fluxerr=dict()
interp_and_errors_array=dict()
interp=dict()

already_run_interp_errors['UVW1']=params.already_run_interp_errors_UVW1
already_run_interp_errors['UVW2']=params.already_run_interp_errors_UVW2
already_run_interp_errors['UVM2']=params.already_run_interp_errors_UVM2
already_run_interp_errors['u_swift']=params.already_run_interp_errors_u_swift
already_run_interp_errors['b_swift']=params.already_run_interp_errors_b_swift
already_run_interp_errors['v_swift']=params.already_run_interp_errors_v_swift
already_run_interp_errors['r_sdss']=params.already_run_interp_errors_r_sdss
already_run_interp_errors['g_sdss']=params.already_run_interp_errors_g_sdss
already_run_interp_errors['i_sdss']=params.already_run_interp_errors_i_sdss
already_run_interp_errors['r_p48']=params.already_run_interp_errors_r_p48
already_run_interp_errors['g_p48']=params.already_run_interp_errors_g_p48


### Interpolate all data on the interpolation dates dates ####
for k in data_dicts.keys():
    print('I am looking at the {0} band of the data'.format(k))
    print('the days of {0} are {1}'.format(k,data_dicts[k]['jd']-explosion_date))
    condition[k]=np.logical_and(interpolation_dates>=np.min(data_dicts[k]['jd']),interpolation_dates<=np.max(data_dicts[k]['jd']))
    JD_basis_interpolation[k] = interpolation_dates[condition[k]]
    print('JD basis interp. of {0} over the interpolation days are {1}'.format(k,JD_basis_interpolation[k]-explosion_date))
    print('and interpolating on {0}'.format(interpolation_dates-explosion_date))
    a=np.array(list(zip(data_dicts[k]['jd'],data_dicts[k]['flux'],data_dicts[k]['fluxerr'])))
    jd_flux_fluxerr[k]=a[a[:, 0].argsort()]
    output_path=output_file_interpolation+'/errors_interpolation_results_{0}'.format(k)
    if os.path.exists(output_path):
        print(output_path+'exists')
    else:
        os.mkdir(output_path)
    interp_and_errors_array[k]=interpolate_errors.interpolate_errors\
    (jd_flux_fluxerr[k],JD_basis_interpolation[k],output_path=output_file_interpolation+'/errors_interpolation_results_{0}'.format(k),
     already_run=already_run_interp_errors[k],show_plot=False,title='{0} on the interpolation dates'.format(k))#array with days_nuv, interp1d P48, interp P48 with another method, lower limit 1sigma, upper limit 1 sigma
    interp[k]=dict()
    for i,j in enumerate(JD_basis_interpolation[k]):
        interp[k][str(round(j-explosion_date,8))]=interp_and_errors_array[k][i,:]

Spectra=[]

for i,j in enumerate(interpolation_dates):
    print('********')
    print('i is',i)
    print('********')
    epoch=dict()
    epoch['time']=round(interpolation_dates[i]-explosion_date,8)#jd_flux_fluxerr_UVW2[2,0]#
    #print('epoch[time] is',str(epoch['time']))
    for k in data_dicts.keys():
        if str(epoch['time']) in interp[k].keys():
            epoch[k]=interp[k][str(epoch['time'])][1] #flux
            epoch[k+'_err']=[interp[k][str(epoch['time'])][3],interp[k][str(epoch['time'])][4]] #fluxerr
    Spectra.append(epoch)

########################################### Filters definition ###########################################

# Filter_vector for all the filters
Filter_vector = np.empty([11, 2], dtype=object)
Filter_vector[0] = [str('Swift'), str('UVW1')]
Filter_vector[1] = [str('Swift'), str('UVW2')]
Filter_vector[2] = [str('Swift'), str('UVM2')]
Filter_vector[3] = [str('Swift'), str('u_swift')]
Filter_vector[4] = [str('Swift'), str('b_swift')]
Filter_vector[5] = [str('Swift'), str('v_swift')]
Filter_vector[6]=[str('ptf_p48'),str('r_p48')]
Filter_vector[7]=[str('ptf_p48'),str('g_p48')]
Filter_vector[8]=[str('sdss'),str('g_sdss')]
Filter_vector[9]=[str('sdss'),str('r_sdss')]
Filter_vector[10]=[str('sdss'),str('i_sdss')]

lib = pyphot.get_library()
print("Library contains: ", len(lib), " filters")

[P,wavelengths_filter]=get_filter.make_filter_object(Filter_vector,filters_directory=filters_directroy)
print('wavelengths are',wavelengths_filter)

################################## plot the x points spectra with error bars ###########################################
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
color['r_p48']='red'
color['g_p48']='green'

symbol=dict()
symbol['UVW1']='d'
symbol['UVW2']='d'
symbol['UVM2']='d'
symbol['u_swift']='d'
symbol['b_swift']='d'
symbol['v_swift']='d'
symbol['r_sdss']='*'
symbol['i_sdss']='*'
symbol['g_sdss']='*'
symbol['r_p48']='d'
symbol['g_p48']='d'


for i,j in enumerate(Spectra):#j is an epoch, epoch['UVW1']=flux, epoch['UVW1_err']=fluxerr
    print('the keys of this epoch are',j.keys())
    pylab.figure()
    for k in data_dicts.keys():
        if k in j.keys():
            print('the filter name {0} is therefore in epoch'.format(k))
            print('epoch[{0}] is {1}'.format(k,epoch[k]))
            pylab.errorbar([wavelengths_filter[k]],[j[k]],yerr=[j[k+'_err'][1]-j[k+'_err'][0]],label=k,fmt=symbol[k],color=color[k])
    pylab.legend()
    pylab.title('flux spectrum on JD={0}'.format(j['time']))
    pylab.grid()
    pylab.ylabel('flux $F\; [erg/s/cm^2/\AA]$', fontsize=20)
    pylab.xlabel(r'wavelength [$\AA$]', fontsize=20)
    pylab.tight_layout()
#pylab.show()
#pdb.set_trace()

################################################# Fit black body using mcmc ###############################

already_run_mcmc=params.already_run_mcmc
already_run_fit=[False]*len(Spectra)

already_run=True
if mcmc==True:
    fast=[False]*len(Spectra)

if already_run==False:
    Best=np.zeros((np.shape(Spectra)[0],8))
    #print(Best)
    #pdb.set_trace()
    for i,j in enumerate(Spectra): #goes through all epochs
        print(i)
        print(j)
        #pdb.set_trace()
        Spectrum = []
        if 'UVW2' in j.keys():
            Spectrum.append(np.array(['Swift', 'UVW2',j['UVW2'],j['UVW2_err'][1]-j['UVW2_err'][0]],dtype=object))
        if 'UVW1' in j.keys():
            Spectrum.append(np.array(['Swift', 'UVW1',j['UVW1'],j['UVW1_err'][1]-j['UVW1_err'][0] ],dtype=object))
        if 'UVM2' in j.keys():
            Spectrum.append(np.array(['Swift', 'UVM2',j['UVM2'],j['UVM2_err'][1]-j['UVM2_err'][0] ],dtype=object))
        if 'u_swift' in j.keys():
            Spectrum.append(np.array(['Swift', 'u_swift',j['u_swift'],j['u_swift_err'][1]-j['u_swift_err'][0] ],dtype=object))
        if 'v_swift' in j.keys():
            Spectrum.append(np.array(['Swift', 'v_swift',j['v_swift'],j['v_swift_err'][1]-j['v_swift_err'][0] ],dtype=object))
        if 'b_swift' in j.keys():
            Spectrum.append(np.array(['Swift', 'b_swift',j['b_swift'],j['b_swift_err'][1]-j['b_swift_err'][0] ],dtype=object))
        if 'r_p48' in j.keys():
            Spectrum.append(np.array(['ptf_p48', 'r_p48',j['r_p48'],j['r_p48_err'][1]-j['r_p48_err'][0]],dtype=object))
        if 'g_p48' in j.keys():
            Spectrum.append(np.array(['ptf_p48', 'g_p48',j['g_p48'],j['g_p48_err'][1]-j['g_p48_err'][0]],dtype=object))
        if 'r_sdss' in j.keys():
            Spectrum.append(np.array(['sdss', 'r_sdss',j['r_sdss'],j['r_sdss_err'][1]-j['r_sdss_err'][0]],dtype=object))
        if 'g_sdss' in j.keys():
            Spectrum.append(np.array(['sdss', 'g_sdss',j['g_sdss'],j['g_sdss_err'][1]-j['g_sdss_err'][0]],dtype=object))
        if 'i_sdss' in j.keys():
            Spectrum.append(np.array(['sdss', 'i_sdss',j['i_sdss'],j['i_sdss_err'][1]-j['i_sdss_err'][0]],dtype=object))
        '''
        if 'LCOGT_g' in j.keys():
            Spectrum.append(np.array(['sdss', 'g_sdss',j['LCOGT_g'],j['LCOGT_g_err']],dtype=object))
        if 'LCOGT_r' in j.keys():
            Spectrum.append(np.array(['sdss', 'r_sdss',j['LCOGT_r'],j['LCOGT_r_err']],dtype=object))
        if 'LCOGT_i' in j.keys():
            Spectrum.append(np.array(['sdss', 'i_sdss',j['LCOGT_i'],j['LCOGT_i_err'] ],dtype=object))
        if 'LCOGT_z' in j.keys():
            Spectrum.append(np.array(['sdss', 'z_sdss',j['LCOGT_z'],j['LCOGT_z_err'] ],dtype=object))
        if 'LCOGT_Uc' in j.keys():
            Spectrum.append(np.array(['johnson', 'u_johnson',j['LCOGT_Uc'],j['LCOGT_Uc_err'] ],dtype=object))
        if 'LCOGT_Vc' in j.keys():
            Spectrum.append(np.array(['johnson', 'v_johnson',j['LCOGT_Vc'],j['LCOGT_Vc_err'] ],dtype=object))
        if 'LCOGT_Bc' in j.keys():
            Spectrum.append(np.array(['johnson', 'b_johnson',j['LCOGT_Bc'],j['LCOGT_Bc_err']],dtype=object))
        if 'LCOGT_Ic' in j.keys():
            Spectrum.append(np.array(['cousin', 'i_cousin',j['LCOGT_Ic'] ,j['LCOGT_Ic_err']],dtype=object))
        if 'LCOGT_Rc' in j.keys():
            Spectrum.append(np.array(['cousin', 'r_cousin',j['LCOGT_Rc'],j['LCOGT_Rc_err'] ],dtype=object))
        if 'RATIR_H' in j.keys():
            Spectrum.append(np.array(['2mass', 'h_2mass',j['RATIR_H'] ,j['RATIR_H_err']],dtype=object))
        if 'RATIR_J' in j.keys():
            Spectrum.append(np.array(['2mass', 'j_2mass',j['RATIR_J'] ,j['RATIR_J_err']],dtype=object))
        if 'RATIR_Z' in j.keys():
            Spectrum.append(np.array(['sdss', 'z_sdss',j['RATIR_Z'] ,j['RATIR_Z_err']],dtype=object))
        if 'RATIR_i' in j.keys():
            Spectrum.append(np.array(['sdss', 'i_sdss',j['RATIR_i'] ,j['RATIR_i_err']],dtype=object))
        '''
        Spectrum_right_format=np.array(Spectrum)
        print('the Spectrum in the format to be fit is',Spectrum_right_format)
        if (j['time']< 1):
            hitemp = 3e5
            hirad=1e15
        elif(j['time'] < 2):
            hitemp = 5e4
            hirad = 1e15
        else:
            hitemp = 2e4
            hirad = 2e15

        math.log10(3000), math.log10(hitemp)


        if mcmc==True:
            if already_run_fit[i]==False:
                if fast[i]==False:
                    [best_temp, best_radius, best_luminosity,best_coeff,chi_square,chi_square_dof]=fit_black_body_flux_filters_mcmc.fit_black_body_flux_filters_mcmc\
                    (Spectrum_right_format,triangle_plot_title=r'$JD-t_{ref}=$'+str(round(j['time'],2)),nwalkers=100,num_steps=300,num_winners=20,
                     already_run_mcmc=already_run_mcmc,already_run_calc_all_chis=already_run_calc_all_chis,Temp_prior=np.array([3000,hitemp]),
                     Radius_prior=np.array([1e14,hirad]),initial_conditions=np.array([15000,3e14]),distance_pc=distance_pc,Ebv=EBV,ndof=None,show_plot=False,output_mcmc=output+'/day_'+str(round(j['time'],3)),show_mag_AB=True,z=z,path_to_txt_file=None,fast=fast[i],dilution_factor=10)
                else:
                    [best_temp, best_radius, best_luminosity, best_coeff] = fit_black_body_flux_filters_mcmc.fit_black_body_flux_filters_mcmc \
                        (Spectrum_right_format,
                         triangle_plot_title=r'$JD-t_{ref}=$' + str(round(j['time'] - explosion_date, 2)), nwalkers=100,
                         num_steps=350, num_winners=20, already_run_mcmc=already_run_mcmc,
                         already_run_calc_all_chis=already_run_calc_all_chis, Temp_prior=np.array([7000, 35000]),
                         Radius_prior=np.array([1e12, 1e15]), initial_conditions=np.array([15000, 8e13]),
                         distance_pc=distance_pc,
                         Ebv=EBV, ndof=None, show_plot=False,
                         output_mcmc=output+'/day_' + str(round(j['time'] - explosion_date, 3)),
                         show_mag_AB=True, z=z, path_to_txt_file=None, fast=fast[i], dilution_factor=10)

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
            if already_run_fit[i]==False:
                [Xi_array, best_temp, index_min, best_coeff, best_radius,best_luminosity]=fit_black_body_flux_filters.fit_black_body_flux_filters(
                    Spectrum_right_format,TempVec=np.logspace(math.log10(3000),math.log10(hitemp),100),num_temp_iterations=None,
                    distance=distance_pc,uncertainties=Spectrum_right_format[:,3],Ebv=EBV,z=redshift,output_file_path=output+'/day_'+str(round(j['time'],3))
                    ,path_to_txt_file=output+'/day_'+str(round(j['time'],3))+'/best_fit_result.txt',filters_directory=params.filters_directory)
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

    Best = np.genfromtxt(output+'/Results.txt', skip_header=1)

#TO DO : put outside of the code
#####################################   Plot best T and best R with error bars, compare with other methods ###############

#Best_Ofer=np.genfromtxt('./comparision_with_others/13dqy_BBtemp_rad.txt')#[4:21,:]


pylab.figure()
pylab.plot(Best[:, 0], Best[:, 1], 'ro', label='best fit')
for i, j in enumerate(Best):
    pylab.vlines(Best[i, 0] , Best[i, 2], Best[i, 3], color='red')
pylab.xlabel(r'time $[days]$')
pylab.ylabel(r'$T_{BB}\; [K]$')
pylab.grid()
pylab.legend()
# pylab.savefig('./results_fit_sed_mat/r_bb_evo.pdf',
#              facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format='pdf', transparent=False,
#              bbox_inches=None, pad_inches=0.5)
#pylab.show()

pylab.figure()
pylab.plot(Best[:, 0], Best[:, 4], 'bo', label='best fit')
for i, j in enumerate(Best):
    pylab.vlines(Best[i, 0] , Best[i, 5], Best[i, 6], color='blue')
pylab.xlabel(r'time $[days]$')
pylab.ylabel(r'$r_{BB}\; [cm]$')
pylab.grid()
pylab.legend()
# pylab.savefig('./results_fit_sed_mat/r_bb_evo.pdf',
#              facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format='pdf', transparent=False,
#              bbox_inches=None, pad_inches=0.5)

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
#pylab.savefig(output+'/L_bb_evo_log.pdf',
#              facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format='pdf', transparent=False,
#              bbox_inches=None, pad_inches=0.5)

################ L #######

if mcmc==True:
    path_to_errors=Path(output+'/results_errors_lumi')
    if path_to_errors.is_dir() == False:
        os.mkdir(output+'/results_errors_lumi')
    ####################### Calc errors on L ################

    errors_luminosity_early=np.zeros((np.shape(Spectra)[0],3))

    errors_lum_ran=False
    if errors_lum_ran==False:
        for i, j in enumerate(Spectra):
            print('i is',i)
            #print('file is output+/day_' + str(round(j['time'], 3)) + '/flatchain.txt')
            output_txt_file_path=output+'/results_errors_lumi/day_' + str(round(j['time'], 3))
            output_pdf_file_path=output_txt_file_path
            my_path=Path(output_txt_file_path)
            if my_path.is_dir()==False:
                os.mkdir(output_txt_file_path)
            temperatures_early=np.genfromtxt(output+'/day_' + str(round(j['time'], 3))+'/flatchain.txt')[:,0]
            radii_early=np.genfromtxt(output+'/day_' + str(round(j['time'], 3))+'/flatchain.txt')[:,1]
            luminosities_early=energy_conversions.convert_energy(4*math.pi*np.power(radii_early*1e-2,2)*sigma_Boltzmann*np.power(temperatures_early,4),'J','erg')
            #histos de luminosities
            histos=fitter_general.plot_1D_marginalized_distribution(luminosities_early, bests=None, output_pdf_file_path=output_pdf_file_path, output_txt_file_path=output_txt_file_path,
                                              parameters_labels=None, number_bins=1000)
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
    pylab.legend(prop={'size': 12})
    ax=pylab.gca()
    #ax.set_yscale("log")
    #pylab.grid()
    pylab.xticks(fontsize=16)
    pylab.yticks(fontsize=16)
    for i,j in enumerate(Spectra):
        pylab.vlines(j['time'],errors_luminosity_early[i,1],errors_luminosity_early[i,2],color='blue')
    pylab.grid(True, which="both")
    pylab.tight_layout()
    pylab.savefig(output+'/L_bb_evo.pdf',
                  facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format='pdf', transparent=False,
                  bbox_inches=None, pad_inches=0.5)


############## Dan's plot ##########################



fig, axes2d = plt.subplots(nrows=3, ncols=3,
                           sharex=True, sharey=True,
                           figsize=(10, 10))
Spectra2D = np.empty((3, 3), dtype=object)
Spectra2D[0:3, 0] = Spectra[0:3]
Spectra2D[0:3, 1] = Spectra[3:6]
Spectra2D[0:3, 2] = Spectra[6:9]


Result_T2D = np.empty((3, 3), dtype=object)
Result_T2D[0:3, 0] = Best[0:3, 1]
Result_T2D[0:3, 1] = Best[3:6, 1]
Result_T2D[0:3, 2] = Best[6:9, 1]


Result_R2D = np.empty((3, 3), dtype=object)
Result_R2D[0:3, 0] = Best[0:3, 4]
Result_R2D[0:3, 1] = Best[3:6, 4]
Result_R2D[0:3, 2] = Best[6:9, 4]


for i, row in enumerate(Spectra2D):
    for k, j in enumerate(row):
        # print(i)
        print(j['time'])
        print(j.keys())
        axes2d[k, i].plot(np.arange(1e-7, 3e-6, 1e-9) * 1e10,
                          black_body_flux_density.black_body_flux_density(Result_T2D[i, k], np.arange(1e-7, 3e-6, 1e-9),
                                                                          'P',
                                                                          distance_pc=distance_pc,
                                                                          Radius=distances_conversions.cm_to_solar_radius(
                                                                              Result_R2D[i, k]), redshift=z, Ebv=EBV)[
                              2][:, 1], 'k-')
        if 'UVW1' in j.keys():
            axes2d[k, i].plot([wavelengths_filter['UVW1']], [j['UVW1']], color='purple', marker='o',label='UW1')
            axes2d[k, i].vlines(wavelengths_filter['UVW1'], j['UVW1_err'][0], j['UVW1_err'][1])
            # axes2d[k,i].errorbar([wavelengths_filter['UVW1']], [j['UVW1']],yerr=[j['UVW1_err']], color='purple',fmt='o', label='UVW1')
        if 'UVW2' in j.keys():
            axes2d[k, i].plot([wavelengths_filter['UVW2']], [j['UVW2']], color='k', marker='o',label='UW2')  # , label='UVW2')
            axes2d[k, i].vlines(wavelengths_filter['UVW2'], j['UVW2_err'][0], j['UVW2_err'][1])
            # axes2d[k,i].errorbar([wavelengths_filter['UVW2']], [j['UVW2']],yerr=[j['UVW2_err']], color='k',fmt='o', label='UVW2')
        if 'UVM2' in j.keys():
            axes2d[k, i].plot([wavelengths_filter['UVM2']], [j['UVM2']], color='b', marker='o',label='UM2')  # , label='UVM2')
            axes2d[k, i].vlines(wavelengths_filter['UVM2'], j['UVM2_err'][0], j['UVM2_err'][1])

            # axes2d[k,i].errorbar([wavelengths_filter['UVM2']], [j['UVM2']],yerr=[j['UVM2_err']], color='b',fmt='o', label='UVM2')
        if 'r_p48' in j.keys():
            axes2d[k, i].plot([wavelengths_filter['r_p48']], [j['r_p48']], 'ro',label='P48 r')  # ,label='P48 r')
            axes2d[k, i].vlines(wavelengths_filter['r_p48'], j['r_p48_err'][0], j['r_p48_err'][1])
        if 'g_p48' in j.keys():
            axes2d[k, i].plot([wavelengths_filter['g_p48']], [j['g_p48']], 'go',label='P48 g')  # ,label='P48 g')
            axes2d[k, i].vlines(wavelengths_filter['g_p48'], j['g_p48_err'][0], j['g_p48_err'][1])
        if 'g_sdss' in j.keys():
            axes2d[k, i].plot([wavelengths_filter['g_sdss']], [j['g_sdss']], 'g*',label='P60 g')  # ,label='P60 g')
            axes2d[k, i].vlines(wavelengths_filter['g_sdss'], j['g_sdss_err'][0], j['g_sdss_err'][1])
        if 'i_sdss' in j.keys():
            axes2d[k, i].plot([wavelengths_filter['i_sdss']], [j['i_sdss']], 'm*',label='P60 i')  # ,label='P60 i')
            axes2d[k, i].vlines(wavelengths_filter['i_sdss'], j['i_sdss_err'][0], j['i_sdss_err'][1])
        if 'r_sdss' in j.keys():
            # print(j['SDSS_r_err'])
            axes2d[k, i].plot([wavelengths_filter['r_sdss']], [j['r_sdss']], 'r*',label='P60 r')  # ,label='P60 r')
            axes2d[k, i].vlines(wavelengths_filter['r_sdss'], j['r_sdss_err'][0], j['r_sdss_err'][1])
            #    if 'COUSIN'
        if 'u_swift' in j.keys():
            # print(j['u_err'])
            # pdb.set_trace()
            axes2d[k, i].plot([wavelengths_filter['u_swift']], [j['u_swift']], color='orange', marker='o',label='u')  # ,label='u')
            axes2d[k, i].vlines(wavelengths_filter['u_swift'], j['u_swift_err'][0], j['u_swift_err'][1])
        if 'v_swift' in j.keys():
            axes2d[k, i].plot([wavelengths_filter['v_swift']], [j['v_swift']], color='magenta', marker='o',label='V')  # ,label='v_swift')
            axes2d[k, i].vlines(wavelengths_filter['v_swift'], j['v_swift_err'][0], j['v_swift_err'][1])
        if 'b_swift' in j.keys():
            axes2d[k, i].plot([wavelengths_filter['b_swift']], [j['b_swift']], color='cyan', marker='o',label='B')  # ,label='b_swift')
            axes2d[k, i].vlines(wavelengths_filter['b_swift'], j['b_swift_err'][0], j['b_swift_err'][1])
 
        axes2d[k, i].set_xscale("log", nonposx='clip')
        axes2d[k, i].set_yscale("log", nonposx='clip')
        #axes2d[k, i].legend(loc='top right')
        axes2d[k, i].grid()
        axes2d[k, i].set_title(r'JD-t$_0$={0}'.format(round(j['time'],2)))
        axes2d[k,i].set_ylim(1e-17,1e-15)
        # axes2d[k,i].set_ylabel('flux $F\; [erg/s/cm^2/\AA]$', fontsize=20)
        # axes2d[k,i].set_xlabel(r'wavelength [$\AA$]', fontsize=20)
        # axes2d[k,i].savefig('results_fit_sed_mat/day_'+str(round(j['time'],3))+'/spectrum.pdf',
        #              facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format='pdf', transparent=False,
        #              bbox_inches=None, pad_inches=0.5)
        # axes2d[k,i].show()
    axes2d[1, 1 ].legend(loc='lower right')
    axes2d[2,1].set_xlabel(r'wavelength [$\AA$]', fontsize=20)
    axes2d[1,0].set_ylabel('flux $F\; [erg/s/cm^2/\AA]$', fontsize=20)
    pylab.savefig(output + '/2D_SEDs.pdf',
                  facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format='pdf', transparent=False,
                  bbox_inches=None, pad_inches=0.5)

plt.tight_layout()
plt.show()



