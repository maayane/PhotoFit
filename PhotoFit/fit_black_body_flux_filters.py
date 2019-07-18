#! //anaconda/bin/python



import numpy as np
import sys
import os
import numpy as np
from . import black_body_flux_density
from . import calc_black_body_flux_filters
import pylab
from . import fitting_tools
from . import distances_conversions
from . import energy_conversions
import math
import shutil
import logging
import pyphot
from . import get_filter
import pdb
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def fit_black_body_flux_filters(Spectrum,TempVec=None,num_temp_iterations=None,distance=None,Ebv=None,ndof=None,uncertainties=None,show_plot=False,output_file_path=None,show_mag_AB=True,z=0,path_to_txt_file=None,filters_directory=None):
	"""Description: Fit a black-body spectrum alpha*F_bb to an observed spectrum, given as fuxes through a set of filters.
	Input  :- Observed spectrum [Filter familly, filter_name, Flux]. flux in [m,erg/cm^/s/A]
			- Vector of temperatures to test [K].If empty matrix the use default. Default is logspace(3,6,num_temp_iterations). For each of these temperature, the code fits A*Model(T) to the data. then the best temperature is chosen
			- num_temp_iterations: optionnal number of temperatures in the defaukte TempVec. Default is 100
			- Optional distance [pc], default is 10pc.
			- Extinction correction (E_{B-V}) to apply to the spectrum before the fit. Default is 0.
			The correction is applied as: f_true=f_obs*10^(+0.4*A), where A is defined as mag_obs=mag_true+A.
			SInce the extinction is wavelength dependant, the wavelength used to correct the spectrum is the central wavelength of each filter.
			- show_plot: if set to True, show the plot. create it in anycase.
			- output_file_path: optional output file name, (if you run this function in a loop and do not want the output directories to be removed), default is outputs_from_fit_black_body_flux_filters_function
			-z: redshift to correct the spectrum for
	Output :- numpy array with the following fields:
				*temperature
				*for this temperature, best fit multiplication factor alpha. alpha is (R/d)^2, and therefore gives you the radius if you know the distance.
	 		(see my lecture notes in observationnal astronomy for an explanation of alpha).
	 			*for this temperature, value of best fit chi2
	 			*for this temperature, value of the radius
	 		-value of best temperature
	 		-value of best temperature index in Xi_array
	 		-value of best temperature coeff1
	 		-value of best radius
	 		-value of best luminosity
	Plots and output files: -plot best fit in .outputs_from_fit_black_body_flux_filters_function/fit_result.png in flux
							- if show_mag_AB is True, same in mag AB
							-plot chi2, coeff, and R in outputs_from_fit_black_body_flux_filters_function/fit_grid_results.png
							-save sorted_Xi_array in outputs_from_fit_black_body_flux_filters_function/Sorted_fir_results.txt
	Tested : ?
	    By : Maayane T. Soumagnac Nov 2016
	   URL :
	Example: [Xi_array,best_temp,index_min,best_coeff,best_radius]=fit_black_body.fit_black_body_spec(black_body_with_errors,distance=1.,show_plot=True)
	Reliable:
	 TO DO: give an option for the speed: decrease the length of TempVec, give the spec units as options,
	 and change the location /Users/maayanesoumagnac/Maayane_Astro_python_library/data to the user one"""
	#print 'the spectrum is',Spectrum
	#pdb.set_trace()
	if output_file_path is None:
		if os.path.exists('./utputs_from_fit_black_body_flux_filters_function'):
			logger.info('output_path/txt_files did exist')
			#shutil.rmtree('./outputs_from_fit_black_body_flux_filters_function')
		else:
			logger.info('the output file file did not exist yet. I am creating it now')
			os.makedirs('./outputs_from_fit_black_body_flux_filters_function')
		output_file='./outputs_from_fit_black_body_flux_filters_function'
	else:
		if os.path.exists(str(output_file_path)):
			logger.info('output_path did exist')
			#shutil.rmtree(output_file_path)
		else:
			logger.info('the output file file did not exist yet. I am creating it now')
			os.makedirs(output_file_path)
		output_file=output_file_path

	#******************** convert waveunits into meters ***********************
	Spectrum_good_units=np.zeros(np.shape(Spectrum))
	Filter_vector=np.empty([np.shape(Spectrum)[0], 2],dtype=object)
	for i,j in enumerate(Spectrum):
		Filter_vector[i]=[str(Spectrum[i,0]),str(Spectrum[i,1])]
	wavelengths_min=np.zeros(np.shape(Filter_vector)[0])
	wavelengths_max=np.zeros(np.shape(Filter_vector)[0])
	central_wavelength=np.zeros(np.shape(Filter_vector)[0])
	mags=np.zeros(np.shape(Filter_vector)[0])
	mags_correct = np.zeros(np.shape(Filter_vector)[0])
	lib = pyphot.get_library()
	#print('The filters library contains {0} filters'.format(len(lib)))
	#print('from which the following filters are used here:')
	#print('Filter_vector:',Filter_vector)
	[P_vectorx, effective_wavelength] = get_filter.make_filter_object(Filter_vector,filters_directory=filters_directory)
	#print('the shape of Filter vector is', np.shape(Filter_vector)[0])
	P_vector = np.empty((np.shape(Filter_vector)[0], 3), dtype=object)
	#print(np.shape(P_vector))
	#print(P_vectorx['filter_family'])
	#print(np.shape(P_vectorx['filter_family'][:]))
	P_vector[:, 0] = P_vectorx['filter_family'][:]
	P_vector[:, 1] = P_vectorx['filter_name'][:]
	P_vector[:, 2] = P_vectorx['filter_object'][:]
	if np.shape(Spectrum)[1]>3:
		spectrum_flux = np.zeros((np.shape(Filter_vector)[0], 3))
	else:
		spectrum_flux = np.zeros((np.shape(Filter_vector)[0], 2))
	for i,s in enumerate(Filter_vector):
		#wavelengths_min[i]=np.min(P_vector[i,2].wavelength.magnitude)
		#wavelengths_max[i]=np.max(P_vector[i,2].wavelength.magnitude)
		#print(type(Spectrum[i,2]))
		mags[i] = -2.5 * np.log10(Spectrum[i, 2]) - P_vector[i,2].AB_zero_mag
		#spectrum_flux[i,0]=P_vector[i,2].leff.item()
		spectrum_flux[i, 0] = P_vector[i, 2].cl.item()#le eff_wl de Eran (auquel je veux pouvoir comparer) est le cl de pyphot
		spectrum_flux[i,1]=Spectrum[i,2]
	if np.shape(Spectrum)[1] > 3:
		spectrum_flux[i,2]=Spectrum[i,3]
	min_wavelength=np.min(wavelengths_min)
	max_wavelength=np.max(wavelengths_max)
	wavelengths=np.linspace(min_wavelength,max_wavelength,num=1000)
	#print('wavelengths are',wavelengths)
	#pdb.set_trace()
	wavelengths_in_meters=wavelengths*1e-10
	#print(P_vector)
	'''
	for i, s in enumerate(Filter_vector):

		if s[0].lower() == 'ptf_p48':
			if s[1].lower() == 'r':
				print 'You gave the R filter of the PTF_P48 family'
				Transmission = np.genfromtxt('/Users/maayanesoumagnac/Maayane_Astro_python_library/data/filters/P48_R_T.rtf', delimiter=None)
				P = pyphot.Filter(Transmission[:, 0], Transmission[:, 1], name='PTF_P48_R', dtype='photon',
								  unit='Angstrom')
		else:
			f = lib.find(s[0].lower())  # filter family
			for name_filter in f:
				lib[name_filter].info(show_zeropoints=True)
			P = lib[s[1]]  # filter name
		wavelengths_min[i]=np.min(P.wavelength.magnitude)
		wavelengths_max[i]=np.max(P.wavelength.magnitude)
		central_wavelength[i]=P.cl.item()
		print type(Spectrum[i,2])
		mags[i] = -2.5 * np.log10(Spectrum[i, 2]) - P.AB_zero_mag

	min_wavelength=np.min(wavelengths_min)
	max_wavelength=np.max(wavelengths_max)
	wavelengths=np.linspace(min_wavelength,max_wavelength,num=1000)
	wavelengths_in_meters=wavelengths*1e-10
	'''
	# ******************** Fit ************************
	if TempVec is None:
		if num_temp_iterations is None:
			print('I am exploring 100 temperatures')
			Temp=np.logspace(3.,6.,500)
		else:
			print('I am exploring {0} temperatures'.format(num_temp_iterations))
			Temp = np.logspace(3., 6., num_temp_iterations)
	if distance is None:
		dist=distances_conversions.pc_to_cm(10) # dist in cm
	else:
		Temp=TempVec
		dist=distances_conversions.pc_to_cm(distance)
	coeff1=np.zeros(np.shape(Temp))
	Xi_array=np.empty([np.shape(Temp)[0],4], dtype=object)
	if Ebv is None:
		Ebv=0
	#Spectrum_correct=np.empty(np.shape(Spectrum),dtype=object)
	#Spectrum_correct[:,0]=Spectrum[:,0]
	#Spectrum_correct[:,1]=Spectrum[:,1]
	#Spectrum_correct[:,2]=correct_spectrum.correct_spectrum_for_redshift_and_extinction(np.array(zip(central_wavelength, Spectrum[:,2])),comments="#",WaveUnits='A',Ebv=Ebv,z=z,show_plot=False,title=None,output_png_file=None)[:,1]
	'''
	if show_mag_AB==True:
			for i, s in enumerate(Filter_vector):
				if s[0].lower() == 'ptf_p48':
					if s[1].lower() == 'r':
						print 'You gave the R filter of the PTF_P48 family'
						Transmission = np.genfromtxt(
							'/Users/maayanesoumagnac/Maayane_Astro_python_library/data/filters/P48_R_T.rtf',
							delimiter=None)
						P = pyphot.Filter(Transmission[:, 0], Transmission[:, 1], name='PTF_P48_R', dtype='photon',
										  unit='Angstrom')
				else:
					P = lib[s[1]]  # filter name
				mags_correct[i] = -2.5 * np.log10(Spectrum_correct[i, 2]) - P.AB_zero_mag
	'''
	if uncertainties is None:
		for i, j in enumerate(Temp):
			#print(wavelengths_in_meters)
			#pdb.set_trace()
			print('*********************')
			print('i=',i)
			print('*********************')
			A = np.array(calc_black_body_flux_filters.calc_black_body_flux_filters(j, np.arange(1e-7, 3e-6, 1e-9), Filter_vector=Filter_vector,Radius=None, distance_pc=None, output_plot=False,show_plots=False,output_txt=False, lib=lib,z=z,Ebv=Ebv)[:, 3])
			matrix_solution = np.dot(1. / (np.dot(A.transpose(), A)), (np.dot(A.transpose(), spectrum_flux[:,1])))
			coeff1[i] = matrix_solution
			Xi_array[i, 2] = fitting_tools.objective_no_cov(
					coeff1[i] * calc_black_body_flux_filters.calc_black_body_flux_filters(j, np.arange(1e-7, 3e-6, 1e-9),
																						  Filter_vector=Filter_vector, Radius=None,
																						  distance_pc=None,
																						  output_txt=False,
																						  output_plot=False, lib=lib,z=z, Ebv=Ebv)[:, 3],spectrum_flux[:, 0:2]).chi_square_value()
			Xi_array[i, 0] = j
			Xi_array[i,1]=coeff1[i]
			if coeff1[i]<=0:
				Xi_array[i,3]=0
			else:
				Xi_array[i,3]= math.sqrt(coeff1[i]*dist**2)#radius in cm
	else:
		invcov=np.diag(np.power(uncertainties, -2))
		#print('invcov is',invcov)
		#print('shape is',np.shape(invcov))
		#pdb.set_trace()
		for i, j in enumerate(Temp):
			index = i  #
			print('*********************')
			print('Linear fit: iteration {0}/{1}'.format(i+1,len(Temp)))
			print('*********************')
			A = np.array(
				calc_black_body_flux_filters.calc_black_body_flux_filters(j, np.arange(1e-7, 3e-6, 1e-9), Filter_vector=Filter_vector,
																		  Radius=None, distance_pc=None, output_txt=False,
																		  output_plot=False, lib=lib,z=z,Ebv=Ebv,filters_directory=filters_directory)[:, 3])
			matrix_solution=np.dot(1./(np.dot(np.dot(A.transpose(),invcov),A)),(np.dot(np.dot(A.transpose(),invcov),spectrum_flux[:,1])))
			coeff1[i] = matrix_solution
			#print('coeff1[i] is',coeff1[i])
			#print('j is',j)
			Xi_array[i, 2] = fitting_tools.objective_with_cov(
					coeff1[i] * calc_black_body_flux_filters.calc_black_body_flux_filters(j, np.arange(1e-7, 3e-6, 1e-9),
																						  Filter_vector=Filter_vector, Radius=None,
																						  distance_pc=None,
																						  output_txt=False,
																						  output_plot=False, lib=lib,z=z,Ebv=Ebv,filters_directory=filters_directory)[:,
								3], spectrum_flux[:,0:2], invcov).chi_square_value()
			coeff1[i] = matrix_solution
			Xi_array[i, 0] = j
			Xi_array[i, 1] = coeff1[i]
			if coeff1[i]<=0:
				Xi_array[i,3]=0
			else:
				Xi_array[i,3]= math.sqrt(coeff1[i]*dist**2)#radius in cm
	#print('Xi_array[:,2] is',Xi_array[:,2])
	#pdb.set_trace()
	index_min=np.argmin(Xi_array[:,2])
	best_coeff1=coeff1[index_min]
	best_temp=Temp[index_min]
	best_radius=Xi_array[index_min,3]
	#mettre best radius en m: 4piR^2sigmaT^4
	best_luminosity=energy_conversions.convert_energy(4*math.pi*(best_radius*1e-2)**2*5.7e-8*best_temp**4,'J','erg') #en w
	sorted_Xi_array = Xi_array[Xi_array[:,2].argsort()]
	if distance is None:
		np.savetxt(output_file+'/Sorted_fit_result.txt', sorted_Xi_array,
				   header='best temperature, best coefficient A (in model A*blackbody), sorted chi2, best fit radius (cm) at distance 10pc ')
	else:
		np.savetxt(output_file+'/Sorted_fit_result.txt', sorted_Xi_array,
				   header='best temperature, best coefficient A (in model A*blackbody), sorted chi2, best fit radius (cm) at distance '+str(distance)+'pc')

	#******************** Plots ************************

	#Plot the fitting results
	pylab.figure()
	#print('effective wave',effective_wavelength)
	#print(Spectrum[:,2])
	print('spectrum_flux[:,0]',spectrum_flux[:,0])
	print('spectrum_flux[:,1]', spectrum_flux[:, 1])
	pylab.plot(spectrum_flux[:,0], spectrum_flux[:,1], 'mo',label=r'data')
	'''
	if Ebv!=0 and z!=0:
		pylab.plot(central_wavelength, Spectrum_correct[:, 2], 'bo', label=r'data corrected for extinction')
	if Ebv==0 and z!=0:
		pylab.plot(central_wavelength, Spectrum_correct[:, 2], 'bo', label=r'data corrected for redshift')
	if Ebv!=0 and z==0:
		pylab.plot(central_wavelength, Spectrum_correct[:, 2], 'bo', label=r'data corrected for extinction')
	'''
	#This is in case there are 2 param in model l1*black_body+l2: pylab.plot(Spectrum_corrected[:,0]*1e9,best_coeff1*black_body_flux_density.black_body_flux_density(best_temp, Spectrum_corrected[:, 0], 'P')[0][:,1]+best_coeff2,'-b',label='best fit')
	#print(best_coeff1*calc_black_body_flux_filters.calc_black_body_flux_filters(best_temp,np.arange(1e-7, 3e-6, 1e-9),Filter_vector=Filter_vector,filters_directory=filters_directory,Radius=None,distance_pc=None,output_plot=False,output_txt=False,z=z,Ebv=Ebv)[:,3])
	#pdb.set_trace()
	#print(Filter_vector)
	#pdb.set_trace()
	if z!=0:
		pylab.plot(spectrum_flux[:,0],
			   best_coeff1 * calc_black_body_flux_filters.calc_black_body_flux_filters(best_temp,np.arange(1e-7, 3e-6, 1e-9),Filter_vector=Filter_vector,filters_directory=filters_directory,Radius=None,distance_pc=None,output_plot=False,output_txt=False,z=z,Ebv=Ebv)[:,3]
,
			   'ro', label='synphot(best fit redshifted+extincted)')
	#else:
	#	print 'central wavelength is', central_wavelength
		#pdb.set_trace()
	#	pylab.plot(central_wavelength,best_coeff1 * calc_black_body_flux_filters.calc_black_body_flux_filters(best_temp,
	#																			wavelengths_in_meters * (z + 1),
	#																			Filter_vector, Radius=None,
	#																			distance_pc=None, output_plot=False,
	#																			output_txt=False)[:, 3],'go', label = 'best fit')
	#print('bb')
	bb=black_body_flux_density.black_body_flux_density(best_temp, np.arange(1e-7, 3e-6, 1e-9), 'P',
															   distance_pc=distance,
															   Radius=distances_conversions.cm_to_solar_radius(
																   best_radius),redshift=z,Ebv=Ebv)[2]



	#print('after bb')
	#pylab.figure()
	pylab.plot(bb[:,0]*1e10,bb[:,1],'-g', label='best fit redshifted+extincted')
	pylab.xlabel('wavelength $[\AA]$')
	pylab.ylabel('Flux density $[erg/cm^2/s/A]$')
	pylab.title('Spectrum (flux)')
	pylab.grid()
	pylab.legend(loc=1)
	pylab.savefig(output_file+'/fit_result_FLux.png', facecolor='w', edgecolor='w',
				  orientation='portrait', papertype=None, format='png', transparent=False, bbox_inches=None,
				  pad_inches=0.1)
	#pylab.show()
	'''
	if show_mag_AB==True:
		pylab.figure()
		pylab.plot(central_wavelength, mags[:], 'ro', label=r'data')
		if Ebv != 0:
			pylab.plot(central_wavelength, mags_correct[:], 'bo', label=r'data corrected for extinction')
		# This is in case there are 2 param in model l1*black_body+l2: pylab.plot(Spectrum_corrected[:,0]*1e9,best_coeff1*black_body_flux_density.black_body_flux_density(best_temp, Spectrum_corrected[:, 0], 'P')[0][:,1]+best_coeff2,'-b',label='best fit')
		#pylab.plot(central_wavelength,
		#		   best_coeff1 * calc_black_body_flux_filters.calc_black_body_flux_filters(best_temp,
		#																				   wavelengths_in_meters,
		#																				   Filter_vector, Radius=None,
		#																				   distance=None,
		#																				   output_plot=False,
		#																				   output_txt=False)[:, 3]
		#		   ,
		#		   'go', label='best fit')
		pylab.xlabel('wavelength $[\AA]$')
		pylab.ylabel('mag_AB')
		pylab.grid()
		pylab.title('Spectrum (mag)')
		pylab.legend(loc=1)
		pylab.savefig(output_file + '/fit_result_Mag.png', facecolor='w', edgecolor='w',
					  orientation='portrait', papertype=None, format='png', transparent=False, bbox_inches=None,
					  pad_inches=0.1)

	'''
	#Plot Xi_array fields
	f, axarr = pylab.subplots(3, sharex=True)
	#x_formatter = ScalarFormatter(useOffset=False)
	axarr[0].plot(Xi_array[:,0], Xi_array[:,2], 'b',
				  label=r'$\chi^2$')
	axarr[0].plot(best_temp, Xi_array[index_min, 2], 'ro',markersize=5)
	axarr[0].set_title(r'Results of fit')
	axarr[0].set_ylabel(r'best fit $\chi^2$')
	#axarr[1].set_xlabel(r'checked temperature')
	axarr[0].grid()
	axarr[0].axvline(x=best_temp,color='k',linestyle='--')
	if TempVec is None:
		axarr[0].set_xscale('log')
	axarr[1].plot(Xi_array[:,0], Xi_array[:,1],'r',
				  label=r'multiplicative factor')
	axarr[1].set_ylabel(r'best fit multiplicative factore')
	axarr[1].plot(best_temp, Xi_array[index_min, 1], 'ro',markersize=5)
	#pylab.xlim(min(julian_days_NUV), max(julian_days_NUV))
	#pylab.legend(label=['common measurement days'])
	# pylab.gca().invert_yaxis()
	axarr[1].grid()
	axarr[1].axvline(x=best_temp,color='k',linestyle='--')
	if TempVec is None:
		axarr[1].set_xscale('log')
	#axarr[0].xaxis.set_major_formatter(x_formatter)
	axarr[2].plot(Xi_array[:,0], Xi_array[:,3], 'm',
				  label=r'Radius')
	axarr[2].plot(best_temp, Xi_array[index_min, 3], 'ro',markersize=5)
	axarr[2].set_ylabel(r'best fit radius')
	axarr[2].set_xlabel(r'checked temperature')
	axarr[2].grid()
	axarr[2].axvline(x=best_temp,color='k',linestyle='--')
	if TempVec is None:
		axarr[2].set_xscale('log')
	#axarr[1].xaxis.set_major_formatter(x_formatter)
	pylab.savefig(output_file+'/fit_grid_results.png', facecolor='w', edgecolor='w', orientation='portrait',
				  papertype=None, format='png', transparent=False, bbox_inches=None, pad_inches=0.5)

	print('the best fit temperature is {0}'.format(best_temp))
	if show_plot == True:
		pylab.show()
	if distance is None:
		print('the best fit radius has been calculated for a distance of 10pc and is {0} cm ({1} solar radii)'.format(best_radius,distances_conversions.cm_to_solar_radius(best_radius)))
	else:
		print('the best fit radius has been calculated the provided distance of {0}pc and is {1} cm ({2} solar radii)'.format(distance,best_radius,distances_conversions.cm_to_solar_radius(best_radius)))
		print('the best fit coeff is {0}'.format(best_coeff1))
	result_array=np.zeros((1,6))
	result_array[0,0]=0
	result_array[0,1]=best_temp
	result_array[0,2]=index_min
	result_array[0,3]=best_coeff1
	result_array[0,4]=best_radius
	result_array[0,5]=best_luminosity
	#print result_array
	#pdb.set_trace()
	if path_to_txt_file  is not None:
		np.savetxt(path_to_txt_file,result_array,header='place_holder, best temperature, index, best coeff, best radius, best luminosity')
	else:
		np.savetxt(output_file+'/best_fit_results.txt',np.array(result_array),header='place_holder,best temperature, index, best coeff, best radius, best luminosity')
	return Xi_array,best_temp,index_min, best_coeff1, best_radius,best_luminosity
