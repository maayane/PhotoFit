#! //anaconda/bin/python



import os
import numpy as np
import pdb
from . import calc_black_body_flux_filters
import pylab
from . import distances_conversions
from . import energy_conversions
import math
import logging
import pyphot
from . import fitter_general
from . import get_filter
from . import black_body_flux_density

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
sigma_Boltzmann=5.670367e-8

def fit_black_body_flux_filters_mcmc(Spectrum,nwalkers=100,num_steps=350,num_winners=20,already_run_mcmc=False,already_run_calc_all_chis=False,Temp_prior=None,Radius_prior=None,initial_conditions=None,distance_pc=None,Ebv=None,R_ext=None,ndof=None,show_plot=False,output_mcmc=None,show_mag_AB=True,z=0,triangle_plot_title=None,path_to_txt_file=None,fast=False,dilution_factor=10,filters_directory=None,verbose=False):
	"""Description: Fit a black-body spectrum alpha*F_bb to an observed spectrum, given as fuxes through a set of filters. Uses emcee, as opposed to fit_black_body_flux_filters
	that uses a grid of T and a linear fit for alpha.
	Input  :- Spectrun N-4 np array filter family, filter name, flux, error on flux
			-
			- if fast=true, the best fit will be taken s the median posterior, not as the result of an optimization algoryhtm
			- setting dilution factor is using the optimization algorythm, but on less combonations of parameters. It is an intermediate solution before moving to fast=true

	Plots and output files: -plot best fit in .outputs_from_fit_black_body_flux_filters_function/fit_result.png in flux
							- if show_mag_AB is True, same in mag AB
							-plot chi2, coeff, and R in outputs_from_fit_black_body_flux_filters_function/fit_grid_results.png
							-save sorted_Xi_array in outputs_from_fit_black_body_flux_filters_function/Sorted_fir_results.txt
	Tested : ?
	    By : Maayane T. Soumagnac Nov 2016
	   URL :
	Example:
	Reliable:
	 TO DO: give an option for the speed: decrease the length of TempVec, give the spec units as options,
	 and change the location /Users/maayanesoumagnac/Maayane_Astro_python_library/data to the user one"""

	#******************** convert waveunits into meters ***********************
	#print('I am running fit_bb_flux_filters_mcmc')
	if np.shape(Spectrum)[1]<4:
		print('ERROR:Spectrum is in a wrong format, it is supposed to be a N-4 np array filter family, filter name, flux, error on flux')
		pdb.set_trace()
	#Spectrum_good_units=np.zeros(np.shape(Spectrum))
	Filter_vector=np.empty([np.shape(Spectrum)[0], 2],dtype=object)
	for i,j in enumerate(Spectrum):
		Filter_vector[i]=[str(Spectrum[i,0]),str(Spectrum[i,1])]
	wavelengths_min=np.zeros(np.shape(Filter_vector)[0])
	wavelengths_max=np.zeros(np.shape(Filter_vector)[0])
	effective_wavelength=np.zeros(np.shape(Filter_vector)[0])
	mags=np.zeros(np.shape(Filter_vector)[0])
	spectrum_flux = np.zeros((np.shape(Filter_vector)[0],3)) # centrl wavelength, flux, flux err
	#mags_correct = np.zeros(np.shape(Filter_vector)[0])
	lib = pyphot.get_library()
	if verbose==True:
		print('The filters library contains {0} filters,'.format(len(lib)))
		print('from which the following filters are used here:')
		print('Filter vector is',Filter_vector)
	#P_vector=np.empty((np.shape(Filter_vector)[0],np.shape(Filter_vector)[1]+1),dtype=object)
	#P_vector[:,0]=Filter_vector[:,0]
	#P_vector[:,1]=Filter_vector[:,1]
	#print(P_vector)
	#pdb.set_trace()
	[P_vectorx,effective_wavelength]=get_filter.make_filter_object(Filter_vector,filters_directory=filters_directory)
	if verbose == True:
		print('the shape of Filter vector is',np.shape(Filter_vector)[0])
	P_vector=np.empty((np.shape(Filter_vector)[0],3),dtype=object)
	if verbose == True:
		print(np.shape(P_vector))
		print(P_vectorx['filter_family'])
		print(np.shape(P_vectorx['filter_family'][:]))
	P_vector[:,0]=P_vectorx['filter_family'][:]
	P_vector[:,1] = P_vectorx['filter_name'][:]
	P_vector[:,2]=P_vectorx['filter_object'][:]
	#P_vector[i, 2] = P
	'''
	for i, s in enumerate(Filter_vector):
		if s[0].lower() == 'ptf_p48':
			if s[1].lower() == 'r':
				print('You gave the R filter of the PTF_P48 family')
				Transmission = np.genfromtxt('/Users/maayanesoumagnac/Maayane_Astro_python_library/data/filters/P48_R_T.rtf', delimiter=None)
				P = pyphot.Filter(Transmission[:, 0], Transmission[:, 1], name='PTF_P48_R', dtype='photon',
								  unit='Angstrom')
				#print('the central wavelength is',P.cl)
				#print('the effective wavelength is', P.leff)
				#pdb.set_trace()
		elif s[0].lower() == 'galex':
			if s[1].lower() == 'galex_nuv':
				print('You gave the nuv filter of the galex family')
				Transmission = np.genfromtxt(
					'/Users/maayanesoumagnac/Maayane_Astro_python_library/data/filters//GALEX_NUV.dat', delimiter=None)
				P = pyphot.Filter(Transmission[:, 0], Transmission[:, 1], name='GALEX_NUV', dtype='photon',
								  unit='Angstrom')
				#print('the central wavelength is',P.cl)
				#print('the effective wavelength is', P.leff)
				#pdb.set_trace()


		else:
			f = lib.find(s[0].lower())  # filter family
			for name_filter in f:
				lib[name_filter].info(show_zeropoints=True)
			P = lib[s[1]]  # filter name
		P_vector[i,2]=P
		'''
	for i,s in enumerate(Filter_vector):
		#wavelengths_min[i]=np.min(P_vector[i,2].wavelength.magnitude)
		#wavelengths_max[i]=np.max(P_vector[i,2].wavelength.magnitude)
		#print(type(Spectrum[i,2]))
		mags[i] = -2.5 * np.log10(Spectrum[i, 2]) - P_vector[i,2].AB_zero_mag
		spectrum_flux[i,0]=P_vector[i,2].leff.item()
		spectrum_flux[i,1]=Spectrum[i,2]
		spectrum_flux[i,2]=Spectrum[i,3]
	min_wavelength=np.min(wavelengths_min)
	max_wavelength=np.max(wavelengths_max)
	wavelengths=np.linspace(min_wavelength,max_wavelength,num=1000)
	wavelengths_in_meters=wavelengths*1e-10
	#print(P_vector)
	#pdb.set_trace()
	'''
	#print wavelengths_in_meters
	#pdb.set_trace()
	#	wavelengths[i] = P.leff.item()
	#print wavelengths
	#print Spectrum_good_units
	#******************** correct flux with extinction ************************
	#if Ebv!=0:
		#convert spectrum into microns, as it is the units supported by the extinction functions
	#	Spectrum_microns = np.zeros(np.shape(Spectrum))
	#	Spectrum_corrected=np.zeros(np.shape(Spectrum))
	#	Spectrum_microns[:,0]=Spectrum_good_units[:,0]*1e6 #wavelength in microns
	#	Spectrum_microns[:,1]=Spectrum_good_units[:,1]
	#	Spectrum_corrected[:,1]=extinction.correct_obs_flux_for_extinction(Spectrum_microns,Ebv,Model=None,R=None)[:,1]
	#	Spectrum_corrected[:,0]=Spectrum_good_units[:,0]
	#	#print('SPectrum microns',Spectrum_microns
	#	#print('Spectrum corrected',Spectrum_corrected
	#	# Spectrum_corrected is now in meters and flux
	#else:
	#	#print('There is no extinction'
	#	Spectrum_corrected=Spectrum_good_units
	'''
	# ******************** Fit ************************
	if Temp_prior is None:
		if verbose == True:
			print('the prior for the temperature is [10^3,10^6]')
		Temp_prior=np.array([1e3,1e6])
	else:
		if verbose == True:
			print('the prior for the temperature is',Temp_prior)
	if Radius_prior is None:
		if verbose == True:
			print('the prior for the temperature is [1e12,1e17]')
		Radius_prior=np.array([1e12,1e17])
	else:
		if verbose == True:
			print('the prior for the temperature is',Temp_prior)
	if initial_conditions is None:
		if verbose == True:
			print('the default initial conditions are',  [1e5,1e14])
		init=[1e5,1e14]
	else:
		init=initial_conditions
	if distance_pc is None:
		dist=distances_conversions.pc_to_cm(10) # dist in cm
	else:
		dist=distances_conversions.pc_to_cm(distance_pc)

	if os.path.exists(output_mcmc):
		logger.info('the file exists')
			#shutil.rmtree('./outputs_from_fit_black_body_flux_filters_function')
	else:
			#logger.info('the output file file did not exist yet. I am creating it now')
		os.makedirs(output_mcmc)
		#output_file='./outputs_from_fit_black_body_flux_filters_function'

	#Xi_array=np.empty([np.shape(Temp)[0],4], dtype=object)
	#if Ebv is None:
	#	Ebv=0
	#Spectrum_flux_correct=np.empty(np.shape(spectrum_flux),dtype=object)
	#print(np.shape(Spectrum_flux_correct))
	##pdb.set_trace()
	#Spectrum_flux_correct[:,0] = correct_spectrum.correct_spectrum_for_redshift_and_extinction(spectrum_flux,comments="#",WaveUnits='A',Ebv=Ebv,z=z,show_plot=False,title=None,output_png_file=None)[:,0]
	##spectrum_flux[:,0] #central wavelength in AA
	#Spectrum_flux_correct[:,1] =correct_spectrum.correct_spectrum_for_redshift_and_extinction(spectrum_flux,comments="#",WaveUnits='A',Ebv=Ebv,z=z,show_plot=False,title=None,output_png_file=None)[:,1]
	#Spectrum_flux_correct[:,2]=correct_spectrum.correct_spectrum_for_redshift_and_extinction(spectrum_flux,comments="#",WaveUnits='A',Ebv=Ebv,z=z,show_plot=False,title=None,output_png_file=None)[:,2]

	#print('**********')
	#print(Spectrum_flux_correct[:, 0])
	#print(Spectrum_flux_correct[:,0:2])

	#print('***********')

	#print('I am plotting the imput spectrum of fitter_general in fit_black_body_flux_filters_mcmc ')
	#pylab.figure()
	#pylab.plot(Spectrum_flux_correct[:,0],Spectrum_flux_correct[:,1],'ro')
	#pylab.title('imput spectrum of fitter_general in fit_black_body_flux_filters_mcmc')
	#pylab.show()
	#pdb.set_trace()


	# np.empty(np.shape(spectrum_flux)[0],dtype=object)
	#print np.shape(correct_spectrum.correct_spectrum_for_redshift_and_extinction(spectrum_flux,comments="#",WaveUnits='A',Ebv=Ebv,z=z,show_plot=False,title=None,output_png_file=None))
	#pdb.set_trace()
	#print spectrum_flux
	#print correct_spectrum.correct_spectrum_for_redshift_and_extinction(spectrum_flux,comments="#",WaveUnits='A',Ebv=Ebv,z=z,show_plot=False,title=None,output_png_file=None)
	#pdb.set_trace()
	#Spectrum_error_correct[:] = correct_spectrum.correct_spectrum_for_redshift_and_extinction(spectrum_flux,comments="#",WaveUnits='A',Ebv=Ebv,z=z,show_plot=False,title=None,output_png_file=None)[:,2]
	#print Spectrum_error_correct
	#pdb.set_trace()
	class model_black_body_dist_fixed(object):  #
		def __init__(self, T, r,lam): #T in K, R in cm, lam in \AA
			self.T = T
			self.r = r
			self.lam=lam

		def model_array(self):
			g = np.zeros((len(self.lam), 2))
			g[:, 0] = self.lam
			#print(g[:,0])
			g[:, 1] = calc_black_body_flux_filters.calc_black_body_flux_filters(self.T,np.arange(1e-7, 3e-6, 1e-9),
																				Filter_vector=None, P_vector=P_vector,
																				Radius=self.r, distance_pc=distance_pc,
																				output_plot=False,
																				show_plots=False, output_txt=False,
																				lib=lib, z=z,Ebv=Ebv,R_ext=R_ext)[:, 3]
			return g


	#Important comment: Here the spectrum is not corrected. in calc_black_body_flux_filter, E and z are applied to a theoretical bb of which the synthetic photo is compared to the data

	samples = fitter_general.emcee_n_param(ndim=2, model_nparam=model_black_body_dist_fixed,
										   prior_param=[Temp_prior, Radius_prior], data=spectrum_flux[:,0:2], uncertainties=spectrum_flux[:,2],
										   initial_conditions=init,nwalkers=nwalkers,num_steps=num_steps,
										   flatchain_path=output_mcmc+'/flatchain.txt',
										   already_run=already_run_mcmc)


	if fast==False:
		best = fitter_general.calc_best_fit_n_param(ndim=2, model_nparam=model_black_body_dist_fixed,
												flatchain_path=output_mcmc+'/flatchain.txt',
												data=spectrum_flux[:,0:2], uncertainties=spectrum_flux[:,2], winners=num_winners,
												output_file_path=output_mcmc,
												bounds=[Temp_prior, Radius_prior], already_run_calc_all_chis=already_run_calc_all_chis,
												show_plots=False,dilution_factor=dilution_factor)
		bests = best[:-1]
	else:
		bests= fitter_general.calc_best_fit_n_param_fast(ndim=2,flatchain_path=output_mcmc+'/flatchain.txt')

	print('bests are',bests)

	def plot_opt_fit_bb_flux_filter(model_nparam, bests, data, flatchain_path=None, uncertainties=None,
							 output_file_path='.', xlabel=None, ylabel=None):

		print('*** Plots of the best fit solution ***')
		print('bests are',bests)

		best_fit = model_nparam(bests[0], bests[1], data[:, 0]).model_array()


		'''
		a_tcheck = calc_black_body_flux_filters.calc_black_body_flux_filters(bests[0], np.arange(1e-7, 3e-6, 1e-9),
																			 Filter_vector=None, P_vector=P_vector,
																			 Radius=bests[1], distance_pc=distance_pc,
																			 output_plot=True,
																			 show_plots=False, output_txt=False,
																			 lib=lib, z=z,Ebv=Ebv,R_ext=R_ext,output_file=output_mcmc)
		'''
		fig = pylab.figure()
		pylab.plot(best_fit[:, 0], best_fit[:, 1], 'bo', label=r'best-fit model') # black body with the wavelength of corrected spectrum through filters
		# pylab.plot(emcee_fit[:,0],emcee_fit[:,1],'b-')
		if uncertainties is not None:
			pylab.errorbar(data[:, 0], data[:, 1], yerr=0.5 * uncertainties, fmt='r.', label=r'data') #spectrum corrected
		else:
			pylab.plot(data[:, 0], data[:, 1], 'r.', label=r'data') # corrected spectrum
		# pylab.title(r'Data and best-fit model, $\chi^2/dof={0}$'.format(round(class_chi2.objective_with_uncertainties(best_fit,data,uncertainties).chi_square_value() / (np.shape(data)[0] - ndim), 3)))
		ax = pylab.gca()
		ax.set_xlabel(xlabel)
		ax.set_ylabel(ylabel)
		pylab.legend(loc=1)
		pylab.savefig(output_file_path + '/best_fit.png', facecolor='w', edgecolor='w',
					  orientation='portrait', papertype=None, format='png', transparent=False, bbox_inches=None,
					  pad_inches=0.1)

		if flatchain_path is not None:
			samples = np.genfromtxt(flatchain_path, comments='#')

			fig = pylab.figure()

			for a, b, in samples[np.random.randint(len(samples), size=15)]:
				emcee_fit_indata = model_nparam(a, b, data[:, 0])
				pylab.plot(emcee_fit_indata.model_array()[:, 0], emcee_fit_indata.model_array()[:, 1], 'b-',
						   alpha=0.5)

			best_fit_full=black_body_flux_density.black_body_flux_density(bests[0], np.arange(1e-7, 3e-6, 1e-9), 'P', distance_pc=distance_pc,
															Radius=distances_conversions.cm_to_solar_radius(bests[1]),
															Ebv=Ebv,R_ext=R_ext, redshift=z)[2]

			#print('best_fit_full is',best_fit_full)
			#print('T is', bests[0])
			#print('R is', bests[1])
			#print('distance is', distance_pc)
			#print('EBV is', Ebv)
			#pdb.set_trace()


			#best_fit_full=black_body_flux_density.black_body_flux_density(bests[0], np.arange(1e-7, 3e-6, 1e-9), 'P', distance_pc=distance_pc,
			#												Radius=bests[1],Ebv=Ebv, redshift=z)[2]  # in erg/sec/cm^2/Ang
			pylab.plot(best_fit_full[:,0]*1e10,best_fit_full[:,1], '-r',label=r'bb(bestT,bestR), redshifted and extincted')

			#black_body_flux_density.black_body_flux_density(Temp, wavelengths, 'P', distance_pc=distance_pc,
			#												Radius=Radius,
			#												Ebv=Ebv, redshift=z)[2]  # in erg/sec/cm^2/Ang


			#pylab.plot(a_tcheck[:,2],a_tcheck[:,3],'go',label=r'reconstructed filters_flux: filters(best_fit)')

			pylab.plot(best_fit[:, 0], best_fit[:, 1], 'ro', label=r'model(best fit param)')
			if uncertainties is not None:
				pylab.errorbar(data[:, 0], data[:, 1], yerr=0.5 * uncertainties, fmt='g.',
							   label=r'data')  # spectrum corrected

			'''
            if uncertainties != 'None':
                pylab.errorbar(data[:, 0], data[:, 1], yerr=0.5 * uncertainties, fmt='r.', label=r'data')
            else:
                pylab.plot(data[:, 0], data[:, 1], fmt='r.', label=r'data')
            '''
			# pylab.title(r'Data and best-fit model, $\chi^2/dof={0}$'.format(round(
			#    class_chi2.objective_with_uncertainties(best_fit, data, uncertainties).chi_square_value() / (
			#    np.shape(data)[0] - ndim), 3)))
			ax = pylab.gca()
			ax.set_xlabel(xlabel)
			ax.set_ylabel(ylabel)
			pylab.legend(loc=1)
			print('T is',bests[0])
			print('R is',bests[1])
			print('Ebv is', Ebv)
			print('z is', z)
			print('distance is',distance_pc)
			#pdb.set_trace()
			pylab.savefig(output_file_path + '/best_fit_with_chain.png', facecolor='w', edgecolor='w',
						  orientation='portrait', papertype=None, format='png', transparent=False, bbox_inches=None,
						  pad_inches=0.1)
		# pylab.show()


	plots_opt_fit = plot_opt_fit_bb_flux_filter(model_nparam=model_black_body_dist_fixed, bests=bests,
														data=spectrum_flux[:,0:2],
														flatchain_path=output_mcmc+'/flatchain.txt',
														uncertainties=spectrum_flux[:,2],
														output_file_path=output_mcmc, xlabel='x',
														ylabel='y')


	#pylab.show()


	triangle = fitter_general.plot_2D_distributions(
		flatchain_path=output_mcmc+'/flatchain.txt', bests=bests, title=triangle_plot_title,
		output_file_path=output_mcmc, parameters_labels=['T', 'R'])
	histos = fitter_general.plot_1D_marginalized_distribution(
		flatchain_path=output_mcmc+'/flatchain.txt', bests=bests,
		output_pdf_file_path=output_mcmc, output_txt_file_path=output_mcmc,parameters_labels=['T', 'R'], number_bins=500)

	#print('I got here 0')
	best_temp = bests[0]
	#print('I got here 1')
	best_radius = bests[1]
	#print('I got here 2')
	# mettre best radius en m: 4piR^2sigmaT^4
	best_luminosity = energy_conversions.convert_energy(
		4 * math.pi * (best_radius * 1e-2) ** 2 * sigma_Boltzmann * best_temp ** 4, 'J', 'erg')  # en w
	#print('I got here 0')
	#pylab.show()

	'''
	A = np.array(calc_black_body_flux_filters.calc_black_body_flux_filters(j, wavelengths_in_meters, Filter_vector,Radius=None, distance=None, output_plot=False,show_plots=False,output_txt=False, lib=lib,z=0)[:, 3])
	matrix_solution = np.dot(1. / (np.dot(A.transpose(), A)), (np.dot(A.transpose(), Spectrum_correct[:, 2])))
	coeff1[i] = matrix_solution
	Xi_array[i, 2] = fitting_tools.objective_no_cov(
			coeff1[i] * calc_black_body_flux_filters.calc_black_body_flux_filters(j, wavelengths_in_meters,
																				  Filter_vector, Radius=None,
																				  distance=None,
																				  output_txt=False,
																				  output_plot=False, lib=lib,z=0)[:, 3],Spectrum_correct[:, 2]).chi_square_value()
	Xi_array[i, 0] = j
	Xi_array[i,1]=coeff1[i]
	if coeff1[i]<=0:
		Xi_array[i,3]=0
	else:
		Xi_array[i,3]= math.sqrt(coeff1[i]*dist**2)#radius in cm

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
	'''
	#******************** Plots ************************

	'''
	#Plot the fitting results
	pylab.figure()
	pylab.plot(effective_wavelength, Spectrum[:,2], 'ro',label=r'data')
	if Ebv!=0 and z!=0:
		pylab.plot(effective_wavelength, Spectrum_correct[:, 2], 'bo', label=r'data corrected for extinction')
	if Ebv==0 and z!=0:
		pylab.plot(effective_wavelength, Spectrum_correct[:, 2], 'bo', label=r'data corrected for redshift')
	if Ebv!=0 and z==0:
		pylab.plot(effective_wavelength, Spectrum_correct[:, 2], 'bo', label=r'data corrected for extinction')
	#This is in case there are 2 param in model l1*black_body+l2: pylab.plot(Spectrum_corrected[:,0]*1e9,best_coeff1*black_body_flux_density.black_body_flux_density(best_temp, Spectrum_corrected[:, 0], 'P')[0][:,1]+best_coeff2,'-b',label='best fit')
	if z!=0:
		#print('z is ',z
		#pdb.set_trace()
		pylab.plot(effective_wavelength,
			   best_coeff1 * calc_black_body_flux_filters.calc_black_body_flux_filters(best_temp,wavelengths_in_meters,Filter_vector, Radius=None,distance=None,output_plot=False,output_txt=False,z=z)[:,3]
,
			   'go', label='best fit redshifted')
	else:
		print('central wavelength is', effective_wavelength
		pdb.set_trace()

		pylab.plot(effective_wavelength,best_coeff1 * calc_black_body_flux_filters.calc_black_body_flux_filters(best_temp,
																				wavelengths_in_meters * (z + 1),
																				Filter_vector, Radius=None,
																				distance=None, output_plot=False,
																				output_txt=False)[:, 3],'go', label = 'best fit')

	pylab.plot(np.linspace(1e-10, 1e-6, num=1000) * 1e10,black_body_flux_density.black_body_flux_density(best_temp, np.linspace(1e-10, 1e-6, num=1000), 'P',
															   distance=distance,
															   Radius=distances_conversions.cm_to_solar_radius(
																   best_radius))[2][:,1],'-g', label='best fit')
	pylab.xlabel('wavelength $[\AA]$')
	pylab.ylabel('Flux density $[erg/cm^2/s/A]$')
	pylab.title('Spectrum (flux)')
	pylab.grid()
	pylab.legend(loc=1)
	pylab.savefig(output_file+'/fit_result_FLux.png', facecolor='w', edgecolor='w',
				  orientation='portrait', papertype=None, format='png', transparent=False, bbox_inches=None,
				  pad_inches=0.1)
	#pylab.show()
	if show_mag_AB==True:
		pylab.figure()
		pylab.plot(effective_wavelength, mags[:], 'ro', label=r'data')
		if Ebv != 0:
			pylab.plot(effective_wavelength, mags_correct[:], 'bo', label=r'data corrected for extinction')
		# This is in case there are 2 param in model l1*black_body+l2: pylab.plot(Spectrum_corrected[:,0]*1e9,best_coeff1*black_body_flux_density.black_body_flux_density(best_temp, Spectrum_corrected[:, 0], 'P')[0][:,1]+best_coeff2,'-b',label='best fit')
		#pylab.plot(effective_wavelength,
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
	if TempVec == None:
		axarr[2].set_xscale('log')
	#axarr[1].xaxis.set_major_formatter(x_formatter)
	pylab.savefig(output_file+'/fit_grid_results.png', facecolor='w', edgecolor='w', orientation='portrait',
				  papertype=None, format='png', transparent=False, bbox_inches=None, pad_inches=0.5)

	print('the best fit temperature is {0}'.format(best_temp)
	if show_plot == True:
		pylab.show()
	if distance is None:
		print('the best fit radius has been calculated for a distance of 10pc and is {0} cm ({1} solar radii)'.format(best_radius,distances_conversions.cm_to_solar_radius(best_radius))
	else:
		print('the best fit radius has been calculated the provided distance of {0}pc and is {1} cm ({2} solar radii)'.format(distance,best_radius,distances_conversions.cm_to_solar_radius(best_radius))
		print('the best fit coeff is {0}'.format(best_coeff1)
	result_array=np.zeros((1,6))
	result_array[0,0]=0
	result_array[0,1]=best_temp
	result_array[0,2]=index_min
	result_array[0,3]=best_coeff1
	result_array[0,4]=best_radius
	result_array[0,5]=best_luminosity
	#print result_array
	#pdb.set_trace()
	'''
	if fast == False:
		result_array=np.zeros((1,6))
	else:
		result_array = np.zeros((1, 4))
	result_array[0,0]=best_temp
	result_array[0,1]=best_radius
	result_array[0,2]=best_luminosity
	result_array[0,3]=(best_radius/distances_conversions.pc_to_cm(distance_pc))**2
	if fast==False:
		result_array[0, 4]=best[-1]
		result_array[0, 5]=best[-1]/(np.shape(Spectrum)[0]-2)
	#print('I got here 4')
	#result_array[0,3]=0
	#result_array[0,4]=
	#result_array[0,5]=
	#print result_array
	if path_to_txt_file!=None:
		np.savetxt(path_to_txt_file,result_array,header='best temperature, best radius, best luminosity, best coeff (R/d)**2, chi2,chi2/dof')
	else:
		if fast==False:
			np.savetxt(output_mcmc+'/best_fit_results.txt',np.array(result_array),header='best temperature, best radius, best luminosity,best coeff (R/d)**2')
		else:
			np.savetxt(output_mcmc + '/best_fit_results_from_median.txt', np.array(result_array),
					   header='best temperature, best radius, best luminosity,best coeff (R/d)**2')

	#print(best_temp,best_radius,best_luminosity,best[-1],best[-1]/(np.shape(Spectrum)[0]-2),(best_radius/distances_conversions.pc_to_cm(distance_pc))**2)
	#print('I got here 4')
	if fast == False:
		return best_temp,best_radius,best_luminosity,(best_radius/distances_conversions.pc_to_cm(distance_pc))**2,best[-1],best[-1]/(np.shape(Spectrum)[0]-2)
	else:
		return best_temp, best_radius, best_luminosity, (
		best_radius / distances_conversions.pc_to_cm(distance_pc)) ** 2
