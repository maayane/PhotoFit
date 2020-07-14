#! //anaconda/bin/python


import sys 
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
import dynesty
from dynesty import plotting as dyplot
from dynesty import utils as dyfunc
from . import class_chi2
from scipy.interpolate import interp2d


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
sigma_Boltzmann=5.670367e-8



def uniform_prior_transform(u,**kwargs):
	'''Description: Transforms the uniform random variable u ~ Unif[0., 1.) to uniform (M)SW priors. For more information, 
	   see dynesty documentation (https://dynesty.readthedocs.io/en/latest/quickstart.html#getting-started)
	   input: np.arrays of priors in format np.array([param_min,param_max])
	   defaults:
			T_prior=np.array([5e3,5e4])   
			R_prior=np.array([5e13,5e15])

			Tested : ?
				By : Ido Irani June 2020
			   URL :
			Example:
			Reliable: 
	'''
	inputs={'priors':[np.array([5e3,5e4]),
					  np.array([5e13,5e15])]}                            
	inputs.update(kwargs)
	priors=inputs.get('priors')

	def strech_prior(u,pmin,pmax):
		x = (pmax-pmin) * u + pmin
		return x

	prior_list = []
	for i,prior in enumerate(priors):
		pmin=prior[0]
		pmax=prior[1]
		prior_dist=strech_prior(u[i],pmin,pmax)
		prior_list.append(prior_dist)
		
	return prior_list


def plot_opt_fit_bb_flux_filter(model_nparam, bests, data,distance_pc,Ebv,R_ext,z=0, flatchain_path=None, uncertainties=None,
						 output_file_path='.', xlabel=None, ylabel=None,):
	print('*** Plots of the best fit solution ***')
	print('bests are',bests)
	best_fit = model_nparam.model_array()
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
			model_nparam.T=a
			model_nparam.r=b
			pylab.plot(model_nparam.model_array()[:, 0], model_nparam.model_array()[:, 1], 'b-',
					   alpha=0.5)
		wavelens=np.arange(1e-7, 3e-6, 5e-9)
		
		best_fit_full=black_body_flux_density.black_body_flux_density_fast(bests[0], wavelens, 'P', distance_pc=distance_pc,
														Radius=distances_conversions.cm_to_solar_radius(bests[1]),
														Ebv=Ebv,R_ext=R_ext, redshift=z)
		pylab.plot(best_fit_full[:,0]*1e10,best_fit_full[:,1], '-r',label=r'bb(bestT,bestR), redshifted and extincted')
		pylab.plot(best_fit[:, 0], best_fit[:, 1], 'ro', label=r'model(best fit param)')
		if uncertainties is not None:
			pylab.errorbar(data[:, 0], data[:, 1], yerr=0.5 * uncertainties, fmt='g.',
						   label=r'data')  # spectrum corrected
		ax = pylab.gca()
		ax.set_xlabel(xlabel)
		ax.set_ylabel(ylabel)
		pylab.legend(loc=1)
		print('T is',bests[0])
		print('R is',bests[1])
		print('Ebv is', Ebv)
		print('z is', z)
		print('distance is',distance_pc)
		pylab.savefig(output_file_path + '/best_fit_with_chain.png', facecolor='w', edgecolor='w',
					  orientation='portrait', papertype=None, format='png', transparent=False, bbox_inches=None,
					  pad_inches=0.1)


def fit_black_body_flux_filters_dynesty(Spectrum,**kwargs):
	"""
	Description: Fit a black-body spectrum alpha*F_bb to an observed spectrum, given as fuxes through a set of filters. Uses dynesty, as opposed to fit_black_body_flux_filters
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
	and change the location /Users/maayanesoumagnac/Maayane_Astro_python_library/data to the user one
	"""

	inputs = {'sampler':'dynamic',
			 'nlive':250,
			 'z':0,
			 'plot':False, 
			 #'dlogz':0.5,
			 'maxiter':100000,
			 'maxcall':500000,
			 'priors':[np.array([5e3,5e4]), np.array([5e13,5e15])],
			 'verbose':False,
			 'distance_pc':10e6,
			 'Ebv':0,
			 'output_dynesty':'./results_dynesty/',
			 'num_winners':20,
			 'already_run_dynesty':False,
			 'already_run_calc_all_chis':False,
			 'dilution_factor':10, 

			 }


	inputs.update(kwargs)
	z=inputs.get('z')
	plot=inputs.get('plot')
	priors=inputs.get('priors')
	nlive=inputs.get('nlive')
	dlogz=inputs.get('dlogz')
	maxiter=inputs.get('maxiter')
	maxcall=inputs.get('maxcall')
	num_winners=inputs.get('num_winners')
	distance_pc=inputs.get('distance_pc')
	dist=distances_conversions.pc_to_cm(distance_pc)
	verbose=inputs.get('verbose')
	output_dynesty=inputs.get('output_dynesty')
	flatchain_path=output_dynesty+'/flatchain.txt'
	outpath_weights=output_dynesty+'/'+'weights.txt'
	already_run_dynesty=inputs.get('already_run_dynesty')
	already_run_calc_all_chis=inputs.get('already_run_calc_all_chis')
	dilution_factor=inputs.get('dilution_factor')
	filters_directory=inputs.get('filters_directory')
	Ebv=inputs.get('Ebv')
	R_ext=inputs.get('R_ext')
	triangle_plot_title=inputs.get('triangle_plot_title')
	path_to_txt_file=inputs.get('path_to_txt_file')
	if os.path.exists(output_dynesty):
		logger.info('the file exists')
	else:
		os.makedirs(output_dynesty)
	#******************** convert waveunits into meters ***********************
	if np.shape(Spectrum)[1]<4:
		print('ERROR:Spectrum is in a wrong format, it is supposed to be a N-4 np array filter family, filter name, flux, error on flux')
		pdb.set_trace()
	Filter_vector=np.empty([np.shape(Spectrum)[0], 2],dtype=object)
	for i,j in enumerate(Spectrum):
		Filter_vector[i]=[str(Spectrum[i,0]),str(Spectrum[i,1])]
	wavelengths_min=np.zeros(np.shape(Filter_vector)[0])
	wavelengths_max=np.zeros(np.shape(Filter_vector)[0])
	effective_wavelength=np.zeros(np.shape(Filter_vector)[0])
	mags=np.zeros(np.shape(Filter_vector)[0])
	spectrum_flux = np.zeros((np.shape(Filter_vector)[0],3)) # centrl wavelength, flux, flux err
	lib = pyphot.get_library()
	if verbose==True:
		print('The filters library contains {0} filters,'.format(len(lib)))
		print('from which the following filters are used here:')
		print('Filter vector is',Filter_vector) 
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

	# ******************** Fit ************************ 
	class BBObservations(object):  
		def __init__(self, T, r, dat,uncert,priors): #T in K, R in cm, lam in \AA
			self.T = T
			self.r = r
			self.lam=dat[:,0]
			self._flux_grid = {}
			self._data = dat
			self._uncert = uncert
			self._priors = priors
		def myloglike(self, theta):
			self.T, self.r = theta 
			my_model = self.model_array()
			my_objective = class_chi2.objective_with_uncertainties(my_model, self._data,sigmas=self._uncert)
			chi2 = my_objective.chi_square_value()
			return -0.5 * chi2 
		def model_array(self):
				
				g = np.zeros([len(self.lam), 2])
				g[:,0]=self.lam
				priors=self._priors
				#wl=np.arange(1e-7, 3e-6, 5e-9)
				T_lims=priors[0]*np.array([0.5,1.5])
				R_lims=priors[1]*np.array([0.5,1.5])
			   #out = calc_black_body_flux_filters.calc_black_body_flux_filters(self.T,wl,Filter_vector=None, P_vector=P_vector,
			   #                                                                Radius=self.r, distance_pc=distance_pc,
			   #                                                                output_plot=False,
			   #                                                                show_plots=False, output_txt=False,
			   #                                                                lib=lib, z=z,Ebv=Ebv,R_ext=R_ext)
			   #
				

				#out = calc_black_body_flux_filters.calc_black_body_flux_filters_clean(self.T, Radius=self.r,Ebv=Ebv,distance_pc=distance_pc, P_vector=P_vector,
				#                                                            wavelengths=wl,output_plot=False,show_plots=False, output_txt=False,
				#                                                            lib=lib, z=z,R_ext=R_ext)
				wl=np.arange(1e-7, 3e-6, 1e-8)

				out,flux_grid = interp_black_body_flux_filters(self.T, Radius=self.r,Ebv=Ebv,distance_pc=distance_pc, P_vector=P_vector,
																			wavelengths=wl,output_plot=False,show_plots=False, output_txt=False,
																			lib=lib, z=z,R_ext=R_ext,f_grid=self._flux_grid,T_lims=T_lims, R_lims=R_lims)
				self._flux_grid=flux_grid	
				#print('''relative error is {}'''.format(((out[:,3]-out2[:,3])/out[:,3]<0.01).all()))
				#if not ((out[:,3]-out2[:,3])/out[:,3]<0.02).all(): 
				#    import ipdb; ipdb.set_trace()
				fluxes=out[:, 3]						    		    																	
				g[:, 1]=fluxes
				return g     
	#Important comment: Here the spectrum is not corrected. in calc_black_body_flux_filter, E and z are applied to a theoretical bb of which the synthetic photo is compared to the data
	ndim=2
	def prior_transform(u):
		x=uniform_prior_transform(u,priors=priors)
		return x
	data=spectrum_flux[:,0:2]
	uncertainties=spectrum_flux[:,2]
	#class BBObservations:
	#	def __init__(self, dat,uncert):
	#		self._flux_grid = {}
	#		self._data = dat
	#		self._uncert = uncert
	#		self._priors = priors
#
	#	def myloglike(self, theta):   
	#		param1,param2 = theta   
	#		my_model, flux_grid = model_black_body_dist_fixed(param1,param2, data[:, 0]).model_array(self._flux_grid,self._priors)
	#		self._flux_grid=flux_grid
	#		my_objective = class_chi2.objective_with_uncertainties(my_model, self._data,sigmas=self._uncert)
	#		chi2 = my_objective.chi_square_value()
	#		return -0.5 * chi2  
	bb_obs_instance = BBObservations(5e3,1e14,data,uncertainties,priors)


	if already_run_dynesty != True:
		dsampler = dynesty.DynamicNestedSampler(bb_obs_instance.myloglike, prior_transform, ndim,nlive=nlive)
		dsampler.run_nested(maxiter=maxiter, maxcall=maxcall)
		dresults = dsampler.results
		samples = dresults.samples
		weights = np.exp(dresults.logwt - dresults.logz[-1])  
		np.savetxt(flatchain_path, samples)
		np.savetxt(outpath_weights, weights)    
	else:
		samples=np.genfromtxt(flatchain_path)
		weights = np.genfromtxt(outpath_weights)    
	bests=samples[weights==np.max(weights)][0]
	rchi2=-2*bb_obs_instance.myloglike(bests)/(np.shape(data)[0] - ndim)
	#Temp_prior=priors[0]
	#Radius_prior=priors[1]
	best =  np.append(bests,rchi2)
		#fitter_general.calc_best_fit_2_param(ndim=ndim, model_nparam=bb_obs_instance,
			#								flatchain_path=flatchain_path,
			#								data=spectrum_flux[:,0:2], uncertainties=spectrum_flux[:,2], winners=num_winners,
			#								output_file_path=output_dynesty,
			#								bounds=[Temp_prior, Radius_prior], already_run_calc_all_chis=already_run_calc_all_chis,
			#								show_plots=False,dilution_factor=dilution_factor)
			#
	print('bests are',bests)
	plot_opt_fit_bb_flux_filter(model_nparam=bb_obs_instance, bests=bests,
														data=spectrum_flux[:,0:2],distance_pc=distance_pc,Ebv=Ebv,R_ext=R_ext,z=z,
														flatchain_path=output_dynesty+'/flatchain.txt',
														uncertainties=spectrum_flux[:,2],
														output_file_path=output_dynesty, xlabel='x',
														ylabel='y') 
	triangle = fitter_general.plot_2D_distributions(
		flatchain_path=output_dynesty+'/flatchain.txt', bests=bests, title=triangle_plot_title,
		output_file_path=output_dynesty, parameters_labels=['T', 'R'])
	histos = fitter_general.plot_1D_marginalized_distribution(
		flatchain_path=output_dynesty+'/flatchain.txt', bests=bests,
		output_pdf_file_path=output_dynesty, output_txt_file_path=output_dynesty,parameters_labels=['T', 'R'], number_bins=500) 
	best_temp = bests[0]
	best_radius = bests[1]  
	best_luminosity = energy_conversions.convert_energy(
		4 * math.pi * (best_radius * 1e-2) ** 2 * sigma_Boltzmann * best_temp ** 4, 'J', 'erg')  # en w 
	#******************** Plots ************************    

	result_array=np.zeros((1,6))

	result_array[0,0]=best_temp
	result_array[0,1]=best_radius
	result_array[0,2]=best_luminosity
	result_array[0,3]=(best_radius/distances_conversions.pc_to_cm(distance_pc))**2
	result_array[0, 4]=best[-1]
	result_array[0, 5]=best[-1]/(np.shape(Spectrum)[0]-2)   
	if path_to_txt_file!=None:
		np.savetxt(path_to_txt_file,result_array,header='best temperature, best radius, best luminosity, best coeff (R/d)**2, chi2,chi2/dof')
	else:
		np.savetxt(output_dynesty+'/best_fit_results.txt',np.array(result_array),header='best temperature, best radius, best luminosity,best coeff (R/d)**2')   
	best_coeff=(best_radius/distances_conversions.pc_to_cm(distance_pc))**2
	return best_temp,best_radius,best_luminosity,best_coeff,best[-1],best[-1]/(np.shape(Spectrum)[0]-2)



def interp_black_body_flux_filters(Temp,Radius,Ebv,distance_pc,grid_path='./TR_interp_grid.txt',f_grid={},Filter_vector=None,P_vector=None,R_lims=[3e13,4e15],T_lims=[1e3,6e4],wavelengths=np.arange(1e-7, 3e-6, 5e-9),output_txt=True,filters_directory=None,output_plot=True,lib=None,show_plots=False,R_ext=3.08,z=0,output_file=None,**kwargs):
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
	#if output_txt==True:
	#    if os.path.exists('./outputs_from_calc_black_body_flux_filters_function'):
	#        logger.info('output_path/txt_files did exist, I am removing it and creating a new one')
	#        shutil.rmtree('./outputs_from_calc_black_body_flux_filters_function')
	#    else:
	#        logger.info('the output file file did not exist yet. I am creating it now')
	#    os.makedirs('./outputs_from_calc_black_body_flux_filters_function')
	if lib is None:
		lib = pyphot.get_library()
	else:
		lib=lib
	
	
	if f_grid=={}:
		print('''\nConstructing T,R synthetic photometry grid in all filters for the first time in this epoch. This might take a few moments...''')
		Tv, Rv, grid=calc_black_body_flux_filters.generate_grid(T_lims,R_lims,Ebv,z,grid_path,filters_directory,distance_pc,lib,R_ext,P_vector=P_vector)
		for i in  range(len(P_vector)):
			filt=P_vector[i,0]+'+'+P_vector[i,1]
			f=interp2d(Tv, Rv, grid[filt],kind='linear',fill_value=np.nan)
			f_grid[filt]=f

	#plt.figure()
	#plt.contour(Tv,Rv,grid[filt])
	#plt.show()
	output_array = np.empty([np.shape(P_vector)[0], 4], dtype=object)
	string = np.empty(np.shape(P_vector)[0], dtype=object)
	for i in range(len(P_vector)):
		filt=P_vector[i,0]+'+'+P_vector[i,1]
		fluxes=f_grid[filt](Temp,Radius)
		output_array[i, 0] = P_vector[i, 0]
		output_array[i, 1] = P_vector[i, 1]
		output_array[i, 2] = P_vector[i, 2].cl.item()# le eff_wl de Eran (auquel je veux pouvoir comparer) est le cl de pyphot
		output_array[i, 3] = fluxes
		string[i] = output_array[i, 1] + '\n' + r'$(\lambda_c=$' + str(round(output_array[i, 2], 3)) + '$)$'        
			
	return output_array,f_grid #N*4 array where each line is: [filter family, filtername,effective wavelength of P,\bar{f}(P)], f is in [erg/sec/cm^2/\AA ]

	   
